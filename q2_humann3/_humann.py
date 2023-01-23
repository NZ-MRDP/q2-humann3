import os
import subprocess
import tempfile
from glob import glob

import biom
from q2_types.feature_table import BIOMV210Format
from q2_types.per_sample_sequences import (
    FastqGzFormat, SingleLanePerSampleSingleEndFastqDirFmt)

from q2_humann3._format import (Bowtie2IndexDirFmt2, HumannDbDirFormat,
                                HumannDBSingleFileDirFormat,
                                HumannDBSingleReferenceFileDirFormat)


def _single_sample(
    sequence_sample_path: str,
    nucleotide_database_path: str,
    protein_database_path: str,
    pathway_database_path: str,
    pathway_mapping_path: str,
    threads: int,
    memory_use: str,
    metaphlan_options: str,
    output: str,
) -> None:
    cmd = [
        "humann3",
        "-i",
        sequence_sample_path,
        "-o",
        output,
        "--threads",
        str(threads),
        "--memory-use",
        memory_use,
        "--output-format",
        "biom",
        "--remove-column-description-output",
        "--nucleotide-database",
        nucleotide_database_path,
        "--protein-database",
        # TODO: Fix this nonsense
        protein_database_path,
        "--pathways-database",
        "{},{}".format(
            os.path.join(pathway_mapping_path, "mapping.gz"),
            os.path.join(pathway_database_path, "mapping.gz"),
        ),
        "--metaphlan-options",
        metaphlan_options,
    ]
    subprocess.run(cmd, check=True)


def _join_taxa_tables(input_dir_path: str, output_path: str):
    # I was not able to get custom output names to work from metaphlan
    # so we're looking for the default names
    taxa_tables = glob(f"{input_dir_path}/**/*metaphlan_bugs*", recursive=True)
    cmd = ["merge_metaphlan_tables.py"] + taxa_tables
    with open(output_path, "w") as outfile:
        subprocess.run(cmd, stdout=outfile)
        subprocess.run(cmd, check=True)


def _join_tables(table: str, output: str, name: str) -> None:
    """Merge multiple sample output into single tables"""
    tmp_output = output + "-actual"
    cmd = [
        "humann_join_tables",
        "-i",
        table,
        "-o",
        tmp_output,
        "--file_name",
        "%s" % name,
    ]

    subprocess.run(cmd, check=True)

    # doing convert manually as we need to filter out the leading comment as
    # humann2_renorm_table cannot handle comment lines
    for_convert = biom.load_table(tmp_output)
    lines = for_convert.to_tsv().splitlines()
    lines = lines[1:]  # drop leading comment
    with open(output, "w") as fp:
        fp.write("\n".join(lines))
        fp.write("\n")


def _renorm(table: str, method: str, output: str) -> None:
    """Renormalize a table"""
    cmd = [
        "humann_renorm_table",
        "-i",
        "%s" % table,
        "-o",
        "%s" % output,
        "-u",
        "%s" % method,
    ]
    subprocess.run(cmd, check=True)


def _metaphlan_options(bowtie2db: str, stat_q: float) -> str:
    """
    Takes the parameters needed for MetaPhlAn4 and combines them
    into a valid string.

    Parameters
    ----------
    bowtie2db : str
        directory containing the bowtie2 executable
    stat_q : float
        Quantile value for the robust average
    """
    # TODO: The index needs to be set programmatically

    # Calling bowtie in this way requires an index name, which is actually just the filenames in
    # The bowtie2db directory without the extensions
    index_name = {e.split(".")[0] for e in os.listdir(bowtie2db)}
    if len(index_name) != 1:
        raise ValueError(
            "The index files in the Bowtie database are not named in a"
            " consistent fashion. Check that all files in the bowtie"
            " database have the same base name."
        )

    (index_name,) = index_name
    # index = os.path.join(bowtie2db, index_name)
    return f"--offline --bowtie2db {bowtie2db} --index {index_name} --stat_q {stat_q} --add_viruses --unclassified_estimation"


def run(
    demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
    nucleotide_database: HumannDbDirFormat,
    protein_database: HumannDbDirFormat,
    pathway_database: HumannDBSingleFileDirFormat,
    pathway_mapping: HumannDBSingleFileDirFormat,
    bowtie_database: Bowtie2IndexDirFmt2,
    threads: int = 1,
    memory_use: str = "minimum",
    metaphlan_stat_q: float = 0.2,
) -> (biom.Table, biom.Table, biom.Table, biom.Table):  # type:  ignore
    """
    Run samples through humann3.

    Parameters
    ----------
    samples : SingleLanePerSampleSingleEndFastqDirFmt
        Samples to process
    threads : int, optional
        The number of threads that humann3 should use
    memory_use : str, optional
        The amount of memory to use, default is minimum
    metaphlan_stat_q : float, optional
        Quantile value for the robust average, for Metaphlan, default is 0.2

    Notes
    -----
    This command consumes per-sample FASTQs, and takes those data through
    "humann3", then through "humann3_join_tables" and finalizes with
    "humann3_renorm_table".

    Returns
    -------
    biom.Table
        A gene families table normalized using "cpm"
    biom.Table
        A pathway coverage table normalized by relative abundance
    biom.Table
        A pathway abundance table normalized by relative abundance
    """
    with tempfile.TemporaryDirectory() as tmp:
        iter_view = demultiplexed_seqs.sequences.iter_views(FastqGzFormat)  # type: ignore

        for _, view in iter_view:
            metaphlan_options = _metaphlan_options(
                str(bowtie_database),
                metaphlan_stat_q,
            )
            _single_sample(
                str(view),
                nucleotide_database_path=str(nucleotide_database),
                protein_database_path=str(protein_database),
                pathway_database_path=str(pathway_database),
                pathway_mapping_path=str(pathway_mapping),
                threads=threads,
                memory_use=memory_use,
                metaphlan_options=metaphlan_options,
                output=tmp,
            )

        final_tables = {}
        for (name, method) in [
            ("genefamilies", "relab"),
            ("pathcoverage", "relab"),
            ("pathabundance", "relab"),
            ("taxonomy", "relab"),
        ]:
            joined_path = os.path.join(tmp, "%s.biom" % name)
            result_path = os.path.join(tmp, "%s.%s.biom" % (name, method))

            if name == "taxonomy":
                _join_taxa_tables(input_dir_path=tmp, output_path=result_path)
                final_tables[name] = biom.load_table(result_path)
            else:
                _join_tables(table=tmp, output=joined_path, name=name)
                if name != "pathcoverage":
                    _renorm(joined_path, method, result_path)
                    final_tables[name] = biom.load_table(result_path)
                else:
                    final_tables[name] = biom.load_table(joined_path)

    return (
        final_tables["genefamilies"],
        final_tables["pathcoverage"],
        final_tables["pathabundance"],
        final_tables["taxonomy"],
    )


def rename_pathways(
    table: BIOMV210Format,
    name: str = None,
    reference_mapping: HumannDBSingleReferenceFileDirFormat = None,
    simplify: bool = False,
) -> biom.Table:  # type: ignore
    """rename_pathways.

    Parameters
    ----------
    table : BIOMV210Format
        table
    name : str
        name
    reference_mapping : HumannDBSingleFileDirFormat
        reference_mapping
    simplify : bool
        simplify

    Returns
    -------
    biom.Table

    """
    return _rename_table(table, name, reference_mapping, simplify)


def rename_gene_families(
    table: BIOMV210Format,
    name: str = None,
    reference_mapping: HumannDBSingleReferenceFileDirFormat = None,
    simplify: bool = False,
) -> biom.Table:  # type: ignore
    """rename_pathways.

    Parameters
    ----------
    table : BIOMV210Format
        table
    name : str
        name
    reference_mapping : HumannDBSingleFileDirFormat
        reference_mapping
    simplify : bool
        simplify

    Returns
    -------
    biom.Table

    """
    return _rename_table(table, name, reference_mapping, simplify)


# def rename_pathways()


def _rename_table(
    table: BIOMV210Format,
    name: str = None,
    reference_mapping: HumannDBSingleReferenceFileDirFormat = None,
    simplify: bool = False,
) -> biom.Table:  # type: ignore
    """rename_table.

    Parameters
    ----------
    table : BIOMV210Format
        table
    name : str
        name
    simplify : bool
        simplify

    Returns
    -------
    biom.Table

    """
    table_path = str(table)
    with tempfile.TemporaryDirectory() as tmp:
        output_path = os.path.join(tmp, "renorm.biom")

        cmd = [
            "humann_rename_table",
            "-i",
            table_path,
            "-o",
            output_path,
        ]
        if simplify:
            cmd.append("--simplify")

        if name:
            cmd.extend(["n", name])

        if reference_mapping:
            cmd.extend(
                ["-c", os.path.join(str(reference_mapping), "reference.txt.bz2")]
            )

        subprocess.run(cmd, check=True)

        return biom.load_table(output_path)
