import os
import subprocess
import tempfile

import biom
# from q2_types.feature_table import FeatureTable, Frequency
from q2_types.per_sample_sequences import (
    FastqGzFormat, SingleLanePerSampleSingleEndFastqDirFmt)

from q2_humann3._format import (Bowtie2IndexDirFmt2, HumannDbDirFormat,
                                HumannDBSingleFileDirFormat)

# from q2_types.bowtie2 import Bowtie2IndexDirFmt

# import typing


def _single_sample(
    sequence_sample_path: str,
    nucleotide_database_path: str,
    protein_database_path: str,
    pathway_database_path: str,
    pathway_mapping_path: str,
    bowtie_database_path: str,
    threads: int,
    memory_use: str,
    metaphlan_stat_q: float,
    output: str,
) -> None:
    print(
        "%    %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   "
    )
    print(os.listdir(bowtie_database_path))
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
        # TODO: Do we still need this flag if we're breaking up the
        #       arguments?
        # "--metaphlan-options",
        "--stat-q {} --add-viruses --unclassified-estimation".format(
            metaphlan_stat_q
        ),
        # --offline # Don't check for or install databases
        "--offline --bowtie2db {} --index mpa_vJan21_CHOCOPhlAnSGB_202103".format(
            bowtie_database_path
        ),
    ]
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
    Run samples through humann2.

    Parameters
    ----------
    samples : SingleLanePerSampleSingleEndFastqDirFmt
        Samples to process
    threads : int, optional
        The number of threads that humann2 should use
    memory_use : str, optional
        The amount of memory to use, default is minimum
    metaphlan_stat_q : float, optional
        Quantile value for the robust average, for Metaphlan, default is 0.2

    Notes
    -----
    This command consumes per-sample FASTQs, and takes those data through
    "humann2", then through "humann2_join_tables" and finalizes with
    "humann2_renorm_table".

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
            _single_sample(
                str(view),
                nucleotide_database_path=str(nucleotide_database),
                protein_database_path=str(protein_database),
                pathway_database_path=str(pathway_database),
                pathway_mapping_path=str(pathway_mapping),
                bowtie_database_path=str(bowtie_database),
                threads=threads,
                memory_use=memory_use,
                metaphlan_stat_q=metaphlan_stat_q,
                output=tmp,
            )

        final_tables = {}
        for (name, method) in [
            ("genefamilies", "cpm"),
            ("pathcoverage", "relab"),
            ("pathabundance", "relab"),
        ]:

            joined_path = os.path.join(tmp, "%s.biom" % name)
            result_path = os.path.join(tmp, "%s.%s.biom" % (name, method))

            _join_tables(tmp, joined_path, name)
            _renorm(joined_path, method, result_path)

            final_tables[name] = biom.load_table(result_path)

    return (
        final_tables["genefamilies"],
        final_tables["pathcoverage"],
        final_tables["pathabundance"],
    )
