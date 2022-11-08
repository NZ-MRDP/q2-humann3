import qiime2.plugin
from q2_types.feature_table import FeatureTable, Frequency, RelativeFrequency
from q2_types.per_sample_sequences import SequencesWithQuality
from q2_types.sample_data import SampleData
from qiime2.plugin import Bool, Choices, Float, Int, Range, SemanticType, Str

import q2_humann3
from q2_humann3._format import (Bowtie2IndexDirFmt2, HumannDbDirFormat,
                                HumannDbFileFormat,
                                HumannDBSingleFileDirFormat)
from q2_humann3._types import (HumannDB, Nucleotide, Pathway, PathwayMapping,
                               Protein, ReferenceNameMapping)

plugin = qiime2.plugin.Plugin(
    name="humann3",
    version="0.0.0",
    website="https://huttenhower.sph.harvard.edu/humann/",
    package="q2_humann3",
    user_support_text=(
        "To get help with HUMAnN3, please post a question to "
        "the HUMAnN Google Group form: "
        "https://groups.google.com/forum/#!forum/humann-users"
    ),
    citation_text=None,
)

plugin.register_semantic_types(
    HumannDB, Nucleotide, Pathway, Protein, ReferenceNameMapping
)

plugin.register_formats(
    HumannDbDirFormat,
    HumannDbFileFormat,
    HumannDBSingleFileDirFormat,
    Bowtie2IndexDirFmt2,
)

plugin.register_semantic_type_to_format(
    HumannDB[Nucleotide | Protein], HumannDbDirFormat
)
Bowtie2Index2 = SemanticType("Bowtie2Index2")

plugin.register_semantic_types(Bowtie2Index2)
plugin.register_semantic_type_to_format(Bowtie2Index2, Bowtie2IndexDirFmt2)

# TODO: Add pathways and investigate what the other "databases" look like
plugin.register_semantic_type_to_format(
    HumannDB[
        PathwayMapping | Pathway | ReferenceNameMapping,
    ],
    HumannDBSingleFileDirFormat,
)
plugin.methods.register_function(
    function=q2_humann3.run,
    inputs={
        "demultiplexed_seqs": SampleData[SequencesWithQuality],
        "nucleotide_database": HumannDB[Nucleotide],
        "protein_database": HumannDB[Protein],
        "pathway_database": HumannDB[Pathway],
        "pathway_mapping": HumannDB[PathwayMapping],
        "bowtie_database": Bowtie2Index2,
    },
    parameters={
        "threads": Int,
        "memory_use": Str % Choices({"minimum", "maximum"}),
        "metaphlan_stat_q": Float % Range(0, 1, inclusive_end=True),
    },
    outputs=[
        ("genefamilies", FeatureTable[Frequency]),  # type: ignore
        ("pathcoverage", FeatureTable[RelativeFrequency]),  # type: ignore
        ("pathabundance", FeatureTable[RelativeFrequency]),  # type: ignore
    ],
    input_descriptions={
        "demultiplexed_seqs": (
            "sequence files that you wish to profile,"
            " in fastq (or fastq.gz) format. Multiple"
            " sequence files per sample need to first be"
            " concatenated into 1 file."
        ),
        "nucleotide_database": "directory containing the nucleotide database",
        "protein_database": "directory containing the protein database",
        "pathway_database": "directory providing a tab-delimited mapping",
        "pathway_mapping": "directory providing the pathways mapping",
        "bowtie_database": "directory containing the bowtie2 executable",
    },
    parameter_descriptions={
        "threads": "number of threads/processes",
        "memory_use": "the amount of memory to use",
        "metaphlan_stat_q": "Quantile value for the robust average",
    },
    output_descriptions={
        "genefamilies": (
            "This file details the abundance of each gene family" " in the community."
        ),
        "pathcoverage": (
            "Pathway coverage provides an alternative description"
            " of the presence (1) and absence (0) of pathways in"
            " a community, independent of their quantitative"
            " abundance."
        ),
        "pathabundance": (
            "This file details the abundance of each pathway in"
            " the community as a function of the abundances of"
            " the pathway's component reactions, with each"
            " reaction's abundance computed as the sum over"
            " abundances of genes catalyzing the reaction."
        ),
    },
    name="Characterize samples using HUMAnN3",
    description="Execute the HUMAnN3",
)

plugin.methods.register_function(
    function=q2_humann3.rename_table,
    inputs={
        "table": FeatureTable[Frequency | RelativeFrequency],
        "reference_mapping": HumannDB[ReferenceNameMapping],
    },
    parameters={
        "name": Int % None
        | Str  # type: ignore
        % Choices(
            {
                "kegg-orthology",
                "kegg-pathway",
                "kegg-module",
                "ec",
                "metacyc-rxn",
                "metacyc-pwy",
                "pfam",
                "eggnog",
                "go",
                "infogo1000",
            }
        ),
        "simplify": Bool,
    },
    name="Rename Table",
    outputs=[
        ("rename_table", FeatureTable[Frequency | RelativeFrequency]),  # type: ignore
    ],
    description="Rename the feature table IDs",
    input_descriptions={
        "table": (
            "Utility for renormalizing TSV files Each level of a stratified."
            " Table will be normalized using the desired scheme."
        ),
        "reference_mapping": (
            "Pass an explicit database to use for renaming."
            " Use if name option is not available."
        ),
    },
    parameter_descriptions={
        "name": "Name of the reference database to use for renaming files",
        "simplify": "Remove non-alphanumeric characters from names",
    },
    output_descriptions={
        "rename_table": "The modified output table",
    },
)
