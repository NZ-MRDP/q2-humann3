import qiime2.plugin
from q2_types.feature_table import FeatureTable, Frequency, RelativeFrequency
from q2_types.per_sample_sequences import SequencesWithQuality
from q2_types.sample_data import SampleData
from qiime2.plugin import Bool, Choices, Float, Int, Range, SemanticType, Str

import q2_humann3
from q2_humann3._format import (Bowtie2IndexDirFmt2, HumannDbDirFormat,
                                HumannDbFileFormat,
                                HumannDBSingleFileDirFormat,
                                HumannDBSingleReferenceFileDirFormat)
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
    HumannDBSingleReferenceFileDirFormat,
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
        PathwayMapping | Pathway,
    ],
    HumannDBSingleFileDirFormat,
)

plugin.register_semantic_type_to_format(
    HumannDB[ReferenceNameMapping],
    HumannDBSingleReferenceFileDirFormat,
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
        "n_parallel_samples": Int,
        "humann3_threads": Int,
        "memory_use": Str % Choices({"minimum", "maximum"}),
        "metaphlan_stat_q": Float % Range(0, 1, inclusive_end=True),
    },
    outputs=[
        ("genefamilies", FeatureTable[Frequency]),  # type: ignore
        ("pathcoverage", FeatureTable[Frequency]),  # type: ignore
        ("pathabundance", FeatureTable[RelativeFrequency]),  # type: ignore
        ("taxonomy", FeatureTable[RelativeFrequency]),  # type: ignore
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
        "bowtie_database": "directory containing the bowtie2 database reference files",
    },
    parameter_descriptions={
        "n_parallel_samples": (
            "Humann3 runs explicitly on a per-sample basis however q2-humann3 runs on a table of samples, i.e. "
            "multiple samples. The sample-thread specifies how many samples should be run in parallel, this should "
            "be a maximuim of n-1 processors. It is important to note that the memory required will scale "
            "with threads. If 8GB of ram is expected for a single sample and 4 threads are selected you will need a "
            "minimum of 32GB of ram. Memory use is highly dataset dependent and will change based on reference "
            "databases, however, as a general run Humann3 will consume ~16GB of ram per sample. Use with caution"
        ),
        "humann3_threads": (
            " The number of threads humann3 will use when processing a single sample. This will"
            " for example call metaphlan with the number of threads specified where metaphlan will"
            " implement its own multithreading. This should not be used in conjunction with the"
            " sample_threads parameter. One or both of these paramers should be 1 unless you"
            " are absolutely certain you have enough processors available. You should expect"
            "  the number of processors required to be the number of"
            " sample_threads * humann3_threads + 1, and the required memory to be a large multiple"
            " of that, though it is highly dependent on your data set"
        ),
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
        "taxonomy": (
            "Taxonomic profile of microbial community of samples,"
            " generated using clade-specific marker genes."
        ),
    },
    name="Characterize samples using HUMAnN3",
    description="Execute the HUMAnN3",
)

_rename_params = {
    "parameters": {
        "name": Str  # type: ignore
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
                "uniref90",
                "uniref50"

            }
        ),
        "simplify": Bool,
    },
    "description": "Rename the feature table IDs",
    "input_descriptions": {
        "table": (
            "Utility for renormalizing TSV files Each level of a stratified."
            " Table will be normalized using the desired scheme."
        ),
        "reference_mapping": (
            "Pass an explicit database to use for renaming."
            " Use if name option is not available."
        ),
    },
    "parameter_descriptions": {
        "name": "Name of the reference database to use for renaming files",
        "simplify": "Remove non-alphanumeric characters from names",
    },
    "output_descriptions": {
        "rename_table": "The modified output table",
    },
}

# qiime only allows one type output per function
# Because pathways and gene families have different types
# we have to a a function for each
plugin.methods.register_function(
    function=q2_humann3.rename_pathways,
    inputs={
        "table": FeatureTable[RelativeFrequency],
        "reference_mapping": HumannDB[ReferenceNameMapping],
    },
    name="Rename Pathways Table",
    outputs=[
        ("rename_table", FeatureTable[RelativeFrequency]),  # type: ignore
    ],
    **_rename_params
)

plugin.methods.register_function(
    function=q2_humann3.rename_gene_families,
    inputs={
        "table": FeatureTable[Frequency],
        "reference_mapping": HumannDB[ReferenceNameMapping],
    },
    name="Rename Pathways Table",
    outputs=[
        ("rename_table", FeatureTable[RelativeFrequency]),  # type: ignore
    ],
    **_rename_params
)
