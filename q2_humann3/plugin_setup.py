import qiime2.plugin
from q2_types.feature_table import FeatureTable, Frequency, RelativeFrequency
from q2_types.per_sample_sequences import SequencesWithQuality
from q2_types.sample_data import SampleData
from qiime2.plugin import Choices, Int, SemanticType, Str

import q2_humann3
from q2_humann3._format import (Bowtie2IndexDirFmt2, HumannDbDirFormat,
                                HumannDbFileFormat,
                                HumannDBSingleFileDirFormat)
from q2_humann3._types import (HumannDB, Nucleotide, Pathway, PathwayMapping,
                               Protein)

# from q2_types.bowtie2 import Bowtie2Index


plugin = qiime2.plugin.Plugin(
    name="humann3",
    version="0.0.0",
    website="http://huttenhower.sph.harvard.edu/humann2",
    package="q2_humann3",
    user_support_text=(
        "To get help with HUMAnN2, please post a question to "
        "the HUMAnN Google Group form: "
        "https://groups.google.com/forum/#!forum/humann-users"
    ),
    citation_text=None,
)

plugin.register_semantic_types(HumannDB, Nucleotide, Pathway, Protein)

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
    HumannDB[PathwayMapping | Pathway], HumannDBSingleFileDirFormat
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
        "metaphlan_options": Str,
    },
    outputs=[
        ("genefamilies", FeatureTable[Frequency]),  # type: ignore
        ("pathcoverage", FeatureTable[RelativeFrequency]),  # type: ignore
        ("pathabundance", FeatureTable[RelativeFrequency]),  # type: ignore
        ("taxonomy", FeatureTable[RelativeFrequency]),  # type: ignore
    ],
    input_descriptions={
        "demultiplexed_seqs": "--- UPDATE ---",
        "nucleotide_database": "directory containing the nucleotide database",
        "protein_database": "directory containing the protein database",
        "pathway_database": "directory providing a tab-delimited mapping",
        "pathway_mapping": "directory providing the pathways mapping",
        "bowtie_database": "directory containing the bowtie2 executable",
    },
    parameter_descriptions={
        "threads": "number of threads/processes",
        "memory_use": "the amount of memory to use",
        "metaphlan_options": "options to be provided to the MetaPhlAn software",
    },
    output_descriptions={
        "genefamilies": ("This file details the abundance of each gene family"
                         " in the community."),
        "pathcoverage": ("Pathway coverage provides an alternative description"
                         " of the presence (1) and absence (0) of pathways in"
                         " a community, independent of their quantitative"
                         " abundance."),
        "pathabundance": ("This file details the abundance of each pathway in"
                          " the community as a function of the abundances of"
                          " the pathway's component reactions, with each"
                          " reaction's abundance computed as the sum over"
                          " abundances of genes catalyzing the reaction."),
        "taxonomy": ("--- UPDATE ---"),
    },
    name="Characterize samples using HUMAnN2",
    description="Execute the HUMAnN2",
)
