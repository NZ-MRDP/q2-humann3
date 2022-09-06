import qiime2.plugin
from q2_types.feature_table import FeatureTable, Frequency, RelativeFrequency
from q2_types.per_sample_sequences import SequencesWithQuality
from q2_types.sample_data import SampleData
from q2_types.bowtie2 import Bowtie2Index

import q2_humann3
from q2_humann3._types import HumannDB, Nucleotide, Protein, Pathway
from q2_humann3._format import HumannDbFileFormat, HumannDbDirFormat

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
plugin.register_formats(HumannDbDirFormat, HumannDbFileFormat)

plugin.register_semantic_type_to_format(
    HumannDB[Nucleotide | Protein | Pathway], HumannDbDirFormat
)

plugin.methods.register_function(
    function=q2_humann3.run,
    inputs={
        "demultiplexed_seqs": SampleData[SequencesWithQuality],
        "nucleotide_database": HumannDB[Nucleotide],
        "protein_database": HumannDB[Protein],
        "pathway_database": HumannDB[Pathway],
    },
    parameters={"threads": qiime2.plugin.Int},
    name="Characterize samples using HUMAnN2",
    outputs=[
        ("genefamilies", FeatureTable[Frequency]),  # type: ignore
        ("pathcoverage", FeatureTable[RelativeFrequency]),  # type: ignore
        ("pathabundance", FeatureTable[RelativeFrequency]),  # type: ignore
        ("taxonomy", FeatureTable[RelativeFrequency]),  # type: ignore
    ],
    description="Execute the HUMAnN2",
)
