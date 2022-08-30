import qiime2.plugin
from q2_types.feature_table import FeatureTable, Frequency, RelativeFrequency
from q2_types.per_sample_sequences import SequencesWithQuality
from q2_types.sample_data import SampleData

import q2_humann3

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


plugin.methods.register_function(
    function=q2_humann3.run,
    inputs={"demultiplexed_seqs": SampleData[SequencesWithQuality]},  # type: ignore
    parameters={"threads": qiime2.plugin.Int},
    name="Characterize samples using HUMAnN2",
    outputs=[
        ("genefamilies", FeatureTable[Frequency]),  # type: ignore
        ("pathcoverage", FeatureTable[RelativeFrequency]),  # type: ignore
        ("pathabundance", FeatureTable[RelativeFrequency]),  # type: ignore
    ],
    description="Execute the HUMAnN2",
)
