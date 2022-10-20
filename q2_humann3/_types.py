from qiime2.plugin import SemanticType

HumannDB = SemanticType("HumannDB", field_names=["annotation"])

Nucleotide = SemanticType("Nucleotide", variant_of=HumannDB.field["annotation"])
Protein = SemanticType("Protein", variant_of=HumannDB.field["annotation"])
Pathway = SemanticType("Pathway", variant_of=HumannDB.field["annotation"])
PathwayMapping = SemanticType("PathwayMapping", variant_of=HumannDB.field["annotation"])
ReferenceNameMapping = SemanticType(
    "ReferenceNameMapping", variant_of=HumannDB.field["annotation"]
)
