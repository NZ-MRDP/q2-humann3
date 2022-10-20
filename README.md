## Installing q2-humann3

q2-humann3 must be installed in a QIIME2 environment, for help installing it see the [documentation](https://docs.qiime2.org/2022.8/install/native/#install-qiime-2-within-a-conda-environment)

Inside the QIIME2 conda environment

#### install Humann3

```bash
pip install humann --no-binary :all:
```

#### Install q2-humann3

from the root q2-humann3 directory

```bash
pip install .
qiime dev refresh-cache
```

####

Ensure it is installed

```bash
qiime humann3 run --help
```

### Helpful commands for testing q2-humann3

There is data in the `.assets` directory that can be used for testing

The equivalent humann3 command that we will be running in QIIME is the following:

```
humann3 -i assets/humann3-test-data/test-fastq_0_L001_R1_001.fastq.gz -o foo-out --output-format biom --remove-column-description-output --pathways-database assets/humann3-test-data/dbs/pathways_DEMO/metacyc_reactions_level4ec_only.uniref.bz2,assets/humann3-test-data/dbs/pathways_DEMO/metacyc_pathways_structured_filtered_v24 --input-format fastq.gz --nucleotide-database assets/humann3-test-data/dbs/chocophlan_DEMO --protein-database assets/humann3-test-data/dbs/uniref_DEMO
```

Notice there are four files outside of the sequence data that is required. These files are provided by default when running humann3, but we want to be explicit about what we are passing.
Surprises are good for birthdays, not for software.

There are two pathways files that we'll pass
metacyc_reactions_level4ec_only.uniref.bz2
metacyc_pathways_structured_filtered_v24

a directory of nucleotide sequences
chocophlan_DEMO

And a directory of protein sequences
uniref_DEMO

### Convert the files into QIIME compatible formats

**Pathway Mapping File**

```
qiime tools import --type HumannDB[PathwayMapping] --input-path assets/humann3-test-data/dbs/pathways_DEMO/metacyc_reactions_level4ec_only.uniref.gz --output-path assets/qiime-test-data/pathway-mapping
```

**Pathways File**

```
qiime tools import --type HumannDB[Pathway] --input-path assets/humann3-test-data/dbs/pathways_DEMO/metacyc_pathways_structured_filtered_v24 --output-path assets/qiime-test-data/pathways
```

**Chocophlan**

```
qiime tools import --type HumannDB[Nucleotide] --input-path assets/humann3-test-data/dbs/chocophlan_DEMO/ --output-path assets/qiime-test-data/nucleotides
```

**Uniref**

```
qiime tools import --type HumannDB[Protein] --input-path assets/humann3-test-data/dbs/pathways_DEMO/metacyc_pathways_structured_filtered_v24 --output-path assets/qiime-test-data/proteins
```

## Run the default q2-humann3 command:

```bash
qiime humann3 run --i-demultiplexed-seqs assets/qiime-test-data/trimmed-seqs.qza --i-nucleotide-database assets/qiime-test-data/nucleotides.qza --i-protein-database assets/qiime-test-data/proteins.qza --i-pathway-database assets/qiime-test-data/pathways.qza --i-pathway-mapping assets/qiime-test-data/pathway-mapping.qza --o-genefamilies gene-families --o-pathcoverage coverage --o-pathabundance abundance --o-taxonomy taxonomy
```
