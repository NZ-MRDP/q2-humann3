import itertools

from q2_types.bowtie2 import Bowtie2IndexFileFormat
from q2_types.per_sample_sequences import FastqGzFormat
from qiime2.plugin import model


class HumannDbFileFormat(model.BinaryFileFormat):
    def _validate_(self, *args):
        pass


class HumannDbDirFormat(model.DirectoryFormat):
    data = model.FileCollection(r".+\..+", format=HumannDbFileFormat)

    @data.set_path_maker
    def data_path_maker(self, file):

        return file


HumannDBSingleFileDirFormat = model.SingleFileDirectoryFormat(
    "HumannDBSingleFileDirFormat", "mapping.gz", HumannDbFileFormat
)

# TODO: make this generic for single file humann3


class Bowtie2IndexDirFmt2(model.DirectoryFormat):
    idx1 = model.File(r".+(?<!\.rev)\.1\.bt2l?", format=Bowtie2IndexFileFormat)
    idx2 = model.File(r".+(?<!\.rev)\.2\.bt2l?", format=Bowtie2IndexFileFormat)
    ref3 = model.File(r".+\.3\.bt2l?", format=Bowtie2IndexFileFormat)
    ref4 = model.File(r".+\.4\.bt2l?", format=Bowtie2IndexFileFormat)
    rev1 = model.File(r".+\.rev\.1\.bt2l?", format=Bowtie2IndexFileFormat)
    rev2 = model.File(r".+\.rev\.2\.bt2l?", format=Bowtie2IndexFileFormat)
    pkl = model.File(r".+\.pkl?", format=Bowtie2IndexFileFormat)

    def get_basename(self):
        paths = [str(x.relative_to(self.path)) for x in self.path.iterdir()]
        prefix = _get_prefix(paths)
        return prefix[:-1]  # trim trailing '.'


# SO: https://stackoverflow.com/a/6718380/579416
def _get_prefix(strings):
    def all_same(x):
        return all(x[0] == y for y in x)

    char_tuples = zip(*strings)
    prefix_tuples = itertools.takewhile(all_same, char_tuples)
    return "".join(x[0] for x in prefix_tuples)
