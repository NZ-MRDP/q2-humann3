from qiime2.plugin import model

from q2_types.per_sample_sequences import FastqGzFormat


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
