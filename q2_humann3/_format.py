from qiime2.plugin import model


class HumannDbFileFormat(model.BinaryFileFormat):
    def _validate_(self, *args):
        pass


class HumannDbDirFormat(model.DirectoryFormat):
    data = model.FileCollection(r".+\..+", format=HumannDbFileFormat)

    @data.set_path_maker
    def data_path_maker(self, file):

        return file
