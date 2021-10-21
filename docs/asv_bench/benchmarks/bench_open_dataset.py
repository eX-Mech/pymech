try:
    from pymech.dataset import open_dataset
except ImportError:
    NO_DATASET = True
else:
    NO_DATASET = False

from . import Channel3D


class OpenDatasetSuite(Channel3D):
    def setup(self):
        if NO_DATASET:
            raise NotImplementedError

        super().setup()

    def time_read(self):
        open_dataset(self.file)
