from pymech.neksuite import readnek
from . import Channel3D


class ReadnekSuite(Channel3D):
    def time_read(self):
        readnek(self.file)
