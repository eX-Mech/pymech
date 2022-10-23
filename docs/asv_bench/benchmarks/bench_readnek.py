from pymech.neksuite import readnek
from . import Channel3D


class ReadnekSuite(Channel3D):
    def time_read(self):
        readnek(self.file)

    def time_read_skip_geom(self):
        readnek(self.file, skip_vars=("x", "y", "z"))

    def time_read_skip_ux_uy(self):
        readnek(self.file, skip_vars=("ux", "uy"))
