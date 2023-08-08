# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.
from pathlib import Path

test_data = Path(__file__).parent.parent.parent.parent / "tests" / "data"


class Channel3D:
    def setup(self):
        self.file = test_data / "nek" / "channel3D_0.f00001"
