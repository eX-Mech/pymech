from pathlib import Path

import vaex
from .neksuite import readnek
from .dataset import can_open_nek_dataset


def open_dataframe(path, **kwargs):
    """Helper function for opening a file as an :class:`vaex.dataframe.Dataframe`.

    Parameters
    ----------
    path : str
            Path to a field file (only Nek files are supported at the moment.)

    kwargs : dict
            Keyword arguments passed on to the compatible open function.

    """
    if can_open_nek_dataset(path):
        _open = _open_nek_dataframe
    else:
        raise NotImplementedError(f"Filetype: {Path(path).suffix} is not supported.")

    return _open(path, **kwargs)


def _open_nek_dataframe(path):
    field = readnek(path)
    return vaex.concat(
        vaex.from_arrays(
            x=elem.pos[0].flatten(),
            y=elem.pos[1].flatten(),
            z=elem.pos[2].flatten(),
            ux=elem.vel[0].flatten(),
            uy=elem.vel[1].flatten(),
            uz=elem.vel[2].flatten(),
        )
        for elem in field.elem
    )
