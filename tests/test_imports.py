from types import FunctionType


def test_top_level():
    import pymech as pm

    assert isinstance(pm.readnek, FunctionType)
    assert isinstance(pm.writenek, FunctionType)
    assert isinstance(pm.readre2, FunctionType)
    assert isinstance(pm.readrea, FunctionType)
    assert isinstance(pm.writere2, FunctionType)
    assert isinstance(pm.writerea, FunctionType)
    assert isinstance(pm.readplane, FunctionType)
    assert isinstance(pm.readdns, FunctionType)
    assert isinstance(pm.open_dataset, FunctionType)
    assert isinstance(pm.open_mfdataset, FunctionType)
