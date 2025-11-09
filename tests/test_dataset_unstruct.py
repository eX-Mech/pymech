import pymech as pm


def test_unstruct_dataset(test_data_dir):
    pm.open_dataset(
        test_data_dir / "nek" / "naca" / "naca0.f00001", mesh_type="UNSTRUCT"
    )
