import pymech as pm


def test_unstruct_dataset(test_data_dir):
    pm.open_dataset(
        test_data_dir + "/nek" + "/naca" + "/naca(10x4).e", mesh_type="UNSTRUCT"
    )


test_unstruct_dataset("C:/Users/gaura/Desktop/pymech-test-data")
