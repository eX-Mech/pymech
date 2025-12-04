import pytest

import pymech as pm


@pytest.mark.parametrize(
    "file", ["cbox0.fld", "cbox0_issue_48.fld", "cbox0_issue_56.fld"]
)
def test_cbox(file, test_data_dir):
    pm.open_dataset(test_data_dir / "nek" / file)
