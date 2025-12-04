# Pymech

[![PyPI](https://img.shields.io/pypi/v/pymech)](https://pypi.org/project/pymech/)
[![Build Status](https://img.shields.io/github/actions/workflow/status/eX-Mech/pymech/build.yaml?branch=main)](https://github.com/eX-Mech/pymech/actions)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/eX-Mech/pymech/main.svg)](https://results.pre-commit.ci/latest/github/eX-Mech/pymech/main)
[![Coverage Status](https://coveralls.io/repos/github/eX-Mech/pymech/badge.svg)](https://coveralls.io/github/eX-Mech/pymech)
[![Documentation Status](https://readthedocs.org/projects/pymech/badge/?version=latest)](http://pymech.readthedocs.org/en/stable/)
[![DOI](https://zenodo.org/badge/50511298.svg)](https://zenodo.org/badge/latestdoi/50511298)

A Python suite of routines for [Nek5000] and [SIMSON]. Install with:

```
pip install pymech
```

Read the full documentation at [Pymech doc](http://pymech.readthedocs.io/en/stable).

## Getting started

For some quick wins, download some sample data

```sh
curl -fLO https://raw.githubusercontent.com/eX-Mech/pymech-test-data/main/nek/channel3D_0.f00001
```

Fire up a Python / IPython console and execute:

```py
import matplotlib.pyplot as plt
import pymech as pm

ds = pm.open_dataset('channel3D_0.f00001')
ds.mean(['x', 'z']).ux.plot()
plt.show()
```

You should see something like

![](https://pymech.readthedocs.io/en/stable/_images/usage_37_1.png)

For an overview of how to use pymech, [follow the usage
documentation](https://pymech.readthedocs.io/en/stable/usage.html).

## Contributing

Found something that does not work? Want to add a new feature to pymech? Have a
look at [contributing
guidelines](https://pymech.readthedocs.io/en/stable/contributing.html).

[Nek5000]: https://nek5000.mcs.anl.gov/
[SIMSON]: https://github.com/KTH-Nek5000/SIMSON/
