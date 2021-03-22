[![Build Status](https://travis-ci.org/eX-Mech/pymech.svg?branch=master)](https://travis-ci.org/eX-Mech/pymech/builds)
[![Coverage Status](https://coveralls.io/repos/github/eX-Mech/pymech/badge.png?branch=master)](https://coveralls.io/github/eX-Mech/pymech?branch=master)
[![Documentation Status](https://readthedocs.org/projects/pymech/badge/?version=latest)](http://pymech.readthedocs.org/en/latest/?badge=latest)

# pymech

A Python suite of routines for *Nek5000* and *Simson*.

Install with `pip install pymech`

Read the full documentation at [Pymech doc](http://pymech.readthedocs.io).

## Getting started

For some quick wins, download some sample data

```sh
curl -LO https://github.com/eX-Mech/pymech/raw/master/tests/nek/channel3D_0.f00001
```

Fire up a Python console and execute:

```py
import pymech as pm
ds = pm.open_dataset('channel3D_0.f00001')
ds.mean(['x', 'z']).ux.plot()
```

You should see something like

![](https://pymech.readthedocs.io/en/latest/_images/usage_33_1.png)

## Contributing

Found something that does not work? Want to add a new feature to pymech? Have a
look at [contributing
guidelines](https://pymech.readthedocs.io/en/latest/contributing.html).
