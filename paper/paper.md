---
title: 'Pymech: A Python suite of routines for Nek5000 and Simson'
tags:
- CFD
- fluid dynamics
- post-processing
- Nek5000
- SIMSON
- file I/O
authors:
- name: Ashwin Vishnu Mohanan
  affiliation: 1
  orcid: 0000-0002-2979-6327
- name: Guillaume Chauvat
  affiliation: 2
- name: Nicolo Fabbiane
  affiliation: 3
- name: Jacopo Canton
  affiliation: 4
affiliations:
- name: Stockholm University, Stockholm, Sweden
  index: 1
- name: KTH Royal Institute of Technology, Stockholm, Sweden
  index: 2
- name: ONERA, Paris, France
  index: 3
- name: ETH Zürich, Zürich, Switzerland
  index: 4
date: 17 March 2021
bibliography: paper.bib
---

# Summary

- Nek5000 [@NEK50002019] and SIMSON [@Chevalier.etalSIMSON2007]

- Features and functionalities:
  - Reads raw data as is without interpolation
  - Reasonably quick for small to medium sized files (quantify?)
  - Handles esoteric file formats used in open source solvers

- Describe exadata based on NumPy [@Harris.etalArray2020]
- Dataset interface based on Xarray [@Hoyer.HammanXarray2017]
- Can read for both small and very large files: only limited by available memory in the computing device. (Should we time reading and writing time with CPython and PyPy?)


- Functionality and real world uses:
  - Key features: Easy to understand interface to imports meshes and solution
    data as a Python-NumPy based data structure `exadata`.
  - Generate publication quality figures and visualizations.
  - Bridge between two solver to transfer solutions of Simson to Nek5000
  - A starting point and input data for post-processing algorithms
    written in Python.
  - Generate, (and soon extrude and combine) meshes
  - Initial condition generation
  - Cloning small solution data from simulation with assumed symmetry into full-3D calculations
  - In summary the package is an "enabler" for exploratory research work

- State of the field:

  - Alternatives for mesh:
    - Nek5000 comes with a few meshing tools to generate (`genbox`, `prex`,
      `pretext`), merge (`nekmerge`) and extrude (`n2to3`) meshes
      from the scratch. With exception to `genbox` the others output the legacy
      mesh format `.rea` (NOTE: confirm).
    - Mesh converters with input formats exodus, CGNS and Gmsh are also bundled
      with Nek5000.
  - Alternatives to read and visualize via VTK: Paraview
    [@AyachitParaView2015]
    interface), Visit [@Childs.etalContract2005] -- both provide Python
    interface, but the API is advanced. Alternatively, generate VTK files and
    later visualize in packages such as turbulucid or MayaVi.
  - Similar: [nekmatlab](https://github.com/nfabbiane/nekmatlab)
  - Alternatives to postprocessing:
    - Nek5000 in post-processing mode, can work in parallel. Performant. Re-use
      differentiation and interpolation operators. Disadvantage: must rely on
      existing implementation and examples to "reverse engineer"
      post-processing. Core API is not well-documented.
    - Use *Filters* in Paraview and equivalent (?) in Visit. Cumbersome to script
      by hand and inspect.


# Acknowledgements

We acknowledge ...

# References
