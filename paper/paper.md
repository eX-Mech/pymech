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
date: 16 November 2020
bibliography: paper.bib
---

# Summary

- Features and functionalities:
  - Reads raw data as is without interpolation
  - Reasonably quick for small to medium sized files (quantify?)
  - Handles esoteric file formats used in open source solvers

- Describe exadata based on NumPy [@harrisArrayProgrammingNumPy2020]
- Dataset interface based on Xarray [@hoyerXarrayNDLabeled2017]


- State of the field:
  - Alternatives for mesh
  - Alternatives to read and visualize via VTK: Paraview (mention Python
    interface), Visit, turbulucid
  - Similar: [nekmatlab](https://github.com/nfabbiane/nekmatlab)


# Acknowledgements

We acknowledge ...

# References
