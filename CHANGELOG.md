# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

<!--

### Added
### Changed
### Deprecated
### Removed
### Fixed
### Security

Type of changes
---------------

Added for new features.
Changed for changes in existing functionality.
Deprecated for soon-to-be removed features.
Removed for now removed features.
Fixed for any bug fixes.
Security in case of vulnerabilities.

-->

## [Unreleased]

### Added
- Functions `readre2` and `writere2` to read/write binary .re2 Nek5000 meshes

### Fixed
- indexing bug in the boundary conditions parameters in `readrea`

### Changed
- Ignore some invalid 'E' internal boundary conditions in `readrea`
  since they are written this way by re2torea and ignored by Nek5000.


## [1.3.3] - 2020-09-29

### Fixed
- Various fixes -- including writing element map, correct order for min/max
  metadata -- in `writenek`

### Added
- Function `exadata.merge` to merge meshes together and build proper connectivity

## [1.3.2] - 2020-09-23

### Fixed
- `writenek` detects system endianness and byte-swaps arrays, if needed, before
writing

### Removed
- Warnings while reading/writing scalars

## [1.3.1] - 2020-09-17

### Changed
- use ndarray.tofile() for faster output

## [1.3.0.post0] - 2020-07-16

### Fixed
- Packaging issue of sdist and wheel. Now uses `find_packages` instead of package name.

## [1.3.0] - 2020-07-16

### Added
- Dataset module which extends xarray

### Changed
- Faster `readnek` function uses less for loops
- Lazy load `exadata` limits as properties
- Autogenerate documentation using sphinx extension `autodoc`

## [1.2.0] - 2020-03-18

### Added
- License GPL v3 or later

### Changed
- Miscellaneous improvements in documentation, testing and packaging

[Unreleased]: https://github.com/jcanton/pymech/compare/1.3.3...HEAD
[1.3.3]: https://github.com/jcanton/pymech/compare/1.3.2...1.3.3
[1.3.2]: https://github.com/jcanton/pymech/compare/1.3.1...1.3.2
[1.3.1]: https://github.com/jcanton/pymech/compare/1.3.0.post0...1.3.1
[1.3.0.post0]: https://github.com/jcanton/pymech/compare/1.3.0...1.3.0.post0
[1.3.0]: https://github.com/jcanton/pymech/compare/1.2.0...1.3.0
[1.2.0]: https://github.com/jcanton/pymech/releases/tag/1.2.0
