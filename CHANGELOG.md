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

[Unreleased]: https://github.com/jcanton/pymech/compare/1.3.1...HEAD
[1.3.1]: https://github.com/jcanton/pymech/compare/1.3.0.post0...1.3.1
[1.3.0.post0]: https://github.com/jcanton/pymech/compare/1.3.0...1.3.0.post0
[1.3.0]: https://github.com/jcanton/pymech/compare/1.2.0...1.3.0
[1.2.0]: https://github.com/jcanton/pymech/releases/tag/1.2.0
