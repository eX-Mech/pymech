Release checklist for maintainers
=================================

Via GitHub Actions
------------------

- Update :doc:`changelog` with version comparison links and detailed description
- Update version and date in ``CITATION.cff``.
- Make a `new release`_.
- Inspect upload in TestPyPI_.
- Execute manual `deploy workflow`_ to download from TestPyPI_, run tests and
  publish to PyPI_.

.. _new release: https://github.com/eX-Mech/pymech/releases/new
.. _deploy workflow: https://github.com/eX-Mech/pymech/actions/workflows/deploy.yaml

Manual method (not recommended)
-------------------------------

.. note::

   For demonstration's sake, we assume that the next version is ``$VERSION``
   and the package name is ``$PACKAGE``.

- Install nox::

      pip install nox

- Ensure tests pass locally and on CI::

      nox -s tests

- Compile documentation::

      nox -s docs

- Commit changes and make an annotated tag::

      git commit
      git tag -a $VERSION

- Build and upload to TestPyPI_::

      nox -s testpypi

- Download, test and upload to PyPI_::

      nox -s pypi

- Upload to repository::

      git push --follow-tags --atomic origin main

.. _twine: https://twine.readthedocs.io/en/latest/
.. _TestPyPI: https://test.pypi.org/project/pymech/
.. _PyPI: https://pypi.org/project/pymech/
