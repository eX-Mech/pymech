Release checklist for maintainers
=================================

.. note::

   For demonstration's sake, we assume that the next version is ``$VERSION``
   and the package name is ``$PACKAGE``.

- Install in development mode::

      pip install -e '.[dev]'

- Ensure tests pass locally and on CI::

      pytest --cov --numprocesses=4

- Compile documentation::

      cd docs/
      make html

- Commit changes and tag a release::

      git commit
      git tag $VERSION

- Prepare source distribution package and wheel::

      python setup.py sdist bdist_wheel

- Verify the package with twine_::

      twine check dist/*

- Upload to TestPyPI_ and verify::

      twine upload --repository testpypi dist/*
      cd /tmp
      python -m venv testpypi
      source testpypi/bin/activate
      pip install \
          --index-url https://test.pypi.org/simple \
          --extra-index-url https://pypi.org/simple \
          $PACKAGE

- Upload to PyPI_ with twine_::

      twine upload dist/*

- Upload to repository::

      git push --tags && git push

.. _twine: https://twine.readthedocs.io/en/latest/
.. _TestPyPI: https://packaging.python.org/guides/using-testpypi/
.. _PyPI: https://pypi.org/
