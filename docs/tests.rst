.. _tests:

Tests
=====

To test locally, :ref:`set up a virtual environment as shown here
<get-started>`. After that, install Pymech along with test dependencies and run
the tests::

    pip install -e '.[tests]'
    pytest -s -v


Continuous Integration
----------------------

Tests are run automatically at every `push` and `pull request` to the repository.

The tests are run on GitHub_ for Python versions between |py_min_version| and up.

The code coverage is also logged at Coveralls_.

.. External links:

.. _GitHub: https://github.com/eX-Mech/pymech/actions
.. _Coveralls: https://coveralls.io/github/eX-Mech/pymech
