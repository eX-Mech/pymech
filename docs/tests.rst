.. _tests:

Tests
=====

To test locally, set up a virtual environment and then::

  git clone https://github.com/eX-Mech/pymech
  cd pymech/
  pip install -e '.[tests]'
  pytest -s -v


Continuous Integration
----------------------

Tests are run automatically at every `push` to the repository.

The tests are run on TravisCI_ for Python versions between 3.6 and 3.8.

.. image:: https://travis-ci.org/eX-Mech/pymech.svg?branch=master
   :target: https://travis-ci.org/eX-Mech/pymech

The code coverage is also automatically checked by Coverall_.

.. image:: https://coveralls.io/repos/github/eX-Mech/pymech/badge.svg?branch=master
   :target: https://coveralls.io/github/eX-Mech/pymech


.. External links:

.. _TravisCI: https://travis-ci.org/eX-Mech/pymech
.. _Coverall: https://coveralls.io/github/eX-Mech/pymech
