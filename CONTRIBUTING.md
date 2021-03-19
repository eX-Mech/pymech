# Contributing

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

## Types of Contributions

### Report Bugs

Report bugs at <https://github.com/eX-Mech/pymech/issues>.

If you are reporting a bug, please include:

  - Your operating system name and version.
  - Any details about your local setup that might be helpful in
    troubleshooting.
  - Detailed steps to reproduce the bug.

### Fix Bugs

Look through the GitHub issues for bugs. Anything tagged with "bug" and
"help wanted" is open to whoever wants to implement it.

### Implement Features

Look through the GitHub issues for features. Anything tagged with
"enhancement" and "help wanted" is open to whoever wants to implement
it.

### Write Documentation

Pymech could always use more documentation,
whether as part of the official Pymech docs,
in docstrings, or even on the web in blog posts, articles, and such.

### Submit Feedback

The best way to send feedback is to file an issue at
<https://github.com/eX-Mech/pymech/issues>.

If you are proposing a feature:

  - Explain in detail how it would work.
  - Keep the scope as narrow as possible, to make it easier to
    implement.
  - Remember that this is a volunteer-driven project, and that
    contributions are welcome :)

## Get Started!

Ready to contribute? Here's how to set up Pymech for local development.

1. Fork the `pymech` repo on GitHub.
2. Clone your fork locally:

        $ git clone git@github.com:your_name_here/pymech.git

3. Install your local copy into a virtualenv. This is how you set up your fork for local development:

        $ cd pymech/
        $ python -m venv venv
        $ source venv/bin/activate
        $ pip install -e ".[dev]"
        $ pre-commit install

4. Create a branch for local development:

        $ git checkout -b name-of-your-bugfix-or-feature

    Now you can make your changes locally.

5. When you're done making changes, check that your changes pass the linters
   `flake8` and `darker` as well as the tests:

        $ flake8 pymech
        $ black --check pymech
        $ pytest

    To get `flake8` and `darker`, just pip install them into your virtualenv.
    These checks would be ensured via the pre-commit hooks while committing.

6. Commit your changes and push your branch to GitHub:

        $ git add .
        $ git commit -m "Your detailed description of your changes."
        $ git push origin name-of-your-bugfix-or-feature

7.  Submit a pull request through the GitHub website.

## Pull Request Guidelines

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring.
3. The pull request should pass the CI. Check
   <https://travis-ci.org/eX-Mech/pymech/pull_requests> and
   <https://github.com/eX-Mech/pymech/actions>
   and make sure that the tests pass for all supported Python versions.
