# Developer Guide

(get-started)=

## Get Started!

Ready to contribute? Here's how to set up Pymech for local development.

1. Fork the `pymech` repo on GitHub: <https://github.com/eX-Mech/pymech/fork>
2. Clone your fork locally:
    ```{tab} HTTPS
        git clone --recursive https://github.com/your_name_here/pymech.git
    ```
    ```{tab} SSH
        git clone --recursive git@github.com:your_name_here/pymech.git
    ```
3. To set up for local development, first, create a [virtual
   environment](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments):
    ```{tab} Unix
        cd pymech/
        python3 -m venv venv
        source venv/bin/activate
    ```
    ```{tab} Windows
        cd pymech/
        py -m venv venv
        .\venv\Scripts\activate
    ```
    and then install `pymech` in editable mode along with development dependencies:
    ```
    pip install -e ".[dev]"
    pre-commit install
    ```
4. Create a branch for local development:

       git checkout -b name-of-your-bugfix-or-feature

    Now you can make your changes locally.

5. When you're done making changes, check that your changes pass the tests:

       pytest

6. Commit your changes and push your branch to GitHub. The pre-commit hooks
   installed in step 3 would automatically ensure that the changes are
   formatted:

       git add .
       git commit -m "Your detailed description of your changes."
       git push origin name-of-your-bugfix-or-feature

7.  Submit a pull request through the GitHub website.

## Pull Request Guidelines

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should be covered by unit tests. Add more unit tests if the
   code coverage drops.
2. If the pull request adds functionality, the docs should be updated. Describe
   the new functionality into a function with a docstring.
3. The pull request should pass the CI. Check
   <https://github.com/eX-Mech/pymech/actions>
   and make sure that the tests pass for all supported Python versions.

```{toctree}
internals
tests
release
```
