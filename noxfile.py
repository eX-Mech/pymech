"""Task runner for the developer

Usage
-----

   nox -l

   nox -s <session>

   nox -k <keyword>

"""

import os
import shlex
from functools import partial
from pathlib import Path
from shutil import rmtree

import nox
from nox_pdm import session

PACKAGE = "pymech"
CWD = Path.cwd()

TEST_ENV_VARS = {}
if os.getenv("CI"):
    TEST_ENV_VARS["PYTEST_ADDOPTS"] = "--color=yes"

no_venv_session = partial(session, venv_backend="none")
pytest_cmd = "python -m pytest".split()
mypy_cmd = "python -m mypy".split()

# nox.options.sessions = ["tests", "types"]


def run_ext(session, cmd):
    """Run an external command, i.e. outside a nox managed virtual envionment"""
    session.run(*shlex.split(cmd), external=True)


@session
def tests(session):
    """Execute unit-tests using pytest"""
    session.install(".[tests]")
    session.run(
        *pytest_cmd,
        *session.posargs,
        env=TEST_ENV_VARS,
    )


@session(name="tests-cov-vtk", python=["3.9", "3.10", "3.11", "3.12"])
def tests_cov_vtk(session):
    """Execute unit-tests using pytest+pytest-cov+VTK dependencies"""
    session.install(".[tests,vtk]")
    session.run(
        *pytest_cmd,
        "--cov",
        "--cov-append",
        "--cov-config=pyproject.toml",
        "--no-cov-on-fail",
        "--cov-report=term-missing",
        *session.posargs,
        "tests/test_vtk.py",
        env=TEST_ENV_VARS,
    )


@no_venv_session(name="tests-cov")
def tests_cov(session):
    """Execute unit-tests using pytest+pytest-cov"""
    session.notify(
        "tests",
        [
            "--cov",
            "--cov-config=pyproject.toml",
            "--no-cov-on-fail",
            "--cov-report=term-missing",
            *session.posargs,
        ],
    )


@session
def types(session):
    """Execute type-checking using mypy"""
    if (CWD / "src").exists():
        test_source = "src"
    else:
        test_source = PACKAGE

    session.install(".[types]")
    session.run(*mypy_cmd, *session.posargs, test_source, "tests")


@session(name="coverage-html")
def coverage_html(session, nox=False):
    """Generate coverage report in HTML. Requires `tests-cov` session."""
    report = Path.cwd() / ".coverage" / "html" / "index.html"
    session.install("coverage[toml]")
    session.run("coverage", "html")

    print("Code coverage analysis complete. View detailed report:")
    print(f"file://{report}")


@no_venv_session(name="format")
def format_(session):
    """Run pre-commit hooks on all files to set and lint code-format"""
    run_ext(session, "pre-commit install")
    run_ext(session, "pre-commit run --all-files")


@session
def lint(session):
    """Run pre-commit hooks on files which differ in the current branch from origin/HEAD."""
    remote = "origin/HEAD" if not session.posargs else session.posargs[0]
    session.install("pre-commit")
    session.run("pre-commit", "install")
    session.run("pre-commit", "run", "--from-ref", remote, "--to-ref", "HEAD")


def _prepare_docs_session(session):
    session.install(".[docs]")
    session.chdir("docs")

    build_dir = Path.cwd() / "_build"
    source_dir = "."
    output_dir = str(build_dir.resolve() / "html")
    return source_dir, output_dir


@session
def docs(session):
    """Build documentation using Sphinx."""
    source, output = _prepare_docs_session(session)
    session.run(
        "python", "-m", "sphinx", "-b", "html", *session.posargs, source, output
    )  # Same as sphinx-build
    print("Build finished.")
    print(f"file://{output}/index.html")


@session(name="docs-autobuild")
def docs_autobuild(session):
    """Build documentation using sphinx-autobuild."""
    source, output = _prepare_docs_session(session)
    session.run(
        "python",
        "-m",
        "sphinx_autobuild",
        "--watch",
        "../src",
        "--re-ignore",
        r"(_build|generated)\/.*",
        source,
        output,
    )  # Same as sphinx-autobuild
    print("Build finished.")
    print(f"file://{output}/index.html")


@no_venv_session
def testpypi(session):
    """Release clean, build, upload to TestPyPI"""
    session.notify("release-clean")
    session.notify("release-build")
    session.notify("release-upload", ["--repository", "testpypi"])


@no_venv_session
def pypi(session):
    """Release clean, download from TestPyPI, test, upload to PyPI"""
    session.notify("release-clean")
    # NOTE: parametrizing dist_type ends up in erraneous deduplication of sessions
    # by nox
    for dist_type in ("no-binary", "only-binary"):
        session.notify(f"download-testpypi(dist_type='{dist_type}')")
        session.notify(f"release-tests(dist_type='{dist_type}')")
    session.notify("release-upload", ["--repository", "pypi"])


@session(name="download-testpypi")
@nox.parametrize("dist_type", ["no-binary", "only-binary"])
def download_testpypi(session, dist_type):
    """Download from TestPyPI and run tests"""
    (Path.cwd() / "dist").mkdir(exist_ok=True)
    session.chdir("./dist")

    git_tags = session.run(
        "git", "tag", "--list", "--sort=taggerdate", external=True, silent=True
    )
    latest_version = git_tags.splitlines()[-1]
    spec = f"{PACKAGE}=={latest_version}"
    session.run(
        "python",
        "-m",
        "pip",
        "index",
        "versions",
        "--index-url",
        "https://test.pypi.org/simple",
        "--pre",
        PACKAGE,
    )
    session.run(
        "python",
        "-m",
        "pip",
        "download",
        "--index-url",
        "https://test.pypi.org/simple",
        "--extra-index-url",
        "https://pypi.org/simple",
        "--pre",
        "--no-deps",
        f"--{dist_type}",
        PACKAGE,
        spec,
    )


@session(name="release-tests")
@nox.parametrize("dist_type", ["no-binary", "only-binary"])
def release_tests(session, dist_type):
    """Execute test suite with build / downloaded package in ./dist"""
    if dist_type == "only-binary":
        pattern = "*.whl"
    else:
        pattern = "*.tar.gz"

    dist_packages = [str(p) for p in Path("./dist").glob(pattern)]
    session.install(*dist_packages)
    session.install("pytest")

    session.run(
        *pytest_cmd,
        env=TEST_ENV_VARS,
    )


@no_venv_session(name="release-clean")
def release_clean(session):
    """Remove build and dist directories"""
    session.log("Removing build and dist")
    rmtree("./build/", ignore_errors=True)
    rmtree("./dist/", ignore_errors=True)


@session(name="release-build")
def release_build(session):
    """Build package into dist."""
    session.install("build")
    session.run("python", "-m", "build")


@session(name="release-upload")
def release_upload(session):
    """Upload dist/* to repository testpypi (default, must be configured in ~/.pypirc).
    Also accepts positional arguments to `twine upload` command.

    """
    session.install("twine")
    session.run("twine", "check", "--strict", "dist/*")
    args = session.posargs

    # See
    # https://pypi.org/help/#apitoken and
    # https://twine.readthedocs.io/en/latest/#environment-variables
    env = {"TWINE_USERNAME": "__token__"}

    test_pypi_token = os.getenv("TEST_PYPI_TOKEN")
    pypi_token = os.getenv("PYPI_TOKEN")
    if "testpypi" in args and test_pypi_token:
        env["TWINE_PASSWORD"] = test_pypi_token
    elif pypi_token:
        env["TWINE_PASSWORD"] = pypi_token

    session.run("twine", "upload", *args, "dist/*", env=env)
