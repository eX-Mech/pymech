import setuptools


setuptools.setup(
    use_scm_version={
        'write_to': 'pymech/_version.py',
        'write_to_template': (
            '__version__ = "{version}"\n'
        ),
    }
)
