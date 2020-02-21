import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='pymech',
    version='1.2.0',
    author='Jacopo Canton, Nicolo Fabbiane and Guillaume Chauvat',
    author_email='jacopo.canton@gmail.com',
    description='A Python suite of routines for Nek5000 and Simson.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/jcanton/pymech',
    packages=['pymech'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
)
