# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

intersphinx-python:
	@python -msphinx.ext.intersphinx https://docs.python.org/3/objects.inv

intersphinx-nek:
	@python -msphinx.ext.intersphinx https://nek5000.github.io/NekDoc/objects.inv

intersphinx-xr:
	@python -msphinx.ext.intersphinx https://xarray.pydata.org/en/stable/objects.inv

autobuild:
	@sphinx-autobuild --watch ../pymech --re-ignore '_build\/.*' "$(SOURCEDIR)" "$(BUILDDIR)/html"
