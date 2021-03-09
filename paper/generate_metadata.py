"""Generate pandoc metadata from codemeta.json"""
import sys
from collections import OrderedDict
from json import load
from pathlib import Path

import ruamel.yaml


here = Path.cwd()

with (here.parent / "codemeta.json").open() as fp:
    cm = load(fp)


#  print(cm)


def ordered_set(iterable):
    return OrderedDict.fromkeys(iterable).keys()


uniq_affiliations = list(ordered_set(author["affiliation"] for author in cm["author"]))


def make_author(author):
    return {
        "name": " ".join([author["givenName"], author["familyName"]]),
        "orcid": author.get("@id", ""),
        "affiliation": uniq_affiliations.index(author.get("affiliation", "")) + 1,
    }


def make_affiliations():
    return [{"name": name, "index": i+1} for i, name in enumerate(uniq_affiliations)]


# Prepare data for metadata output
data = {
    "title": cm.get("title", ""),
    "tags": cm.get("keywords", []),
    "authors": [make_author(_) for _ in cm.get("author", [{}])],
    "affiliation": make_affiliations(),
    "date": cm.get("dateModified", ""),
    "bibliography": "paper.bib",
}


# FIXME: New API wraps dictionaries
# yaml = ruamel.yaml.YAML(typ='safe')
# yaml.compact(seq_seq=False, seq_map=False)
# yaml.allow_unicode = True
# yaml.indent(sequence=4, offset=2, mapping=4)
def dump(data, dest):
    # yaml.dump(data, dest)
    ruamel.yaml.dump(
		data,
		dest,
		explicit_start=True,
		allow_unicode=True,
		Dumper=ruamel.yaml.RoundTripDumper
	)


dest = sys.stdout
dump(data, dest)
#  output = here / "metadata.yaml"
#  dest = output.open("w")
#  with output.open("w") as dest:
#      dump(data, dest)
