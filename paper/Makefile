draft: paper.pdf

paper.pdf: paper.md paper.bib
	pandoc -o paper.pdf -V documentclass=scrartcl -V geometry:margin=1in --pdf-engine=xelatex --citeproc paper.md
