hmmgs_gene_targeted_assm.pdf: hmmgs_gene_targeted_assm.dvi
	dvipdf hmmgs_gene_targeted_assm.dvi

hmmgs_gene_targeted_assm.dvi: hmmgs_gene_targeted_assm.tex
	latex hmmgs_gene_targeted_assm.tex && bibtex hmmgs_gene_targeted_assm.aux && latex hmmgs_gene_targeted_assm.tex && latex hmmgs_gene_targeted_assm.tex

clean:
	rm *~ hmmgs_gene_targeted_assm.dvi hmmgs_gene_targeted_assm.bbl hmmgs_gene_targeted_assm.pdf hmmgs_gene_targeted_assm.log hmmgs_gene_targeted_assm.blg hmmgs_gene_targeted_assm.aux
