## This is the Gharouni et al. SIR_testing_model

current: target
-include target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt"

######################################################################

subdirs += codes

colddirs += codes

alldirs += $(subdirs)

Ignore += $(alldirs)

######################################################################

## Manuscript

Sources += $(wildcard *.bib *.R) README.md
Sources += manuscript_BMB.tex

manuscript_BMB.pdf: manuscript_BMB.tex

######################################################################

## Download journal-specific format files
## Probably more elegant to do this with some sort of install
springer += spbasic.bst svjour3.cls svglov3.clo
springer_raw = https://raw.githubusercontent.com/latextemplates/svjour/master/

$(springer):
	wget -O $@ "$(springer_raw)/$@"

manuscript_BMB.pdf: $(springer)

######################################################################

## compartmental flowchart in ipe
## sudo apt-get install ipe
Sources += pix/sir_comp.ipe
Ignore += pix/sir_comp.pdf

Ignore += pix/sir_comp.pdf
pix/sir_comp.pdf: pix/sir_comp.ipe
	ipetoipe -pdf $< $@

######################################################################

### Makestuff

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff

Makefile: makestuff/Makefile
makestuff/Makefile:
	git clone $(msrepo)/makestuff
	ls makestuff/Makefile

-include makestuff/os.mk

-include makestuff/texi.mk
-include makestuff/pipeR.mk
-include makestuff/hotcold.mk

-include makestuff/git.mk
-include makestuff/visual.mk

