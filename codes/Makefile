## This is the codes subdirectory

current: target
-include target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt"

######################################################################

install_packages.Rout: install_packages.R
	$(wrapR)

autopipeR = defined
Sources += $(wildcard *.R)

params.Rout: params.R 

SIRfunctions.Rout: SIRfunctions.R

sir_plot.Rout: sir_plot.R SIRfunctions.rda params.rda

modeldefs.tex: sir_plot.Rout ; 

Ignore += R0contour_TTI.pdf R0contour_random.pdf R0contour_TTI2.pdf R0contour_random2.pdf
R0contour_TTI.pdf R0contour_random.pdf R0contour_TTI2.pdf R0contour_random2.pdf: sir_plot.Rout ;

Ignore += targeted_test_df.csv random_test_df.csv
targeted_test_df.csv random_test_df.csv: sir_plot.Rout ;

######################################################################

### Makestuff

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff

Makefile: makestuff/Makefile
makestuff/Makefile:
	cd .. && $(MAKE) makestuff
	ln -s ../makestuff .

-include makestuff/os.mk

-include makestuff/pipeR.mk

-include makestuff/git.mk
-include makestuff/visual.mk
