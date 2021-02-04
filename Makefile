## This is satpred
### Survival Analysis Trianing and PREDiction (satpred)

current: target
-include target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt notes.md"

######################################################################

Sources += $(wildcard *.R *.md. *.Rnw)
Sources += $(wildcard R/*.R)
Sources += $(wildcard man/*.Rd) NAMESPACE DESCRIPTION

autopipeR = defined

######################################################################
satpred.Rout: R/satpred.R
utilities.Rout: R/utilities.R
methods.Rout: R/methods.R
mod_cv.Rout: R/mod_cv.R
mod_train.Rout: R/mod_train.R
posthocfuns.Rout: R/posthocfuns.R


######################################################################

## Package installation and checks
Ignore += satpred_1*

build-package:
	R CMD build .

install-package:
	R CMD INSTALL satpred_1*

check-package:
	echo "devtools::check('.')" | R --slave

update-doc:
	echo "devtools::document('.')" | R --slave

######################################################################

### Makestuff

Sources += Makefile

## Sources += content.mk
## include content.mk

Ignore += makestuff
msrepo = https://github.com/dushoff

Makefile: makestuff/Makefile
makestuff/Makefile:
	git clone $(msrepo)/makestuff
	ls makestuff/Makefile

-include makestuff/os.mk

-include makestuff/pipeR.mk

-include makestuff/git.mk
-include makestuff/visual.mk

-include rnw.mk
