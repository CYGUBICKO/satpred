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
## Implemented methods
### rfsrc.satpred: random survival forest
### gbm.satpred: gradient boosted trees
### deepsurv.satpred: survival neural network

satpred.Rout: R/satpred.R
utilities.Rout: R/utilities.R
gbm_satpred.Rout: R/gbm_satpred.R
deepsurv_satpred.Rout: R/deepsurv_satpred.R
methods.Rout: R/methods.R
pkgsExport.Rout: R/pkgsExport.R
mod_cv.Rout: R/mod_cv.R
mod_train.Rout: R/mod_train.R
posthocfuns.Rout: R/posthocfuns.R
satpredplots.Rout: R/satpredplots.R

######################################################################

## Examples
data.Rout: data.R

### Survival forest
examples_rfsrc.Rout: examples_rfsrc.R data.rda

### GBM
examples_gbm.Rout: examples_gbm.R data.rda

### deepsurv
examples_deepsurv.Rout: examples_deepsurv.R data.rda
deleteME.Rout: deleteME.R data.rda
deleteME_load.Rout: deleteME_load.R data.rda deleteME.rda

### Compare
examples_compare.Rout: examples_compare.R data.rda examples_rfsrc.rda examples_gbm.rda examples_deepsurv.rda

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
