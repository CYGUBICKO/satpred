## This is satpred
### Survival Analysis Trianing and PREDiction (satpred)

current: target
-include target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt notes.md"

setssh:
	git remote set-url origin git@github.com:cygubicko/satpred.git

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
ibrierscore.Rout: R/ibrierscore.R
tdcroc.Rout: R/tdcroc.R

######################################################################

## Examples
data.Rout: data.R
data_tdc.Rout: data_tdc.R

### Survival forest
examples_rfsrc.Rout: examples_rfsrc.R data.rda

### GBM
examples_gbm.Rout: examples_gbm.R data.rda

### GBM3
examples_gbm3.Rout: examples_gbm3.R data_tdc.rda

### deepsurv
examples_deepsurv.Rout: examples_deepsurv.R data.rda

### Compare
examples_compare.Rout: examples_compare.R data.rda examples_rfsrc.rda examples_gbm.rda examples_deepsurv.rda

## Compare glment and pcoxtime for tdc
### using glmnet example
examples_glmnet_data.Rout: examples_glmnet_data.R

## Compare glment and pcoxtime 
### using veteran data pseudo tdc
examples_glmnet_veteran.Rout: examples_glmnet_veteran.R

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

install:
	make update-doc && make build-package && make install-package

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
