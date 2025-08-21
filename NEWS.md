# sgdGMF 1.0.2
* `sgdgmf.fit` : implemented orthogonality between covariates and latent variables
* `orthogonalize` (new function) : implemented orthogonality between covariates and latent variables
* `sgdgmf.fit` : implemented the possibility to not save a copy of the data and fitted values
* `set.control.airwls` : introduced new argument `savedata` to specify of store a copy of the data or not
* `set.control.newton` : introduced new argument `savedata` to specify of store a copy of the data or not
* `set.control.coord.sgd` : introduced new argument `savedata` to specify of store a copy of the data or not
* `set.control.block.sgd` : introduced new argument `savedata` to specify of store a copy of the data or not
* `storedata` (new function) : implemented ex-post inclusion of data in a generic object
* `storedata.sgdgmf` (new method) : implemented ex-post inclusion of data in a fitted `sgdgmf` object
* `sgdgmf.init` : implemented `method = "light"` and improved the memory usage
* `sgdgmf.init.light` (new function) : implemented an memory efficient version of `sgdgmf.init.ols` with `type = "link"`

# sgdGMF 1.0.1

* `sgdgmf.rank` : changed default method from `"onatski"` to `"evr"` method
* `eigengap.evr` (new function) : implemented the eigenvalue ratio method for rank selection
* `eigengap.onatski` : fixed bug occurring when no optimal rank can be selected
* added option `CXXFLAGS = $(CXXFLAGS) -Os` to `Makevars` and `Makevars.win` files to optimize the memory space usad by compiled C++ files
