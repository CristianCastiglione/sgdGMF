# sgdGMF 1.0.1

* `sgdgmf.rank` : changed default method from `"onatski"` to `"evr"` method
* `eigengap.evr` (new function) : implemented the eigenvalue ratio method for rank selection
* `eigengap.onatski` : fixed bug occurring when no optimal rank can be selected
* added option `CXXFLAGS = $(CXXFLAGS) -Os` to `Makevars` and `Makevars.win` files to optimize the memory space usad by compiled C++ files
