[![Travis-CI Build Status](https://travis-ci.org/mpiktas/retacoro.svg?branch=master)](https://travis-ci.org/mpiktas/retacoro)
[![Coverage Status](https://img.shields.io/codecov/c/github/mpiktas/retacoro/master.svg)](https://codecov.io/github/mpiktas/retacoro?branch=master)

# retacoro

The **retacoro** package provides a method to transform table to another table with given column and row sum totals while preserving the structure of the initial table. 

The mathematical ideas behind the package are described in the following two mathoverflow questions:

http://mathoverflow.net/questions/36174/multinomial-transformation-for-matrices

http://mathoverflow.net/questions/156983/injectivity-of-matrix-fingerprint/157174

To install the package use the following command:

```
library(devtools)
install_github("mpiktas/retacoro", build_vignettes = TRUE)
```

