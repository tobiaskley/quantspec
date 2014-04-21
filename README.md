quantspec: Quantile-Based Spectral Analysis with R
==================================================

The aim of the `quantspec` package is to make methods for quantile-based spectral analysis of time series available to data analysts and researchers in statistics.

You can track (and contribute to) the development of `quantspec` at https://github.com/tobiaskley/quantspec. If you encounter unexpected behavior while using `quantspec`, please write an [email](mailto:tobias.kley@ruhr-uni-bochum.de) or file an [issue](http://github.com/tobiaskley/quantspec/issues).

## Getting started with ``quantspec``

First, if you have not done so already, install R from http://www.r-project.org (click on download R, select a location close to you, and download R for your platform). Once you have the latest version of R installed and started execute the following commands on the R shell:

 ```
 install.packages("devtools")
 devtools::install_github("tobiaskley/quantspec", ref="develop")
 ```

This will first install the R package ``devtools`` and then use it to install the latest (development) version of ``quantspec`` from the GitHub repository. In case you do not have LaTeX installed on your computer you may want to use

 ```
 devtools::install_github("tobiaskley/quantspec", ref="develop", build_vignette = FALSE)
 ```

to skip building the vignette.

Now that you have R and ``quantspec`` installed you can access all the functions available. To load the package and access the help files:

```
library(quantspec)
help("quantspec")
```

Three demos are available. They can be started by

```
demo("sp500")
demo("wheatprices")
demo("qar-simulation")
```

At the bottom of the online help page to the package you will find an index to all the help files available. If you did not skip building the vignette (a preprint of a paper including a tutorial and two worked examples) you can access it via

```
vignette("quantspec")
```
