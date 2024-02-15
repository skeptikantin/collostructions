### {collostructions}

The R package {collostructions} is an R implementation for collostructional analysis. The package contains functions to perform Simple, Distinctive, and Co-Varying Collexeme Analyses, as well as some useful functions for transformation and the creation of frequency lists. There are sample data sets for all functions with documentation.

### Installation

You can install the package from this repository like so (you need the packages `devtools` and `githubinstall` installed):

```
library(githubinstall)
githubinstall("skeptikantin/collostructions")
```

### Usage (quick guide)

The package is discussed in Episode X of my R tutorial, which you can find [here](https://www.youtube.com/playlist?list=PLIZN-827NSIONkLPWpjaFr0mlKacSLTRy). The tutorial was done quite a while ago with a much earlier version of the package that lacked many of the features it now has, but the principle is still the same.

To load the package:
```
library(collostructions)
```

For example, a Simple Collexeme Analysis for the go-V construction, the data of which is provided with the package for illustration:

```
data(goVerb)
collex(goVerb, 616336708)
```


### Older Versions

For the moment, there are no previous versions on here. If you would like to use an older version, please drop me an [email](mailto:susanne.flach@uzh.ch).
