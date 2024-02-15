### {collostructions}

The R package {collostructions} is an R implementation for collostructional analysis. The package contains functions to perform Simple, Distinctive, and Co-Varying Collexeme Analyses, as well as some useful functions for transformation and the creation of frequency lists. There are sample data sets for all functions with documentation. The current version is v.0.2.0 (09-Feb-2021) that is compatible with all R4.x versions on all major operating systems. If you spot an error, [drop me a line](mailto:susanne.flach@uzh.ch) or use the issues section here and I'll fix it asap.

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

### References

**The paper that uses the sample code for `collex()`:**

Flach, Susanne. 2015. Let's go look at usage: A constructional approach to formal constraints on go-VERB. In Thomas Herbst & Peter Uhrig (eds.), *Yearbook of the German Cognitive Linguistics Association* (Volume 3), 231-252. Berlin: De Gruyter Mouton. doi:10.1515/gcla-2015-0013.

**A paper that shows the *use* collostructional attraction in experimental data:**

Flach, Susanne. 2020. Schemas and the frequency/acceptability mismatch: Corpus distribution predicts sentence judgements. *Cognitive Linguistics 31(4)*. 609-645. [Link (free access)](https://www.degruyter.com/view/journals/cogl/31/4/article-p609.xml)



**For an explanation of the method:**

Stefanowitsch, Anatol & Stefan Th. Gries. 2003. Collostructions: Investigating the interaction of words and constructions. *International Journal of Corpus Linguistics* 8(2). 209-243.

