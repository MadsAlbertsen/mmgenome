mmgenome
========

mmgenome: Tools for extracting genomes from metagenomes

## Changelog

### 0.6.0
#### Bugfix
 - mmsubset: The subsetting of essential genes was wrong, which affected the mmstats calculations. Fixed thanks to @russianconcussion.  

### 0.6.0
#### Enhancements
 - mmplot_network: new option to scale line width of links between scaffolds by a constant (`scale.links`).

### 0.5.0
#### Enhancements
 - mmsubset: new function to subset based on normal subset syntax. E.g. `gc > 50`.
 
#### Bugfixes
 - mmstats: Now rounding coverages to 2 decimals instead of one (Kirk3gaard).

### 0.4.3
 - mmstats: Error in summing scaffold lengths of large datasets due to integer overflow. Handled by converting integer to numeric (skarst).

### 0.4.2
 - mmplot: log.x and log.y is now false by default in order to function in mmplot_pairs (skarst).

### 0.4.1
 - mmload: Failed to load data if coverage was missing for some scaffolds. Fixed thanks to Matthew Neave (neavemj).

### 0.4.0
 - mmplot_network: It's now possible to to highlight scaffolds.
 - General: Information on essential genes is no longer required.

### 0.3.1
#### Bugfixes
 - mmload: Column names of additional datasets did not load as intended.
 
#### Enhancements 
 - mmplot: The variable "point.size" can be used to scale the points by a constant instead of the scaffold length.
 - mmimport: New simple function to handle import of the data through a Rmarkdown file.

### 0.3.0
#### Enhancements
 - mmplot: It's now possible to highlight subsets. Either by just suppling the subset dataset or using the names of the scaffolds.
 - mmplot: "color" support the argument "none" for coloring all scaffolds black.
 - mmextract: The function now supports to "include" or "exclude" scaffolds based on a list of scaffold names.

### 0.2.0
First stable release of the mmgenome package

### 0.1.0
Initial release of the mmgenome package
