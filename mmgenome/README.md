mmgenome
========

mmgenome: Tools for extracting genomes from metagenomes

## Changelog

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
Fist stable release of the mmgenome package

### 0.1.0
Initial realease of the mmgenome package