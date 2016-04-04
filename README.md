Source code accompanying the paper "Beyond Cognacy: Historical relations between words"
=======================================================================================

## Requires

* python
* networks (http://networkx.org)
* lingpy (http://lingpy.org, to re-calculate the data files)

## Replicate the analyses

In order to replicate the analyses in the paper from a unix terminal, make sure you have LingPy installed and then use the MAKE file:

```bash
$ sh MAKE
```

If you want to reproduce the whole experiment, including the creation of the transition matrices, enter the following command in the terminal:

```bash
$ sh MAKE prepare
```

This should run all necessary analyses and re-produce the output which you find already in this folder.

If you further want to plot all data for the PCD analysis to retrive html-output in which all state transitions can be inspected on a tree, you may provide the "plot" argument to the shell file:

```bash
$ sh MAKE plot
```

The plots can be found in a folder that will be specifically created. Open the first file (index.html) in the folder. From there, you can easily inspect all patterns.

## Files in this package

In order to distinguish the different files, they are given specific prefixes, as follows:

* L points to the library, the core functions designed for the paper.
* D points to the data, which was partially re-edited, but stems mainly from Hamed and Wang (2006, produced by Feng Wang) in the form provided by List (2015) where it was converted to machine-readable text-formats which can be directly imported into LingPy.
* C points to the active scripts which carry out the analyses discussed in the paper, and which use the libraries compiled for this package.
* T points to templates which are only needed if you want to produce html plots of the detailed processes.
* P points to the reference *phylogenies*.

## Abbreviations for the models

The abbreviations used in this code package and the paper are not identical. Here is a mapping that helps with the navigation:

* pcd (partical cognate distance): this is the DWST approach described in the paper
* edit: this is the SANKOFF approach described in the paper
* simple (multistate): this is the FITCH-approach described in the paper
* simple (binary): this is the BINARY-approach described in the paper

## Results

The results can be found in files preceded by `R_acr_` (acr = ancestral state reconstruction). The following part of the filename lists first the model, and second the reference phylogeny which was used.
