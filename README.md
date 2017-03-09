# BackSPIN algorithm

The BackSPIN biclustering algorithm was developed by Amit Zeisel and is described in Zeisel et al. *Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq* **Science** 2015 (PMID: [25700174](http://www.ncbi.nlm.nih.gov/pubmed/25700174), doi: [10.1126/science.aaa1934](http://dx.doi.org/10.1126/science.aaa1934)). Please cite this paper if you use the BackSPIN algorithm in your work.

Original MATLAB implementation by Amit Zeisel. This repo contains a standalone command-line version of BackSPIN, implemented in Python by Gioele La Manno. 

## Getting started

1. [Suggested, not mandatory] Install Anaconda (https://www.continuum.io/downloads) to have all the dependencies up and running (both Python 2.7/3.6 version will work)

2. Install `backspinpy` in one of the following ways:

    [Suggested] Install using **pip**:

    ```
    pip install backspinpy
    ```

    or install from source:

    ```
    git clone https://github.com/linnarsson-lab/BackSPIN
    cd BackSPIN
    python setup.py install
    ```


3. (a) Run directly your the command line tool

    ```
    backspin -i oligos.cef -o oligos_clustered.cef -f 500 -v
    ```

    (b) OR in alternative use the functions form the library directly in Python/IPython

    ```
    python
    >>> from backspinpy import SPIN, backSPIN, fit_CV, feature_selection, CEF_obj
    ```


BackSPIN takes input in CEF format and produces an annotated CEF file as output. Use [ceftools](https://github.com/linnarsson-lab/ceftools) to create and manipulate CEF files.


## Synopsis

       -i [inputfile]
       --input=[inputfile]
              Path of the cef formatted tab delimited file.
              Rows should be genes and columns single cells/samples.
              For further information on the cef format visit:
              https://github.com/linnarsson-lab/ceftools

       -o [outputfile]
       --output=[outputfile]
              The name of the file to which the output will be written

       -d [int]
              Depth/Number of levels: The number of nested splits that will be tried by the algorithm
       -t [int]
              Number of the iterations used in the preparatory SPIN.
              Defaults to 10
       -f [int]   
              Feature selection is performed before BackSPIN. Argument controls how many genes are seleceted.
              Selection is based on expected noise (a curve fit to the CV-vs-mean plot).
       -s [float]
              Controls the decrease rate of the width parameter used in the preparatory SPIN.
              Smaller values will increase the number of SPIN iterations and result in higher 
              precision in the first step but longer execution time.
              Defaults to 0.05
       -T [int]
              Number of the iterations used for every width parameter.
              Does not apply on the first run (use -t instead)
              Defaults to 8
       -S [float]
              Controls the decrease rate of the width parameter.
              Smaller values will increase the number of SPIN iterations and result in higher 
              precision but longer execution time.
              Does not apply on the first run (use -s instead)
              Defaults to 0.25
       -g [int]
              Minimal number of genes that a group must contain for splitting to be allowed.
              Defaults to 2
       -c [int]
              Minimal number of cells that a group must contain for splitting to be allowed.
              Defaults to 2
       -k [float]
              Minimum score that a breaking point has to reach to be suitable for splitting.
              Defaults to 1.15
       -r [float]
              If the difference between the average expression of two groups is lower than threshold the algorythm 
              uses higly correlated genes to assign the gene to one of the two groups
              Defaults to 0.2
       -b [axisvalue]
              Run normal SPIN instead of backSPIN.
              Normal spin accepts the parameters -T -S
              An axis value 0 to only sort genes (rows), 1 to only sort cells (columns) or 'both' for both
              must be passed
       -v  
              Verbose. Print  to the stdoutput extra details of what is happening

## Tutorial

This tutorial assumes that you have installed BackSPIN and downloaded the sample dataset `oligos.cef` (oligodendrocytes from Zeisel et al., Science 2015) from the [release page](https://github.com/linnarsson-lab/BackSPIN/releases).

Run the following from your terminal (here using Mac binary release; for other platforms, invoke backSPIN.py):

       backspin -i oligos.cef -o oligos_clustered.cef -f 500 -v -d 4

Three options are mandatory:

* The input file must be given: `-i oligos.cef`
* The ouput file must be named: `-o oligos_clustered.cef`
* The depth of clustering must be specified: `-d 4`

The **depth of clustering** `-d` is the number of levels of binary splits that will be attempted. `-d 4` indicates that BackSPIN will attempt four levels of splits, e.g. a maximum of 2<sup>4</sup> = 16 clusters will be created. The actual number of clusters may be smaller than this, because BackSPIN has a stopping rule where it will refuse to split the data further.

The **feature selection** option `-f` is not mandatory, but in practice is almost always used (not supported in the binary release, see below). This option will select a number of features (i.e. genes) based on *expected noise*. That is, genes will be ranked by how large their CV (standard deviation divided by the mean) is, compared to other genes that have similar mean expression (as described in Zeisel et al., Science 2015 and in Islam et al. Nature Methods 2014). `-f 500` will select the 500 most variable genes according to this ranking. BackSPIN runs take O(n<sup>3</sup>), so selecting more genes will quickly lead to long runs. 

For more information about feature selection, check out the [tutorial](tutorial_fselection.md).

The `-v` option makes BackSPIN print a verbose description of what's going on.

### The output

The result of BackSPIN clustering is another CEF file, which contains the same data as the input file but reorganized, as follows. First, only genes that were selected by the `-f` option are included. All other genes are removed. Second, rows and columns are sorted in SPIN order and grouped by BackSPIN cluster. Third, several row and column attributes are added, that indicate the cluster number of each gene (row) and cell (column).

       Level_0_group (always 0)
       Level_1_group (0 - 1) 
       Level_2_group (0 - 3)
       Level_3_group (0 - 7)
       Level_4_group (0 - 15)

Because BackSPIN is a biclustering method, each cluster contains both a set of genes (rows) and a set of cells (columns). At each level, the set of rows and columns that have the same cluster ID form a rectangle in the main matrix. These rectangles are non-overlapping, and every row and every column belong to exactly one cluster at each level.

You can use ceftools to manipulate the output. For example, to extract the genes that belong to cluster 0 at level 2, issue the following command:

       < oligos_clustered.cef cef select --where "Level_2_group=0"

An easy way to see the clustering result is to browse the file interactively:

       < oligos_clustered.cef cef view











