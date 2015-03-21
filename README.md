# BackSPIN algorithm

An implementation of the BackSPIN biclustering algorithm as described in Zeisel et al. *Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq* **Science** 2015 (PMID: 25700174, doi: 10.1126/science.aaa1934). Please cite this paper if you use this algorithm in your work.

Original implementation by Amit Zeisel. This repo contains a standalone command-line version of BackSPIN, implemented in Python by Gioele La Manno. 

## Requirements

This version requires [pandas](http://pandas.pydata.org) and [numpy](http://www.numpy.org).

## Synopsis

       -i [inputfile]
       --input=[inputfile]
              Path of the tab delimited file.
              Rows should be genes and columns single cells/samples
              
       -o [outputfolder]
       --output=[outputfolder]
              The name of the folder where the output will be written (all the output will be 
              enclosed in a folder named runout_ddmmyyhhmmss)
              
       -d [int]
              Depth/Number of levels: The number of nested splits that will be tried by the algorythm
       -t [int]
              Number of the iterations used in the preparatory SPIN.
              Defaults to 10
       -s [float]
              Controls the decrease rate of the wid parameter used in the preparatory SPIN.
              Smaller values will increase the number of SPIN iterations and result in higher 
              precision in the first step but longer execution time.
              Defaults to 0.05
       -T [int]
              Number of the iterations used for every wid parameter.
              Does not apply on the first run (use -t instead)
              Defaults to 8
       -S [float]
              Controls the decrease rate of the wid parameter.
              Smaller values will increase the number of SPIN iterations and result in higher 
              precision but longer execution time.
              Does not apply on the first run (use -s isntead)
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
              Mean Threshold for a gene to be considered low.
              Defaults to 0.2
       -b [[axisvalue]]
              Run normal SPIN instead of backSPIN.
              Normal spin accepts the parameters -T -S
              optionally one can pass an axis value 0 to only sort genes (rows), 1 to only sort cells (columns)
       -v  
              Verbose. Print extra details of what is happening to the stdoutput 


