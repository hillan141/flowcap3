# C1.R Flowcap 3 Challenge 1 submission

This is the top level of the FlowCap3 Challenge 1 submission.

The tarball flowcap3-challenge1-ahill.tar.gz contains code and data from the challenge submission.

To begin, unpack the tarball and enter the top level directory

tar xzvf flowcap3-challenge1-ahill.tar.gz     
cd submission-01

In this directory:
* Directory 'validate-pred-final' contains predicted class vectors for the validation set of FCS files 203.fcs - 405.fcs.  The class vectors are stored as .csv files containing no header and a class indicator (0,1,2) for each event.
* Script c1.R contains the code required to reproduce the analysis, assuming the R packages and versions listed in the R prerequisites below.


To reproduce the Challenge 1 analysis:

1.  Place the files 203.fcs to 405.fcs in a directory of your choice, say '/foo/bar/fcs'
2.  Select a directory (it does not have to exist) where the prediction .csv files will go, say '/myresults'
3.  Unzip the submission zip file, creating top-level directory "submission-01"
4.  Move to the "submission-01" directory, start R and enter the commands:

~~~~
    source("c1.R")
    fcs.root <- "/foo/bar/fcs"
    out.root <- "/myresults"
    ret <- wrap.final.wf1.flowSet(out.root=out.root)
~~~~

When the run completes, the class predictions will be found in '/myresults', in the form of one csv file per test FCS file.

Test peformance on a hold-back partition of the training set (n=101 .fcs files) was:

Class 0: 0.9998 sensitivity; 0.00034 false-positive rate
Class 1: 0.373 sensitivity; 0.64 false-positive rate
Class 2: 0.395 sensitivity; 0.34 false-positive rate

Andrew Hill
22 Oct 2012

R session info
--------------

The required R packages are flowCore and flowViz.

The test environment was as shown below.

R version 2.15.1 (2012-06-22)
Platform: x86_64-unknown-linux-gnu (64-bit)

locale:
[1] C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] flowViz_1.22.0     lattice_0.20-10    flowCore_1.24.0    rrcov_1.3-02      
[5] pcaPP_1.9-48       mvtnorm_0.9-9992   robustbase_0.9-4   Biobase_2.18.0    
[9] BiocGenerics_0.4.0

loaded via a namespace (and not attached):
 [1] IDPmisc_1.1.16      KernSmooth_2.23-8   MASS_7.3-22        
 [4] RColorBrewer_1.0-5  feature_1.2.8       graph_1.36.0       
 [7] grid_2.15.1         hexbin_1.26.0       ks_1.8.10          
[10] latticeExtra_0.6-24 stats4_2.15.1       tools_2.15.1       


