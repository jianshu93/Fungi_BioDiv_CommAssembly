# Fungi_BioDiv_CommAssembly
This is code for analysis of biodiversity patterns and community assembly processes of fungal communities

Please check here for orginal publication: https://www.sciencedirect.com/science/article/pii/S1754504818303878 


RC-Daliang_linux.R is a multiprocess version of RC bray. You can change the worker=8 in the RC.p() function according to the maximum number of cores of you computer.

A quick test using 200 samples and 2500 OTUs showed that it took only a few minutes to run RC.p() function to calculate RC bray using 8 processes on my Macbook Pro.

If you have a very large dataset with a lot of samples and OTUs, it may requires a lot of memory. You can also run it on a linux server by modifiying the code a little bit with standard command line input options in R. Then you can call :


RScript RC-Daliang_linux.R input1 input2 ...



