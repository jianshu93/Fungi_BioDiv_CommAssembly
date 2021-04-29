# Fungi_BioDiv_CommAssembly
This is code for analysis of biodiversity patterns and community assembly processes of fungal communities. It can also be used for bacterial community data or protist community data as long as you can have a sample by taxa table.

RC_bray_multip.r is for calculating Raup Crick metric using Bray-Curtis distance metric. A null model proposed by Stegen. et.al., 2011 is used but with abundance of each taxa considered.

Please check here for orginal publication: https://www.sciencedirect.com/science/article/pii/S1754504818303878 


RC_bray_multip.R is a multiprocess version of RC bray, you can run it like this:

```
git clone https://github.com/jianshu93/Fungi_BioDiv_CommAssembly
cd Fungi_BioDiv_CommAssembly
Rscript ./RC_bray_multip.R --input_file=./demo_input/all_samples_SV-counts.tsv --output_file=./RC_bray.txt --processors=8
```


# Tesing
A quick test using 200 samples and 2500 OTUs showed that it took only a few minutes to run RC.p() function to calculate RC bray using 8 processes on my Macbook Pro.

If you have a very large dataset with a lot of samples and OTUs, it may requires a lot of memory. You can also run it on a linux server.

# Reference
Dini-Andreote, Francisco, James C. Stegen, Jan D. Van Elsas, and Joana F. Salles. 2015. “Disentangling Mechanisms That Mediate the Balance Between Stochastic and Deterministic Processes in Microbial Succession..” Proceedings of the National Academy of Sciences of the United States of America 112(11):E1326–32.

Stegen, James C., Xueju Lin, Allan E. Konopka, and James K. Fredrickson. 2012. “Stochastic and Deterministic Assembly Processes in Subsurface Microbial Communities.” The ISME Journal 6(9):1653–64.

Stegen, James C. et al. 2012. “Stochastic and Deterministic Drivers of Spatial and Temporal Turnover in Breeding Bird Communities.” Global Ecology and Biogeography 22(2):202–12.

Zhao, Jianshu et al. 2019. “The Scale Dependence of Fungal Community Distribution in Paddy Soil Driven by Stochastic and Deterministic Processes.” Fungal Ecology 42:100856.







