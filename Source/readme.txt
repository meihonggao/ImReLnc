==============================================================================================================================================================
ImReLnc: identifying immune-related lncRNA characteristics in human cancers based on heuristic correlation optimization (contact: meihonggao@mail.nwpu.edu.cn)
==============================================================================================================================================================


Contents
--------
1.Folder information
2.Function
3.Output
4.Demo
5.Platform

1.Folder information
--------------------
There are 3 folders in ImReLnc and these are:
Source: This folder contains R program and readme file
Data: This folder contains deault datas
Result: This folder contains the result datas


2.Functions
--------------------
2.1 ImReLnc 
parameters:
matrix_path(Optional Default="../Data/matrix/")	##The dictionary of express matrix for 33 cancers
cancer_index(Optional Default=c(1))	##The indexes of cancers in which immune-related lncRNAs will be identified, and the index 1-33 represent BRCA,GBM,OV,LUAD,UCEC,KIRC,HNSC,LGG,THCA,LUSC,PRAD,SKCM,COAD,STAD,BLCA,LIHC,CESC,KIRP,SARC,LAML,PAAD,ESCA,PCPG,READ,TGCT,THYM,KICH,ACC,MESO,UVM,DLBC,UCS, and CHOL respectively.
example:
ImReLnc(,c(1,5))

2.2 PanCancer
parameters:
lnc_pathway_path(Optional Default="../Result/lnc_pathway/")
lnc_exp_path(Optional Default="../Data/matrix/")
gama(Optional Default=c(0.4,0.3,0.3))
example:
PanCancer()


3.Output
--------
3.1 lncRNA pathway pair is in ".../ImReLnc/Result/lnc_pathway/" folder 		##Identification results of immune-related lncRNA in each cancer 
3.2 Pscore.txt is in "/ImReLnc/Result/" folder 			##Pathogenicity level of immune-related lncRNA 


4.Demo
------
Demo.R in the Source folder is an example of using the ImReLnc program 


5.Platform
----------
Software: R x64 4.0.2
Dependencies:stringr,org.Hs.eg.db,edgeR,fgsea,estimate,xlsx,edgeR


    

