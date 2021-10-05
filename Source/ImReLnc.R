#ImReLnc is used to identify immune-related lncRNAs

ImReLnc <- function(matrix_path="../Data/matrix/",cancer_index=c(1)){
	
	
	library(stringr)
	library(org.Hs.eg.db)
	library(edgeR)
	library(fgsea)
	library(estimate)

	cancer_names=c( # 33 cancers from TCGA project
	"BRCA",
	"GBM",
	"OV",
	"LUAD",
	"UCEC",
	"KIRC",
	"HNSC",
	"LGG",
	"THCA",
	"LUSC",
	"PRAD",
	"SKCM",
	"COAD",
	"STAD",
	"BLCA",
	"LIHC",
	"CESC",
	"KIRP",
	"SARC",
	"LAML",
	"PAAD",
	"ESCA",
	"PCPG",
	"READ",
	"TGCT",
	"THYM",
	"KICH",
	"ACC",
	"MESO",
	"UVM",
	"DLBC",
	"UCS",
	"CHOL"
	)
	if(!(file.exists("../Result"))){dir.create("../Result")}
	if(!(file.exists("../Result/gct"))){dir.create("../Result/gct")} #the dictionary of gct file
	if(!(file.exists("../Result/all"))){dir.create("../Result/all")} #the dictionary of lncRNA-pathway correlations
	if(!(file.exists("../Result/lnc_pathway"))){dir.create("../Result/lnc_pathway")} #the dictionary of lncRNA-pathway pairs 
	
	
	for(i in 1:length(cancer_index))
	{
		print(cancer_names[i])
		#read mRNA expression profile
		print("Reading mRNA expression profile...");
		mrna_name=paste(matrix_path, cancer_names[i], sep = "", collapse = "")
		mrna_name=paste(mrna_name, "_mrna_matrix", sep = "", collapse = "")
		mRNA_express <- read.table(file = mrna_name, 
									sep = "\t", header = TRUE, row.names = 1,
									stringsAsFactors = TRUE)
		#mRNA_express <- subset(mRNA_express, select = -X );	

			
		#normalization of mRNA expression profile
		print("Normalization of mRNA expression profile...")
		dge <- DGEList(counts = mRNA_express)
		dge <- calcNormFactors(dge,method="TMM") #计算标准化因子
		mRNA_express <- cpm(dge, log=TRUE, prior.count=3)#标准化并log化

			
		#convert Ensemble id to EntrezID
		sample_no=ncol(mRNA_express)
		mRNA_express=as.data.frame(mRNA_express)
		ENSEMBLE_ID=rownames(mRNA_express)
		mRNA_express=cbind(mRNA_express,ENSEMBLE_ID)
		mRNA_express$ENSEMBLE_ID=unlist(str_split(mRNA_express$ENSEMBLE_ID,"[.]",simplify=T))[,1]	#Delete the value after the decimal point of Ensemble id
		k=keys(org.Hs.eg.db,keytype = "ENSEMBL")
		list=select(org.Hs.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")	
		ID_list=list[match(mRNA_express$ENSEMBLE_ID,list[,"ENSEMBL"]),]
		mRNA_express$ENTREZID=ID_list$ENTREZID
		mRNA_express=na.omit(mRNA_express)
		mRNA_express <- mRNA_express[,c(sample_no+2,1:sample_no)] 
		
		
		#delete duplicate mRNA expression item
		mRNA_express=mRNA_express[!duplicated(mRNA_express$ENTREZID), ]
		rownames(mRNA_express)=mRNA_express[,1]
		mRNA_express=mRNA_express[,-1]

				
		#compute tumor purity
		print("Computing tumor purity...")
		mrna_express_name=paste(cancer_names[i],"_mrna_express", sep = "", collapse = "")
		mrna_gct_name=paste(cancer_names[i],"_mrna_express.gct", sep = "", collapse = "")
		mrna_score_name=paste(cancer_names[i],"_mrna_express_score.gct", sep = "", collapse = "")
		
		write.table(mRNA_express, file =mrna_express_name,sep="\t",row.names =TRUE, col.names =TRUE, quote =FALSE);
		filterCommonGenes(mrna_express_name, mrna_gct_name, "EntrezID")	
		estimateScore(mrna_gct_name,mrna_score_name, platform="illumina")

		mRNA_express_score <- readLines(mrna_score_name)
		NAME <- unlist(strsplit(mRNA_express_score[grep("NAME",mRNA_express_score)],"\t"))
		NAME <- NAME[3:length(NAME)]

		ESTIMATEScore <- unlist(strsplit(mRNA_express_score[grep("ESTIMATEScore",mRNA_express_score)],"\t"))
		ESTIMATEScore <- as.numeric(ESTIMATEScore[3:length(ESTIMATEScore)])
		tumor_purity=cos(0.6049872018+0.0001467884*ESTIMATEScore)

		Tumor_purity = data.frame(tumor_purity)
		row.names(Tumor_purity)=NAME
		
		
		
		#read lncRNA expression profile
		print("Read lncRNA expression profile...")
		lnc_name=paste(matrix_path, cancer_names[i], sep = "", collapse = "")
		lnc_name=paste(lnc_name, "_lnc_symbol_matrix", sep = "", collapse = "")
		lncRNA_express <- read.table(file = lnc_name, 
									sep = "\t", header = TRUE, row.names = NULL,
									stringsAsFactors = TRUE)
		colnames(lncRNA_express)[1]="ID"
		#lncRNA_express <- subset(lncRNA_express, select = -X );	


		#delete duplicate lncRNA expression item
		lncRNA_express=as.data.frame(lncRNA_express)
		lncRNA_express=lncRNA_express[!duplicated(lncRNA_express$ID), ]
		rownames(lncRNA_express)=lncRNA_express[,1]
		lncRNA_express=lncRNA_express[,-1]


		#normalization of lncRNA expression profile
		print("Normalization of lncRNA expression profile")
		dge2 <- DGEList(counts = lncRNA_express)
		dge2 <- calcNormFactors(dge2,method="TMM") #计算标准化因子
		lncRNA_express <- cpm(dge2, log=TRUE, prior.count=3)#标准化并log化
		
		
		#coupute partial correlation
		print("Couputing partial correlation")
		temp_name=rownames(mRNA_express)
		mRNA_express=as.data.frame(lapply(mRNA_express,as.numeric))
		rownames(mRNA_express)=temp_name
		
		pear_cor_coef=cor(t(mRNA_express),t(lncRNA_express))

		
		cor_lm=pear_cor_coef
		cor_mt=cor(t(mRNA_express),Tumor_purity)
		cor_lt=cor(Tumor_purity,t(lncRNA_express))
		cor_ltm=cor_mt%*%cor_lt
		temp1=cor_lm-cor_ltm
		temp2=(sqrt(1-cor_mt*cor_mt))%*%(sqrt(1-cor_lt*cor_lt))
		par_cor_coef=temp1/temp2    #19712 14805

		for(ii in 1:nrow(par_cor_coef))
		{
			for(j in 1:ncol(par_cor_coef))
			{
				if(is.na(par_cor_coef[ii,j]))
				{
					par_cor_coef[ii,j]=runif(1, -1, 1)
				}
			}
		}


		#compute rank score		#sample number = 492
		print("Computing ranking score...")
		statistic <- par_cor_coef*sqrt((sample_no-2-1)/(1-par_cor_coef^2))   # 19712 14805
		p_value <- 2*pnorm(-abs(statistic)) # 19712 14805	pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)从-无穷到q求正态分布的积分值
		RS = -log10(p_value)*sign(p_value)  # 19016 14805
		RS=t(RS)
		RS=na.omit(RS)

		#gsea enrichment analysis	
		print("Immune-related enrichment analysis...")
		immune_pathway=read.csv(file = "../Data/pathway/GeneList.txt", 
									sep = "\t", header = TRUE, row.names = NULL,
									stringsAsFactors = TRUE)


		list1=subset(immune_pathway,Category=="Antigen_Processing_and_Presentation")[,2]
		list2=subset(immune_pathway,Category=="Antimicrobials")[,2]					
		list3=subset(immune_pathway,Category=="BCRSignalingPathway")[,2]
		list4=subset(immune_pathway,Category=="Chemokines")[,2]	
		list5=subset(immune_pathway,Category=="Chemokine_Receptors")[,2]
		list6=subset(immune_pathway,Category=="Cytokines")[,2]	
		list7=subset(immune_pathway,Category=="Cytokine_Receptors")[,2]
		list8=subset(immune_pathway,Category=="Interferons")[,2]					
		list9=subset(immune_pathway,Category=="Interferon_Receptor")[,2]
		list10=subset(immune_pathway,Category=="Interleukins")[,2]	
		list11=subset(immune_pathway,Category=="Interleukins_Receptor")[,2]
		list12=subset(immune_pathway,Category=="NaturalKiller_Cell_Cytotoxicity")[,2]	
		list13=subset(immune_pathway,Category=="TCRsignalingPathway")[,2]
		list14=subset(immune_pathway,Category=="TGFb_Family_Member")[,2]					
		list15=subset(immune_pathway,Category=="TGFb_Family_Member_Receptor")[,2]
		list16=subset(immune_pathway,Category=="TNF_Family_Members")[,2]	
		list17=subset(immune_pathway,Category=="TNF_Family_Members_Receptors")[,2]

							
		pathway_list=list(Antigen_Processing_and_Presentation=list1,
		Antimicrobials=list2,
		BCRSignalingPathway=list3,
		Chemokines=list4,
		Chemokine_Receptors=list5,
		Cytokines=list6,
		Cytokine_Receptors=list7,
		Interferons=list8,
		Interferon_Receptor=list9,
		Interleukins=list10,
		Interleukins_Receptor=list11,
		NaturalKiller_Cell_Cytotoxicity=list12,
		TCRsignalingPathway=list13,
		TGFb_Family_Member=list14,
		TGFb_Family_Member_Receptor=list15,
		TNF_Family_Members=list16,
		TNF_Family_Members_Receptors=list17)

		fgseaRes_all <- c()
		if(length(RS)==0)
		{
			next()
		}
		for(ii in 1:nrow(RS))
		{
			##remove the rows contain Infs
			if(sum( is.infinite(unlist(RS[ii,]))) != 0){
			  next()
			}
			ranks <- RS[ii,]
			fgseaRes <- fgsea(pathway_list, ranks, minSize=1, maxSize=5000, nperm=1000)
			sigValue <- c()
			for(j in 1:nrow(fgseaRes)){
			  if(fgseaRes$ES[j]>0){
				sig_ij <- 1 - 2*fgseaRes$pval[j]
			  }else{
				sig_ij <- 2*fgseaRes$pval[j] - 1
			  }
			  sigValue <- c(sigValue,sig_ij)
			  
			}
			lncRNA <- rownames(RS)[ii]
			fgseaRes_i <- cbind(lncRNA,fgseaRes,sigValue)
			fgseaRes_all <- rbind(fgseaRes_all,fgseaRes_i)
		}
		k=0.995
		sig_ind <- which(abs(fgseaRes_all$sigValue) >= k) #sigValue文中定义的为显著性指标（lncRES(i,k)）,k=0.995
		sig_pairs <- fgseaRes_all[sig_ind,1:2]	#select lncRNA, pathway(k=0.995) sig_pairs即
		sig_pairs <- as.matrix(sig_pairs)
		gsea.Res <- list(sig_pairs,fgseaRes_all)
		names(gsea.Res) <- c("sig_pairs","fgseaRes_all")

		#Output the result
		print("Outputing the results...")
		all_file_name=paste("../Result/all/", cancer_names[i], sep = "", collapse = "")
		all_file_name=paste(all_file_name,"_all_lnc_path_pairs", sep = "", collapse = "")
		lnc_path_file_name=paste("../Result/lnc_pathway/", cancer_names[i], sep = "", collapse = "")
		lnc_path_file_name=paste(lnc_path_file_name,"_immune_lnc_path_pairs", sep = "", collapse = "")
		write.table(as.matrix(fgseaRes_all), file =all_file_name,sep="\t",row.names =TRUE, col.names =TRUE, quote =TRUE);
		write.table(sig_pairs, file =lnc_path_file_name,sep="\t",row.names =TRUE, col.names =TRUE, quote =TRUE);
		
		all_gct_name=paste("mv ","*.gct", sep = "", collapse = "")
		all_gct_name=paste(all_gct_name, " ../Result/gct/", sep = "",collapse = "")
		system(all_gct_name)	

		all_expr_name=paste("rm ","*_mrna_express", sep = "", collapse = "")
		system(all_expr_name)
	}

	
}

PanCancer<- function(lnc_pathway_path="../Result/lnc_pathway/",lnc_exp_path="../Data/matrix/",gama=c(0.4,0.3,0.3)){
	#Determine whether all paths exist
	if (!file.exists(lnc_pathway_path)){
		print("lnc_pathway_path does not exist")
	}
	
	if (!file.exists(lnc_exp_path)){
		print("lnc_exp_path does not exist")
	}


	library(xlsx)
	library(edgeR)

	cancer_names=c(
	"BRCA",
	"GBM",
	"OV",
	"LUAD",
	"UCEC",
	"KIRC",
	"HNSC",
	"LGG",
	"THCA",
	"LUSC",
	"PRAD",
	"SKCM",
	"COAD",
	"STAD",
	"BLCA",
	"LIHC",
	"CESC",
	"KIRP",
	"SARC",
	"LAML",
	"PAAD",
	"ESCA",
	"PCPG",
	"READ",
	"TGCT",
	"THYM",
	"KICH",
	"ACC",
	"MESO",
	"UVM",
	"DLBC",
	"UCS",
	"CHOL"
	)

	immune_names=c(
	"Antigen_Processing_and_Presentation",
	"Antimicrobials",					
	"BCRSignalingPathway",
	"Chemokines",	
	"Chemokine_Receptors",
	"Cytokines",	
	"Cytokine_Receptors",
	"Interferons",					
	"Interferon_Receptor",
	"Interleukins",	
	"Interleukins_Receptor",
	"NaturalKiller_Cell_Cytotoxicity",	
	"TCRsignalingPathway",
	"TGFb_Family_Member",					
	"TGFb_Family_Member_Receptor",
	"TNF_Family_Members",	
	"TNF_Family_Members_Receptors"
	)
	
	cancer_matrix=c()
	de_matrix=c()


	#read lnc-pathway pairs
	print("read lnc-pathway pairs")
	for(i in 1:length(cancer_names))
	{
		lnc_path_name=paste(lnc_pathway_path, cancer_names[i], sep = "", collapse = "")
		lnc_path_name=paste(lnc_path_name, "_immune_lnc_path_pairs", sep = "", collapse = "")
		
		if (file.exists(lnc_path_name))
		{
			cancer_pair=read.table(file = lnc_path_name, 
					sep = "\t", header = TRUE, 
					row.names = 1, stringsAsFactors = FALSE);
			if(nrow(cancer_pair)>0)
			{
				cancer_temp=cbind(cancer_pair,rep(cancer_names[i],nrow(cancer_pair)))
				colnames(cancer_temp)[3]="cancer_name"
				cancer_matrix <- rbind(cancer_matrix,cancer_temp)
			}
			
		}
		
	}

	#find all immune-related lncRNA
	immune_lnc=cancer_matrix[,1]
	immune_lnc <- immune_lnc[!duplicated(immune_lnc)]

	#read lncRNA expression data and conduct DE analysis
	print("read lncRNA expression data and conduct DE analysis")
	for(i in 1:length(cancer_names))
	{
		lnc_exp_name=paste(lnc_exp_path, cancer_names[i], sep = "", collapse = "")
		lnc_exp_name=paste(lnc_exp_name, "_lnc_symbol_matrix", sep = "", collapse = "")
		
		if (file.exists(lnc_exp_name))
		{
			#mRNA_express[!duplicated(mRNA_express$ENTREZID), ]
			lncRNA_express <- read.table(file = lnc_exp_name, 
									sep = "\t", header = TRUE, row.names = NULL,
									stringsAsFactors = TRUE)
			colnames(lncRNA_express)[1]="ID"
			#lncRNA_express <- subset(lncRNA_express, select = -X );	
			


			#delete duplicate lncRNA expression item and filter disease-related lncRNA
			lncRNA_express=as.data.frame(lncRNA_express)
			lncRNA_express=lncRNA_express[!duplicated(lncRNA_express$ID), ]
			rownames(lncRNA_express)=lncRNA_express[,1]
			lncRNA_express=lncRNA_express[,-1]
			lncRNA_express=lncRNA_express[immune_lnc,]
			lncRNA_express=na.omit(lncRNA_express)
			#DE analysis
			lncRNA_express=lncRNA_express[rowSums(cpm(lncRNA_express)>1) >= 2,]	#filter low expression data 9447*488

			group_list=ifelse(as.numeric(substr(colnames(lncRNA_express),14,15)) < 10,'tumor','normal');	#compute the group information of samples
			#table(group_list);	#
			de_result=NULL;
			if((length(unique(group_list))>1)&&(nrow(lncRNA_express)>0))
			{
				lnc_model <- DGEList(counts=lncRNA_express,group=factor(group_list))	#构建模型，只需要两个因素
				TMM_lnc_express <- calcNormFactors(lnc_model)		#compute norm factor

				design <- model.matrix(~0+factor(group_list))	#把group设置成一个model matrix
				rownames(design)<-colnames(TMM_lnc_express)
				colnames(design)<-levels(factor(group_list))

				dge <- estimateGLMCommonDisp(TMM_lnc_express,design)
				dge <- estimateGLMTrendedDisp(dge, design)
				dge <- estimateGLMTagwiseDisp(dge, design)
				fit <- glmFit(dge, design)  #Fit glm
				lrt <- glmLRT(fit,  contrast=c(-1,1)) #因为是两组所以这么做，按照说明书
				nrDEG=topTags(lrt, n=nrow(dge),adjust.method = "BH",sort.by = "PValue")	#2696
				nrDEG=as.data.frame(nrDEG)

				if( nrow(nrDEG)>0)
				{
					de_result=nrDEG[abs(nrDEG$logFC) > 1,];
					de_result=de_result[abs(de_result$FDR) <= 0.05,];
					de_temp=cbind(rownames(de_result),de_result$logFC,rep(cancer_names[i],nrow(de_result)))
					colnames(de_temp)=c("lncRNA","logFC","cancer_name")
					de_matrix=rbind(de_matrix,de_temp)
				}	
			}
		}
		
	}


	#compute F(rank)
	Frank=matrix(0,length(immune_lnc)*length(immune_names),4)
	colnames(Frank)=c("irlncRNA","pathway","cancer_names","norm_cancer_count")

	#compute the number of immune-lnc pairs in cancers
	print("compute the number of immune-lnc pairs in cancers")
	for(i in 1:length(immune_lnc))
	{
		for(j in 1:length(immune_names))
		{
			for(k in 1:nrow(cancer_matrix))
			{
				i_no=(i-1)*17+j
				Frank[i_no,1]=immune_lnc[i]
				Frank[i_no,2]=immune_names[j]
				if((cancer_matrix[k,1]==immune_lnc[i])&&(cancer_matrix[k,2]==immune_names[j]))
				{
						
					if(Frank[i_no,3]==0)
					{
						Frank[i_no,3]=paste("", cancer_matrix[k,3], sep = "", collapse = "")

					}
					else
					{
						Frank[i_no,3]=paste(Frank[i_no,3], ",", sep = "", collapse = "")
						Frank[i_no,3]=paste(Frank[i_no,3], cancer_matrix[k,3], sep = "", collapse = "")

					}
					Frank[i_no,4]=as.numeric(Frank[i_no,4])+1
				}

			}
		}
	}
	print("normalization of rank i")
	Frank[,4]=(as.numeric(Frank[,4])-min(as.numeric(Frank[,4])))/
	(max(as.numeric(Frank[,4]))-min(as.numeric(Frank[,4])))	#normalization of rank i

	#compute Deg value in cancers
	print("coupute Deg value in cancers")
	Deg=matrix(0,length(immune_lnc),4)
	colnames(Deg)=c("irlncRNA","cancer_names","norm_cancer_count","norm_aver_logFC")
	for(i in 1:length(immune_lnc))
	{
		for(j in 1:nrow(de_matrix))
		{
			Deg[i,1]=immune_lnc[i]
			if(immune_lnc[i]==de_matrix[j,1])
			{	if(Deg[i,2]==0)
				{
					Deg[i,2]=paste("", de_matrix[j,3], sep = "", collapse = "")
				}
				else
				{
					Deg[i,2]=paste(Deg[i,2], ",", sep = "", collapse = "")
					Deg[i,2]=paste(Deg[i,2], de_matrix[j,3], sep = "", collapse = "")
				}
				Deg[i,3]=as.numeric(Deg[i,3])+1	
				Deg[i,4]=as.numeric(Deg[i,4])+as.numeric(de_matrix[j,2])	
			}
		}
	}

	for(i in 1:length(immune_lnc))
	{
		if(Deg[i,3]==0)
		{
			Deg[i,4]=0
		}
		else
		{
			Deg[i,4]=as.numeric(Deg[i,4])/as.numeric(Deg[i,3])
		}
	}

	print("normalization of deg and logFC")	
	Deg[,3]=(as.numeric(Deg[,3])-min(as.numeric(Deg[,3])))/
	(max(as.numeric(Deg[,3]))-min(as.numeric(Deg[,3])))	#normalization of deg

	Deg[,4]=(as.numeric(Deg[,4])-min(as.numeric(Deg[,4])))/
	(max(as.numeric(Deg[,4]))-min(as.numeric(Deg[,4])))	#normalization of deg



	write.table(Deg, file ="../Result/Deg_norm.txt",sep="\t",
			row.names =FALSE, col.names =TRUE, quote =TRUE)
	write.table(Frank, file ="../Result/NC_norm.txt",sep="\t",
			row.names =FALSE, col.names =TRUE, quote =TRUE)



	#"coupute final Fscore"
	print("coupute final Fscore")	
	Frank2=matrix(0,length(immune_lnc),2)
	colnames(Frank2)=c("irlncRNA","Pscore")


	#出现pathway出现次数所占比例，进行初始化	gama=c(0.4,0.3,0.3)
	for(i in 1:length(immune_lnc))
	{
		Frank2[i,1]=immune_lnc[i]
		for(j in 1:length(immune_names))
		{
			Frank2[i,2]=as.numeric(Frank2[i,2])+as.numeric(Frank[(i-1)*17+j,4])
		}
		Frank2[i,2]=as.numeric(Frank2[i,2])/17
		Frank2[i,2]=gama[1]*as.numeric(Frank2[i,2])+gama[2]*as.numeric(Deg[i,3])+gama[3]*as.numeric(Deg[i,4])
		
	}
	Frank2=Frank2[order(Frank2[,2],decreasing=T),]
	write.table(Frank2, file ="../Result/Pscore.txt",sep="\t",
			row.names =FALSE, col.names =TRUE, quote =TRUE)
			
}