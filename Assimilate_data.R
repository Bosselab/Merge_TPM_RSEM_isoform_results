setwd("~/Kris/Esthe")
library(dplyr)
library(biomaRt)
library(org.Hs.eg.db)




combine_files_find_genes <- function(path_name,ensembl110){
  list_files = list.files(path_name, pattern ="isoform")
  df=list()
  for (i in list_files){
    #i <- paste(path_name,i, sep="/")
    file = read.table(i, header=TRUE, sep="\t", check.names=FALSE)
    print (i)
    file<- file %>% mutate(File_name=i)
    print (dim(file))
    df<- rbind(df, file)
  }
  
  all_trans <- unique(df$transcript_id)
  
  x=read.table(list_files[1], sep="\t", header=TRUE)
  
  
  new <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
               filters = "ensembl_gene_id",
               values = x$gene_id,
               mart = ensembl110)
  
  colnames(new)[1] ="gene_id"
  
  df<-df %>% inner_join(new, by=join_by("gene_id"), relationship = "many-to-many")
  df<- df %>% group_by(File_name, hgnc_symbol) %>% summarise(total=sum(TPM), across())
  GPC2<- df[df$hgnc_symbol =="GPC2",] %>% group_by(File_name) %>% summarise(total=sum(TPM), across())
  GPC2_final <- GPC2 %>% group_by(File_name) %>% distinct (total, File_name) %>% select(total, File_name)
  GPC2_target = paste("TARGET_", i, sep="")
  write.table(GPC2_final, file=GPC2_target, quote=FALSE)
   #return (df)
  
}





# Enseml data
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

#start = getLDS(attributes = "ensembl_gene_id", filters = "ensembl_gene_id", values = x$gene_id , mart = human, 
#       attributesL = "start_position", martL = human, uniqueRows=T)
listEnsemblArchives()
ensembl110 <- useEnsembl(biomart = "genes", 
                         dataset = "hsapiens_gene_ensembl", 
                         version=110 )



df1 <- combine_files_find_genes(".", ensembl110)
new_df1<-df1[df1$hgnc_symbol != "",]
#head(Esthe_GPC2)
#TARGET_GPC2 <-combine_files_find_genes(".", ensembl110)
#head(TARGET_GPC2)
colnames(new_df1)[3] ="Final_TPM"
write_xlsx(split(new_df1, new_df1$File_name), "TARGET_merged_TPM_fromRSEM.xlsx")

  