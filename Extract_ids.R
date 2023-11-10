library(biomaRt)
library(org.Hs.eg.db)
library(dplyr)

setwd("/Users/mishrap/Kris/TARGET_RSEM")

files = list.files(pattern="results")
df= data.frame()
y= read.table(files[1], sep="\t", header=TRUE)
df <- y$gene_id
colnames(df) <- "gene_id"

for (i in files){
  print (i)
  
  x=read.table(i, sep="\t", header=TRUE)
  print(dim(x))
  slice <-x %>%  select (gene_id, FPKM)
  #slice <- slice %>%  mutate(Sample = i)
  slice<- as.data.frame(slice)
  sample_name = paste(i,"FPKM", sep="_")
  colnames(slice)[2] = sample_name
  
  df<- cbind(df, slice[2])
  
}


human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      mirror = "useast" )
#start = getLDS(attributes = "ensembl_gene_id", filters = "ensembl_gene_id", values = x$gene_id , mart = human, 
#       attributesL = "start_position", martL = human, uniqueRows=T)
listEnsemblArchives()
ensembl110 <- useEnsembl(biomart = "genes", 
                      dataset = "hsapiens_gene_ensembl", 
                      version=110 )

filters = listFilters(ensembl110)
filters[1:5,]
attributes = listAttributes(ensembl110)
attributes[1:15,]




searchAttributes(mart = ensembl110, pattern = "hgnc")

new<- getBM(attributes = c('ensembl_gene_id','external_gene_name'),
            filters = "ensembl_gene_id",
            values = x$gene_id,
            mart = ensembl110)

new3 <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
              filters = "ensembl_gene_id",
              values = df$gene_id,
              mart = ensembl110)

colnames(x)[2] ="ensembl_gene_id"
colnames(df)[1] = "ensembl_gene_id"
join_id = inner_join(df, new3, by="ensembl_gene_id")
#new_join_id =join_id[!join_id == " ",]

#join_id[join_id$hgnc_symbol == "GPC2",] %>%  
#  group_by(hgnc_symbol) %>%  
#  summarise(final_TPM = sum(TPM))
#dim(final_TPM)
#final_TPM = na.omit(final_TPM)
final_TPM <- join_id
final_TPM  <- final_TPM %>%
  group_by(hgnc_symbol, ensembl_gene_id) %>%
  summarise(across(starts_with("Sample"), list(sum =sum)))
final_TPM <- final_TPM[final_TPM$hgnc_symbol !="",]

write.csv(final_TPM, file="TARGET_merged_FPKM_fromRSEM.csv", quote=FALSE, row.names=FALSE
          )



final_TPM = join_id %>% group_by(hgnc_symbol, ensembl_gene_id) %>% summarise(final_TPM = sum(), across())
final_TPM = join_id %>% group_by(hgnc_symbol, ensembl_gene_id) %>% summarise(summary = sum(), across())


filename= "Sample_SRR7752766.isoforms.TPM.csv"
#new_final = final_TPM %>%  filter(hgnc_symbol =="")
final_TPM <- final_TPM[final_TPM$hgnc_symbol !="",]
write.csv(final_TPM, file=filename, quote=FALSE)

