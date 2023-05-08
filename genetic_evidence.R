###COMBINING INFORMATION FROM OPENTARGETS
setwd("~/work/Opentargets/")

library(dplyr)
library(sparklyr)
library(sparklyr.nested)
library(parqr)

#Opentargets data version 22.06 [28-Jun-2022]

sc <- spark_connect(master = "local")

## read drug dataset
drugPath="molecule"
drug=spark_read_parquet(sc, path = drugPath)

## Browse the drug schema to define variables of interest 
columns_drugs <- drug %>%
  sdf_schema() %>%
  lapply(function(x) do.call(tibble, x)) %>%
  bind_rows()

## Select variables of interest
drugSelect <- drug %>%
  select(id,
         drugType,
         isApproved,
         yearOfFirstApproval,
         name,
         tradeNames,
         parentId,
         linkedTargets, 
         linkedTargets) %>% sdf_unnest(linkedTargets,keep_all = TRUE)

drugSelect_one<- drugSelect %>%
  select(id,
         parentId,
         drugType,
         yearOfFirstApproval,
         name,
         tradeNames,
         isApproved,
         rows) %>% sdf_explode(rows,keep_all = TRUE)

drugSelect_two=drugSelect_one%>% collect()
duplicates=drugSelect_two[duplicated(drugSelect_two), ]; nrow(duplicates) #check in case there are duplicated rows 
head(drugSelect_two);nrow(drugSelect_two)
names(drugSelect_two)[[8]]="ensembl_gene_id"
nrow(drugSelect_two[unique(drugSelect_two$id),]) #12854
nrow(drugSelect_two[unique(drugSelect_two$ensembl_gene_id),]) #1463
nrow(drugSelect_two[is.na(drugSelect_two$ensembl_gene_id),]) #8828 have missing target 
#Keep only parent CHEMLID. the child CHEMLID is a modified version of the same active ingredient 
drugSelect_two=subset(drugSelect_two,(is.na(drugSelect_two$parentId)))
nrow(drugSelect_two[unique(drugSelect_two$id),])

drugSelect_two=subset(drugSelect_two,!is.na(drugSelect_two$isApproved))
nrow(drugSelect_two[unique(drugSelect_two$id),])

investigational=subset(drugSelect_two,drugSelect_two$isApproved==FALSE)
nrow(investigational[unique(investigational$id),])

approved=subset(drugSelect_two,drugSelect_two$isApproved==TRUE)
colnames(approved)[[1]]="chemblIds" 

# Keep only drugs of interest 
cancer_remove=read.table("cancer_drugs.txt",h=F)
set_one=approved[approved$name%in%cancer_remove$name,]
set_one=approved[approved$chemblIds%in%cancer_remove$chemblIds,]
set_one=approved[!(approved$name%in%cancer_remove$name),]
set_one=set_one[!(set_one$chemblIds%in%cancer_remove$chemblIds),]

infections=subset(db_drugs_atc,db_drugs_atc$code_4=="J" |db_drugs_atc$code_4=="P") #ANTIINFECTIVES FOR SYSTEMIC USE 
set_two=set_one[set_one$name%in%infections$name,]
set_two=set_one[!(set_one$name%in%infections$name),]

hormones=subset(db_drugs_atc,db_drugs_atc$code_4=="H" |db_drugs_atc$code_3=="G03") #SYSTEMIC HORMONAL PREPARATIONS AND SEX HORMONES 
set_three=set_two[set_two$name%in%hormones$name,]
set_three=set_two[!(set_two$name%in%hormones$name),]

vitamin_minerals=subset(db_drugs_atc,db_drugs_atc$code_3=="A11" |db_drugs_atc$code_3=="A12") #minerals and vitamines
set_four=set_three[set_three$name%in%vitamin_minerals$name,]

clean_excluded_drugs=set_three[!(set_three$name%in%vitamin_minerals$name),]
nrow(clean_excluded_drugs[unique(clean_excluded_drugs$chemblIds),])

#Include ids approve after July 2021
chemblIds_additional=c("CHEMBL4297741","CHEMBL4297551","CHEMBL1631694","CHEMBL2338675","CHEMBL4298211","CHEMBL3707229","CHEMBL1743081",
                       "CHEMBL3707276","CHEMBL259571","CHEMBL4297839","CHEMBL4297832","CHEMBL4297750","CHEMBL3655081","CHEMBL4297590")

#CHEMBL1509 approved 2001 so it was not included,"CHEMBL827" is missing year

#Updated approval status and year
drugSelect_two$yearOfFirstApproval[drugSelect_two$id=="CHEMBL4297741"]=2021
drugSelect_two$isApproved[drugSelect_two$id =="CHEMBL4297741"]="TRUE"
drugSelect_two$yearOfFirstApproval[drugSelect_two$id=="CHEMBL827"]=2021
drugSelect_two$isApproved[drugSelect_two$id =="CHEMBL827"]="TRUE"
drugSelect_two$yearOfFirstApproval[drugSelect_two$id=="CHEMBL4297551"]=2021
drugSelect_two$isApproved[drugSelect_two$id =="CHEMBL4297551"]="TRUE"
drugSelect_two$yearOfFirstApproval[drugSelect_two$id=="CHEMBL3707229"]=2021
drugSelect_two$isApproved[drugSelect_two$id =="CHEMBL3707229"]="TRUE"
drugSelect_two$yearOfFirstApproval[drugSelect_two$id=="CHEMBL1743081"]=2021
drugSelect_two$isApproved[drugSelect_two$id =="CHEMBL1743081"]="TRUE"
drugSelect_two$yearOfFirstApproval[drugSelect_two$id=="CHEMBL3707276"]=2021
drugSelect_two$isApproved[drugSelect_two$id =="CHEMBL3707276"]="TRUE"
drugSelect_two$yearOfFirstApproval[drugSelect_two$id=="CHEMBL4297590"]=2022
drugSelect_two$isApproved[drugSelect_two$id =="CHEMBL4297590"]="TRUE"
drugSelect_two$yearOfFirstApproval[drugSelect_two$id=="CHEMBL259571"]=2022
drugSelect_two$isApproved [drugSelect_two$id =="CHEMBL259571"]="TRUE"
drugSelect_two$yearOfFirstApproval[drugSelect_two$id=="CHEMBL4297839"]=2022
drugSelect_two$isApproved [drugSelect_two$id =="CHEMBL4297839"]="TRUE"
drugSelect_two$yearOfFirstApproval[drugSelect_two$id=="CHEMBL4297832"]=2022
drugSelect_two$isApproved [drugSelect_two$id =="CHEMBL4297832"]="TRUE"
drugSelect_two$yearOfFirstApproval[drugSelect_two$id=="CHEMBL4297750"]=2022
drugSelect_two$isApproved [drugSelect_two$id =="CHEMBL4297750"]="TRUE"
drugSelect_two$yearOfFirstApproval[drugSelect_two$id=="CHEMBL3655081"]=2022
drugSelect_two$isApproved [drugSelect_two$id =="CHEMBL3655081"]="TRUE"

#Additional missing IDs
MITAPIVAT=c("CHEMBL4299940","Protein","2022","MITAPIVAT"," ","TRUE","ENSG00000143627")
MAVACAMTEN=c("CHEMBL4297517","Small molecule","2022","MAVACAMTEN"," ","TRUE","ENSG00000092054")
VUTRISIRAN=c("CHEMBL4594511","Oligonucleotide","2022","VUTRISIRAN"," ","TRUE","ENSG00000118271")
#CHEMBL1568698 has 10 targets but only two are linked to approved indication [ENSG00000186297,ENSG00000187730] 
GANAXOLONE_one=c("CHEMBL1568698","Small molecule","2022","GANAXOLONE"," ","TRUE","ENSG00000186297")
GANAXOLONE_two=c("CHEMBL1568698","Small molecule","2022","GANAXOLONE"," ","TRUE","ENSG00000187730")

new_approved=drugSelect_two[(drugSelect_two$id%in%chemblIds_additional),]

names(new_approved)[[1]]="chemblIds"
approved_added_new=rbind(clean_excluded_drugs,new_approved,MITAPIVAT,MAVACAMTEN,VUTRISIRAN,GANAXOLONE_one,GANAXOLONE_two)
approved_added_new$id=paste(approved_added_new$chemblIds,approved_added_new$ensembl_gene_id,sep="_")
approved_target=subset(approved_added_new,!(is.na(approved_added_new$ensembl_gene_id)))
nrow(approved_target[unique(approved_target$chemblIds),])
approved_no_targets=subset(approved_added_new,(is.na(approved_added_new$ensembl_gene_id)))
nrow(approved_no_targets[unique(approved_no_targets$chemblIds),])

ERT=c("CHEMBL3990026","CHEMBL2108888","CHEMBL4594320","CHEMBL1201824","CHEMBL2108311","CHEMBL3544921","CHEMBL4297549","CHEMBL1201826","CHEMBL1964120",
      "CHEMBL1201865","CHEMBL1201632","CHEMBL1201595","CHEMBL1201514","CHEMBL4297802","CHEMBL3039537","CHEMBL3707382")
# CHEMBL4298010 manual input

ERT_df=approved_added_new[(approved_added_new$chemblIds%in%ERT),]
ERT_df$id=paste(ERT_df$chemblIds,ERT_df$ensembl_gene_id,sep="_")

## read MoA dataset
MoAPath="mechanismOfAction"
MoA=spark_read_parquet(sc,path = MoAPath)

columns_MoA <- MoA %>%
  sdf_schema() %>%
  lapply(function(x) do.call(tibble, x)) %>%
  bind_rows()

MoA_Set<- MoA %>%
  select(actionType,
         mechanismOfAction,
         chemblIds,
         targetName,
         targetType,
         targets) %>% sdf_explode(chemblIds,keep_all = TRUE)
MoA_Set<- MoA_Set %>% sdf_explode(targets,keep_all = TRUE)
MoA_df=MoA_Set%>% collect()
MoA_df=MoA_df[!duplicated(MoA_df), ]
head(MoA_df)
MoA_df$id=paste(MoA_df$chemblIds,MoA_df$targets,sep="_")
MoA_target=subset(MoA_df,!(is.na(MoA_df$targets)))

#Add MoA to all drugs with gene targets 
approved_MoA=merge(approved_target,MoA_target,by="id",all.x = T)
ERT_MoA=MoA_df[(MoA_df$chemblIds%in%ERT),]
ERT_MoA$id=paste(ERT_MoA$chemblIds,ERT_MoA$targets,sep="_")
ERT_final=merge(ERT_df,ERT_MoA,by="id")
approved_MoA_ready=rbind(approved_MoA,ERT_final)
approved_MoA_ready$rows=NULL
approved_MoA_ready$isApproved=NULL
approved_MoA_ready$chemblIds.y=NULL
names(approved_MoA_ready)[[2]]="chemblIds"

## read indication dataset
indicationPath="indication"

indication<- spark_read_parquet(sc,path = indicationPath)
columns_indication <- indication %>%
  sdf_schema() %>%
  lapply(function(x) do.call(tibble, x)) %>%  
  bind_rows()

indicationSelect= indication %>% sdf_unnest(indications,keep_all = TRUE)
indicationSelect= indicationSelect %>% sdf_explode(approvedIndications,keep_all = TRUE)
indicationSelect= indicationSelect %>%
  select(id,
         approvedIndications,
         indicationCount,
         disease,
         efoName,
         maxPhaseForIndication)
indication_df=indicationSelect%>% collect()
indication_df=indication_df[!duplicated(indication_df), ]
indication_inv=subset(indication_df,indication_df$maxPhaseForIndication==2 |indication_df$maxPhaseForIndication==3)
names(indication_approved)[[1]]="chemblIds"
indication_set=indication_approved
indication_set$approvedIndications=NULL

approved_MoA_indication=merge(approved_MoA_ready,indication_set,by="chemblIds",all.x=T)
approved_MoA_indication_clean=subset(two,!is.na(two$disease))
keep=subset(two,two$yearOfFirstApproval==2022)
approved_MoA_indication_ready=rbind(approved_MoA_indication,keep)
approved_MoA_indication_ready$var=paste(approved_MoA_indication_ready$disease,approved_MoA_indication_ready$targets,sep="_")
approved_MoA_indication_ready$trio=paste(approved_MoA_indication_ready$chemblIds,approved_MoA_indication_ready$disease,approved_MoA_indication_ready$targets,sep="_")

## read diseases dataset
diseasesPath="diseases"
diseases<- spark_read_parquet(sc,path = diseasesPath)
columns_diseases <- diseases %>%
  sdf_schema() %>%
  lapply(function(x) do.call(tibble, x)) %>%
  bind_rows()
#extract id,name,therapeuticAreas
diseaseSelect <- diseases %>%
  select(id,
         id,
         name,
         therapeuticAreas) %>% sdf_explode(therapeuticAreas,keep_all = TRUE)
disease_df=diseaseSelect%>% collect()
names(disease_df)[[1]]="disease"
approved_MoA_indication_disease=merge(approved_MoA_indication_clean,disease_df,by="disease")

## read associationByDataTypeDirect dataset
direct_datatype_df=parquet_readr("associationByDatatypeDirect")
direct_genetic_datatype_df=subset(direct_datatype_df,direct_datatype_df$datatypeId=="genetic_association")
direct_genetic_datatype_df=direct_genetic_datatype_df[!duplicated(direct_genetic_datatype_df), ]
dim(direct_genetic_datatype_df)
direct_genetic_datatype_df$score<-NULL
direct_genetic_datatype_df$evidenceCount<-NULL
head(direct_genetic_datatype_df)
names(direct_genetic_datatype_df)[[2]]="targets"
direct_genetic_datatype_df$var=paste(direct_genetic_datatype_df$diseaseId,direct_genetic_datatype_df$targets,sep="_")


final_trio=merge(approved_MoA_indication_ready,direct_genetic_datatype_df,by="var",all.x = T)


save_data=final_trio %>% rowwise() %>% 
  mutate(tradeNames = paste(tradeNames, collapse=',')) %>%
  ungroup()