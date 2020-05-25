#libraries that we need
library("seqinr")
library("R.utils")
library("rBLAST")

# Question 1

#Download the E.coli gene sequency from the Ensemble FTP page
download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz", destfile = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")
#use gunzip to uncompress the file
gunzip("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",overwrite=TRUE)
#create a blast database with makeblast function
makeblastdb("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa",dbtype = "nucl","-parse_seqids")

# Question 2

#Download sample file
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa", destfile = "sample.fa")
#Read file into R as Ecoli
Ecoli<-read.fasta("sample.fa")
str(Ecoli) #check data structure
#Read my allocated sequence (56) into R as myEcoli
myEcoli<-Ecoli[[56]]
#Lets check if the data has been imported properly
myEcoli
str(myEcoli)

#calculate length for my allocated gene
seqinr::getLength(myEcoli)
#calculate GC content for my allocated gene
seqinr::GC(myEcoli)


# Question 3

#sourced the function that ran blast searches in R
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R",destfile = "mutblast.R")
source("mutblast.R")
#test the function
res<-myblastn_tab(myseq = myEcoli,db= "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
#have a look at the blast results
head(res)
str(res)

#first three hits
hits <-as.character(res$sseqid[1:3])
hits

#identify what E.coil gene matches best
db<-read.fasta("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
str(db)
head(names(db))
head(db)

#Question 4

seqinr::write.fasta(myEcoli,names="myEcoli",file.out = "myEcoli.fa")
makeblastdb("myEcoli.fa", dbtype = "nucl","-parse_seqids")
res <- myblastn(myseq = myEcoli, db= "myEcoli.fa")
res
cat(res,fill=TRUE)

#create a mutated copy with 30 substitutions
myEcolimutator<- mutator(myEcoli,30)
res <- myblastn_tab(myseq = myEcolimutator, db= "myEcoli.fa")
res


#Question 5

#test with mismatches
myEcolimutator <- mutator(myEcoli,20)
res<-myblastn_tab(myseq = myEcolimutator,db= "myEcoli.fa")
res

myEcolimutator <- mutator(myEcoli,30)
res<-myblastn_tab(myseq = myEcolimutator,db= "myEcoli.fa")
res

myEcolimutator <- mutator(myEcoli,40)
res<-myblastn_tab(myseq = myEcolimutator,db= "myEcoli.fa")
res

myEcolimutator <- mutator(myEcoli,50)
res<-myblastn_tab(myseq = myEcolimutator,db= "myEcoli.fa")
res

myEcolimutator <- mutator(myEcoli,60)
res<-myblastn_tab(myseq = myEcolimutator,db= "myEcoli.fa")
res

myEcolimutator <- mutator(myEcoli,65)
res<-myblastn_tab(myseq = myEcolimutator,db= "myEcoli.fa")
res

