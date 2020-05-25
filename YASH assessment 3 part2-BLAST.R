#libraries that we need
library("seqinr")
library("R.utils")
library("rBLAST")

# Question 1
#Download the E.coli gene sequency from the Ensemble FTP page
download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz", destfile = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")
#uncompress the file
gunzip("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",overwrite=TRUE)
#create a blast database
makeblastdb("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa",dbtype = "nucl","-parse_seqids")

# Question 2
#Download sample file
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa", destfile = "sample.fa")
#Read file into R
Ecoli<-read.fasta("sample.fa")
str(Ecoli)
myEcoli<-Ecoli[[56]]
myEcoli
str(myEcoli)

#calculate length
seqinr::getLength(myEcoli)
#calculate GC content
seqinr::GC(myEcoli)

# Question 3
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R",destfile = "mutblast.R")
source("mutblast.R")

#test the function
res<-myblastn_tab(myseq = myEcoli,db= "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
str(res)
head(res)

#first three hits
hits <-as.character(res$sseqid[1:3])
hits

#identify what E.coil gene matches best
db<-read.fasta("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
str(db)
head(names(db))
# extract the names of the top hits
myseqs<-db[which(names(db) %in% hits[1])]
# add my selected sequence
myseqs<-c(myseqs,Ecoli)
seqinr::write.fasta(myseqs,names=names(myseqs),file.out = "myseqs.fa")

#extract the names of the top hit
#need to fix
tophit <- db[which(names(db) %in% hits[1])]
tophit[1:3]

seqinr::write.fasta(tophit,names=names(tophit),file.out = "tophit.fa")
makeblastdb("tophit.fa", dbtype = "nucl","-parse_seqids")
res <- myblastn(myseq = myEcoli, db= "tophit.fa")
cat(res,fill=TRUE)

#percent identity

# Question 4
myEcolimutator<- mutator(myEcoli,30)
res <- myblastn_tab(myseq = myEcolimutator, db= "tophit.fa")
res

#ecoli <-DNAStringSet(myEcoli)
str(ecoli)
ecoli <-toString(ecoli)
str(ecoli)
#ecoli <-s2c(ecoli)#
ecoli

#create a mutated copy 
ecoli_mut <- mutator(myseq = myEcoli, 50)
ecoli_mut
ecoli_mut_ <- DNAString(c2s(ecoli_mut))
ecoli_mut_
#need to fix
aln <- pairwiseAlignment (ecoli,ecoli_mut_)

str(ecoli)
str(ecoli_mut_)
pid(aln)
nmismatch(aln)

#Question 5
#write BLAST index (need to fix)
write.fasta (ecoli,names= "ecoli", file.out= "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
makeblastdb(file="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa", dbtype="nucl")

#test with 100 mismatches
ecoli_mut <- mutator (myseq=ecoli,100)
length(ecoli_mut)
res<-myblastn_tab(myseq = ecoli_mut,db= "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
res

ecoli_mut <- mutator (myseq=ecoli,200)
res<-myblastn_tab(myseq = ecoli_mut,db= "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
res

#Question 6








