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

#identify what E.coil gene matches best and the tophits
db <- read.fasta("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
str(db)
head(names(db))
myseqs <- db[which(names(db) %in% hits)] # extract the names of the top hits
myseqs <- c(myseqs,Ecoli) # add the Ecoli sequence
seqinr::write.fasta(myseqs,names=names(myseqs),file.out = "myseqs.fa")

#extract the names of the top hit percent identity, E-value and bit score
tophits <- db[which(names(db) %in% hits[1])]
tophit[1:3]
seqinr::write.fasta(tophit,names=names(tophit),file.out = "tophit.fa")
makeblastdb("tophit.fa", dbtype = "nucl","-parse_seqids")
res <- myblastn(myseq = myEcoli, db= "tophit.fa")
res
cat(res,fill=TRUE)


#Question 4

#create a mutated copy with 30 substitutions
myEcolimutator<- mutator(myEcoli,30)
res <- myblastn_tab(myseq = myEcolimutator, db= "tophit.fa")
res


#Question 5 

mutator #randomize with mutator
myblastn_tab
#create a function to summerise the result and report as a 0 or 1
myfunc <- function (myseq,nmut){
  mutseq <- mutator(myseq=myseq,nmut=nmut)
  res<-myblastn_tab(myseq=mutseq, db= "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
  if(is.null(res)){myres=0} else {myres=1}
  return(myres)}

#test the function
myfunc(myseq=myEcoli,nmut=50)
#repeat this routine 100 times to grt accurate result
replicate(n=100,expr = myfunc(myseq=myEcoli,nmut=50))
#summerise result to decimal number(between 0 and 1),to get a propostion of how many BLASTs were successful
mean (replicate(n=100,expr = myfunc(myseq=myEcoli,nmut=50)))
#all the values of nmut that would like to evaluate
n<-c(0,10,20,30,40,50,60,70,80,90,100)
#run replicate command to all those values
myfunc_rep <-function (nmut){mean (replicate(n=100,expr = myfunc(myseq=myEcoli,nmut=50)))}
#test the BLASTs result
finalres <-sapply(n,myfunc_rep)
finalres

#Follow same method for full sequence
myfunc(myseq=myEcoli,nmut=276)
replicate(n=100,expr = myfunc(myseq=myEcoli,nmut=276))
mean (replicate(n=100,expr = myfunc(myseq=myEcoli,nmut=276)))
n<-c(0,10,20,30,40,50,60,70,80,90,100)
myfunc_rep <-function (nmut){mean (replicate(n=100,expr = myfunc(myseq=myEcoli,nmut=276)))}
finalres <-sapply(n,myfunc_rep)
finalres

#Question 6

Bdata <- c(0,10,20,30,40,50,60,70,80,90,100)
Bdata <- as.data.frame(Bdata)
Bdata
Bdata$repsite<-c(1,1,1,0.95,0.75,0.49,0.25,0.12,0.05,0.02,0.01)
Bdata
plot(Bdata, type="b", col="red",col.main="blue", col.sub="green",
     main="How increasing no.random bases affects BLAST performance",sub="100 repeates,using sequence number 56",
     xlab="Number of site randomised", ylab="Propotion of successful BLASTs"
     )




