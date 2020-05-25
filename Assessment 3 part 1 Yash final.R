#Part 1 first section

#Question 1
# Downlod the data file of "gene_expression.tsv"
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/gene_expression.tsv",destfile="gene_expression.tsv")

# Read into R
x<-read.table("gene_expression.tsv")
# Lets check if the data has been imported properly
head(x)
str(x)

x <- read.table("gene_expression.tsv", header = TRUE , stringsAsFactors = FALSE , row.names = 1, sep="\t")
head(x)
str(x)

#Question 2
# Making new column which is mean of others 
x$Mean<- rowMeans(x)
head(x)

#Question 3
#list 10 genes with the highest mean expression
head(ord<-x[order(-x$Mean),],n=10)

#Question 4
# Number of genes with mean<10
which(B<-rowMeans(x)<10)
#(42124+1573)=43,697

#Question 5
#histogram plot of mean values
hist(x$Mean,breaks = 20,xlab="Genes",ylab="Mean",main="Mean values of the RNA sequency",col="green",border = "brown")
#save histogram in png format
png(file="Assessment part 1 first section histogram.png")
hist(x$Mean,breaks = 20,xlab="Genes",ylab="Mean",main="Mean values of the RNA sequency",col="green",border = "brown")
dev.off()

#Part 1 second section
#Question 6
# Downlod the data file of "growth_data.csv"
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/growth_data.csv",destfile="growth_data.csv")

# Read into R
y<-read.csv("growth_data.csv")

# Lets check data type and if the data has been imported properly
head(y)
str(y)

y <- read.csv("growth_data.csv", header=TRUE, sep="," ,stringsAsFactors = FALSE)
head(y)

#column names of data set
colnames(y)

#Question 7
#Seperate tree circumference in both sites
NE<-subset(y,Site=="northeast")
SW<-subset(y,Site=="southwest")
head(NE)
head(SW)
#Mean and standard deviation of tree circumference in northeast site at start(2004)
mean(NE$Circumf_2004_cm)
sd(NE$Circumf_2004_cm)
#Mean and standard deviation of tree circumference in southwest site at start(2004)
mean(SW$Circumf_2004_cm)
sd(SW$Circumf_2004_cm)
#Mean and standard deviation of tree circumference in northeast site at end(2019)
mean(NE$Circumf_2019_cm)
sd(NE$Circumf_2019_cm)
#Mean and standard deviation of tree circumference in southwest site at end(2019)
mean(SW$Circumf_2019_cm)
sd(SW$Circumf_2019_cm)

#Question 8
boxplot(NE$Circumf_2004_cm,NE$Circumf_2019_cm,SW$Circumf_2004_cm,SW$Circumf_2019_cm)
boxplot(NE$Circumf_2004_cm,NE$Circumf_2019_cm,SW$Circumf_2004_cm,SW$Circumf_2019_cm,ylab="Circumference(cm)",names=c("NE2004","SW2004","NE2019","SW2019"), main="Growth of Two Plantation Sites",las=2,col=c("yellow","blue"),border = "red",ntch=TRUE)

#save boxplot
png(file="Assessment part 1 second section boxplots.png")
boxplot(NE$Circumf_2004_cm,NE$Circumf_2019_cm,SW$Circumf_2004_cm,SW$Circumf_2019_cm,ylab="Circumference(cm)",names=c("NE2004","SW2004","NE2019","SW2019"), main="Growth of Two Plantation Sites",las=2,col=c("yellow","blue"),border = "red",ntch=TRUE)
dev.off()

#Question 9
# Mean growth of past 10 yeras
NE$growth<-NE$Circumf_2019_cm-NE$Circumf_2009_cm
head(NE)

SW$growth<-SW$Circumf_2019_cm-SW$Circumf_2009_cm
head(SW)

#Question 10
#Use t.test to estimate the p-value for 10 years growth at 2 sites
t.test(SW$growth,NE$growth)
#Use wilcox.test to estimate the p-value for 10 years growth at 2 sites
wilcox.test(SW$growth,NE$growth)