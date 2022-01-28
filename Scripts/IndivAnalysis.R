
options(repos = c(CRAN = "http://cran.rstudio.com"))

install.packages("Biostrings")
install.packages("RMySQL")
install.packages("DBI")
install.packages("car") 
install.packages("tidyr")
install.packages("plyr")
install.packages("reshape")
install.packages("scales")
install.packages("ggdendro")
install.packages("plotly")
install.packages("ggdendro")
install.packages("heatmaply")
install.packages("devtools")
install.packages(c('devtools','curl'))
devtools::install_github('hadley/ggplot2',force=TRUE)
install.packages("phytools")
install.packages("phylogram")
install.packages("Rmisc")
install.packages("phylotools")
install.packages("wesanderson")
install.packages("forcats")
devtools::install_github('slowkow/ggrepel')
install.packages("gplots")

library("Biostrings")
library(ggplot2)
library(RMySQL)
library(DBI)
library(car)
library(tidyr)
library(plyr)
library(reshape)
library(scales)
library(ggdendro)
library(plotly)
library(heatmaply)
library(grid)
library(gridExtra)
library(ape)
library(phytools)
library("phylogram")
library(Rmisc)
library(phylotools)
library(wesanderson)
library(forcats)
library(gridExtra)
library(gtable)
library(grid)
library(devtools)
library(ggrepel)
library(gplots)

source("https://bioconductor.org/biocLite.R")
biocLite("ggtree")


#build tools stuff
Sys.getenv('PATH')
system('g++ -v')
system('where make')

#Connect to Database
TEdb = dbConnect(MySQL(), user='drew', password='drew', host='localhost', dbname="TEdb")
dbSendQuery(TEdb, "drop table if exists tes, species")

dbListTables(TEdb)
dbListFields(TEdb,"TEs")
dbListFields(TEdb,"species")

#Display TEdb data
data = fetch(dbSendQuery(TEdb, "SELECT * FROM TEs"), n=-1)
data = fetch(dbSendQuery(TEdb, "SELECT * FROM species"),n=-1)

# Code to roll back versions. Which is why session_info says URL now for source
devtools::install_version("DBI", version = "0.5-1", repos = "http://cran.us.r-project.org", dependencies = TRUE)
devtools::install_version("RMySQL", version = "0.10.9", repos = "http://cran.us.r-project.org", dependencies = TRUE)

#disconnect all sessions from MySQL
dbDisconnectAll <- function(){
  ile <- length(dbListConnections(MySQL())  )
  lapply( dbListConnections(MySQL()), function(x) dbDisconnect(x) )
  cat(sprintf("%s connection(s) closed.\n", ile))
}
dbDisconnectAll()


#Set the palette for graphs
palette <- scale_fill_hue(l=50, c=100)
nogrid <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank(), axis.line = element_line(colour = "black"))
vertText <- theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.35), axis.text = element_text(face = "bold",size = 10), axis.title = element_text(face="bold"))
nogap <- scale_y_continuous(expand = c(0,0)) 
border <- geom_bar(stat = "identity",color="black") 
darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

#Graph formatting

factorViri <- function(x) {
  x <- trimws(x, which = c("both", "left", "right"))
  
  
  factor(x, levels = c("A.halleri", "A.thaliana", "A.coerulea", "B.distachyon", "B.stricta", "C.annuum", "C.clementina", "C.grandiflora", "C.rubella", "C.quinoa", "E.grandis", "E.salsugineum", "F.vesca", "G.arboreum", "G.hirsutum", "G.raimondii", "G.max", "K.laxiflora", "L.usitatissimum", "M.polymorpha", "O.sativa", "P.persica", "P.virgatum", "R.communis", "S.italica", "S.polyrhiza", "S.lycopersicum", "S.tuberosum", "T.pratense", "U.gibba", "V.vinifera", "Z.marina", "Z.mays"))
}

## RepeatMasker functions ### Landscape plots ############################
#setwd("D:\\Documents\\School\\Current Courses\\BIOL 440\\Data\\RMout")
#setwd("D:\\Documents\\School\\Current Courses\\BIOL 440\\Data\\Arabidopsis\\New\\Final Output #Files")
#setwd("D:\\Documents\\School\\Current Courses\\BIOL 440\\Data\\Novel Finding Data\\out")
setwd("D:\\Documents\\School\\Current Courses\\BIOL 440\\Data\\Arabidopsis\\Repbase")
#setwd("D:\\Documents\\School\\Current Courses\\BIOL 440\\Data\\library comparison")
#setwd("D:\\Documents\\School\\Current Courses\\BIOL 440\\Data\\RMask - Repbase\\Viri Output")
TEdb = dbConnect(MySQL(), user='drew', password='drew', host='localhost', dbname="TEdb")
dbSendQuery(TEdb, "drop table if exists TEs, species")
rm(TEs.df)
rm(species.df)

input.filelist <- list.files()

for (input.file in input.filelist){
  
  filename <- toString(input.file)
  filename <- strsplit(filename,".",fixed = TRUE)[[1]]
  
  species <- filename[1]
  genomesize <- as.numeric(as.character(filename[2]))
  TEsize <- as.numeric(as.character(filename[3]))
  
  # if the merged dataset doesn't exist, create it
  if (!exists("TEs.df")){
    TEs.df <- read.table(input.file, header=TRUE, sep="\t", comment.char = "", fill = FALSE)
    TEs.df <- cbind(species,TEs.df)
  } else {
  # if the merged dataset does exist, append to it
    TEs.temp <-read.table(input.file, header=TRUE, sep="\t", comment.char = "", fill = FALSE)
    TEs.temp <- cbind(species,TEs.temp)
    TEs.df<-rbind(TEs.df, TEs.temp)
    rm(TEs.temp)
  }

  if (!exists("species.df")){
    species.df <- data.frame(species, genomesize, TEsize)
  } else {
    species.temp <- data.frame(species,genomesize, TEsize)
    species.df<-rbind(species.df, species.temp)
    rm(species.temp)
  }
}

TEs.df <- subset(TEs.df, Rclass != "Simple_repeat")
TEs.df <- subset(TEs.df, Rclass != "Satellite")
TEs.df <- subset(TEs.df, Rclass != "Low_complexity")
TEs.df <- subset(TEs.df, Rclass != "ARTEFACT")
TEs.df <- subset(TEs.df, Rclass != "RathE1_cons")
TEs.df <- subset(TEs.df, Rclass != "RathE2_cons")
TEs.df <- subset(TEs.df, Rclass != "RathE3_cons")
TEs.df <- subset(TEs.df, Rclass != "Unassigned")
TEs.df <- subset(TEs.df, Rclass != "MobileElement")
TEs.df <- subset(TEs.df, Rclass != "nonLTR")
TEs.df <- subset(TEs.df, Rclass != "Retroelement")
TEs.df <- subset(TEs.df, Rclass != "Unspecified")


levels(TEs.df$Rclass)
TEs.df$species <- factor(sub("^\(.*)\\?", "\\1", TEs.df$Rclass)) #Remove the ? endings
TEs.df$Rclass <- factor(sub("^(.*)\\?", "\\1", TEs.df$Rclass)) #Remove the ? endings
TEs.df$Rfam <- factor(sub("^(.*)\\?", "\\1", TEs.df$Rfam)) #Remove the ? endings
TEs.df$Rfam <- factor(sub("^(.*)-.*", "\\1", TEs.df$Rfam)) #Remove the subfamily
#TEs.df$Rfamily <- sub("^(.*)?", "\\1", TEs.df$Rfamily)
levels(TEs.df$Rclass)

#TEs.df$species <- factor(sub("^(.)", "\\1.", TEs.df$species)) #Remove the ? endings
#species.df$species <- factor(sub("^(.)", "\\1.", species.df$species)) #Remove the ? endings

TEs.df$Rclass<-recode(TEs.df$Rclass,"'RC'='DNA'")
TEs.df$Rclass<-recode(TEs.df$Rclass,"'Retroposon'='LTR'")

TEs.df <- subset(TEs.df, Rclass != "snRNA")
TEs.df <- subset(TEs.df, Rclass != "tRNA")
TEs.df <- subset(TEs.df, Rclass != "Unknown")
TEs.df <- subset(TEs.df, Rclass != "rRNA")
TEs.df <- subset(TEs.df, Rclass != "Satellite")
#TEs.df <- subset(TEs.df, Rclass != "SINE")
TEs.df <- subset(TEs.df, Rclass != "Other")


TEs.df$Rclass<-recode(TEs.df$Rclass,"'LTR'='LTR'")
TEs.df$Rclass<-recode(TEs.df$Rclass,"'LINE'='LINE'")
TEs.df$Rclass<-recode(TEs.df$Rclass,"'SINE'='SINE'")
TEs.df$Rclass<-recode(TEs.df$Rclass,"'DNA'='DNA'")

dbWriteTable(TEdb,name='TEs',value=TEs.df)
dbWriteTable(TEdb,name='species',value=species.df)


##get single TE
#te.to.get <- "Kolobok"
te.to.get <- "Harbinger"
#te.to.get <- "Merlin"
#te.to.get <- "TRIM"
got.te <- fetch(dbSendQuery(TEdb, sprintf("SELECT * FROM TEs WHERE Rfam=\"%s\"",te.to.get)), n=-1)




##Get distinct
TEs.subfamily.distinct <- fetch(dbSendQuery(TEdb, "SELECT SPECIES, COUNT(DISTINCT Rclassfam) FROM TEs GROUP BY species"), n =-1)
  
  
theme(legend.position="top", legend.title = element_blank(),axis.title.x=element_blank(), axis.title.y=element_blank())

#To reverse order of coord clipped charforcats::fct_rev(reorder(Species,Species))

##### Genome size vs TE size
genomes <- fetch(dbSendQuery(TEdb, "SELECT species, genomesize FROM species"), n=-1)

tes <- fetch(dbSendQuery(TEdb, "SELECT species, SUM(LEN_MASKED) FROM TEs GROUP BY species"), n=-1)
gen.tes <- merge(genomes, tes, by="species", all=TRUE)
names(gen.tes) <- c("Species", "Genomic","TEs")
gen.tes <- melt(gen.tes, id.vars='Species')
names(gen.tes) <- c("Species", "Content","Value")

#gen.tes<- factorViri(gen.tes$Species)

ggplot(gen.tes, aes(Species, Value/1000000)) +
  theme(legend.title = element_blank()) +
  geom_bar(aes(fill = Content), width = .8, position = position_dodge(width=0.8), stat="identity") +  
   vertText + nogrid +  palette + labs(x = "Species", y = "DNA Length (Mbp)") + 
  scale_y_continuous(breaks = pretty((gen.tes$Value)/1000000, n = 20),expand = c(0,0)) 

#####Proportion of TEs
genomes <- fetch(dbSendQuery(TEdb, "SELECT species, genomesize FROM species"), n=-1)

tes <- fetch(dbSendQuery(TEdb, "SELECT species, SUM(LEN_MASKED) FROM TEs GROUP BY species"), n=-1)
gen.tes <- merge(genomes, tes, by="species", all=TRUE)
names(gen.tes) <- c("Species", "Genome","TEs")
gen.tes <- transform(gen.tes, prop = 100*TEs/Genome)






genomes <- fetch(dbSendQuery(TEdb, "SELECT species, genomesize FROM species"), n=-1)
tes <- fetch(dbSendQuery(TEdb, "SELECT species, SUM(LEN_MASKED) FROM TEs GROUP BY species"), n=-1)
names(gen.tes) <- c("Species", "Genome","TEs")

gen.tes$prop <- 100*(gen.tes$TEs/gen.tes$Genome)



gen.tes <- melt(gen.tes, id.vars='Species')
ggplot(gen.tes, aes(Species, (TEs/Genome)*100)) + geom_bar(aes(fill = variable), 
                                                             width = .8, position = position_dodge(width=0.8), stat="identity") +  
  theme(legend.position="top", legend.title = 
          element_blank(),axis.title.x=element_blank(), 
        axis.title.y=element_blank()) +nogap +vertText + nogrid +  palette +labs(x = "Species", y = "Genome Size (Mbp>")

####By length
#By class
TE.len.class <- fetch(dbSendQuery(TEdb, "SELECT species, Rclass, SUM(LEN_MASKED) FROM TEs GROUP BY species, Rclass"), n=-1)
names(TE.len.class) <- c("Species", "Class", "Total.Length")
ggplot(data = TE.len.class, aes(x = Species, y = Total.Length/1000000, fill = Class)) +  geom_bar(stat = "identity",color="black") + labs(x = "Library", y = "Number of Bases (Mbp)") + nogrid  + guides(fill=guide_legend(title="TE Class: ")) + palette +vertText +nogap 

TE.len.class <- edit(TE.len.class)


ggplot(data = TE.len.class, aes(x = Species, y = log(Total.Length), fill = Class)) +  geom_bar(stat = "identity") + coord_flip()

####By copy number
#By class
TE.count.class <- fetch(dbSendQuery(TEdb, "SELECT species, Rclass, SUM(FRG_NB) FROM TEs GROUP BY species, Rclass"), n=-1)
names(TE.count.class) <- c("Species", "Class", "Total.Count")

ggplot(data = TE.count.class, aes(x = Species, y = Total.Count/1000, fill = Class)) +  geom_bar(stat = "identity",color="black") + labs(x = "Library", y = "Total Count of TEs (1000's)") + nogrid + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.35)) + guides(fill=guide_legend(title="TE Class: ")) + scale_y_continuous(expand = c(0,0)) +palette



ggplot(data = TE.count.class, aes(x = Library, y = (Total.Count/1000), fill = Class)) +  geom_bar(stat = "identity")  + nogrid + palette +vertText +border +  labs( x = "Library", y = "Total Count of TEs (1000's)")







#By family
TE.len.family <- fetch(dbSendQuery(TEdb, "SELECT species, CONCAT(`Rclass`,'.', `Rfam`), SUM(LEN_MASKED) FROM TEs GROUP BY species, CONCAT(`Rclass`,'|', `Rfam`)"), n=-1)
names(TE.len.family) <- c("Species", "Class.Family", "Total.Length")

ggplot(data = TE.len.family, aes(x = Species, y = log(Total.Length), fill = Class.Family)) +  geom_bar(stat = "identity",color="black") + labs(x = "Library", y = "Number of Bases (Mbp)") + nogrid + theme(legend.position="top") + guides(fill=guide_legend(title="TE Class: ")) + scale_y_continuous(expand = c(0,0)) +palette + labs(x = "Library", y = "Log( Number of Bases [Mbp] )") 




ggplot(data = TE.len.family, aes(x = Species, y = log(Total.Length), fill = Class.Family)) +  geom_bar(stat = "identity") + coord_flip() + nogrid 

#By total TEs
TE.len.total <- fetch(dbSendQuery(TEdb, "SELECT species, SUM(LEN_MASKED) FROM TEs GROUP BY species"), n=-1)
names(TE.len.total) <- c("Species", "Total.Length")
ggplot(data = TE.len.total, aes(x = Species, y = Total.Length)) +  geom_bar(stat = "identity") + coord_flip() + palette + nogrid + palette
ggplot(data = TE.len.total, aes(x = Species, y = log(Total.Length))) +  geom_bar(stat = "identity") + coord_flip() + palette + nogrid + palette
ggplot(data = TE.len.total, aes(x = Species, y = Total.Length^(1/3))) +  geom_bar(stat = "identity") + coord_flip() + palette + nogrid + palette


#As proportion of TE total length
#By class
TE.len.class <- fetch(dbSendQuery(TEdb, "SELECT species, Rclass, SUM(LEN_MASKED) FROM TEs GROUP BY species, Rclass"), n=-1)
TE.len.total <- fetch(dbSendQuery(TEdb, "SELECT species, SUM(LEN_MASKED) FROM TEs GROUP BY species"), n=-1)
TE.prop.class <- merge(TE.len.class, TE.len.total, by="species", all=TRUE)
names(TE.prop.class) <- c("Species", "Class", "Total.class", "Total.TEs")
#TE.prop.class$Species <- factorViri(TE.prop.class$Species)
p1 <- ggplot(data = TE.prop.class, aes(x = Species, y = (100*Total.class/Total.TEs), width=1, fill = Class)) +  labs(x = "", y = "Relative proportion of TEs") + nogrid + palette + vertText + nogap + border  + theme(plot.margin = unit(c(-1.9,1,1,1), "cm")) + scale_x_discrete(position = "top")  + theme(  axis.text.x = element_blank()) 
p1
#By family


#By all
TE.len.class <- fetch(dbSendQuery(TEdb, "SELECT species, SUM(LEN_MASKED) FROM TEs GROUP BY species"), n=-1)
TE.len.total <- fetch(dbSendQuery(TEdb, "SELECT species, SUM(LEN_MASKED) FROM TEs GROUP BY species"), n=-1)
TE.prop.class <- merge(TE.len.class, TE.len.total, by="species", all=TRUE)
names(TE.prop.class) <- c("Species", "Total.TEs")
TE.prop.class$Species <- factorViri(TE.prop.class$Species)
p1 <- ggplot(data = TE.prop.class, aes(x = Species, y = (100*Total.class/Total.TEs), width=1)) +  labs(x = "", y = "Relative proportion of TEs") + nogrid + palette + vertText + nogap + border  + theme(plot.margin = unit(c(-1.9,1,1,1), "cm")) + scale_x_discrete(position = "top")  + theme(  axis.text.x = element_blank()) 
p1

####By proportion of genome length
#By class
TE.genlen.class <- fetch(dbSendQuery(TEdb, "SELECT species, Rclass, SUM(LEN_MASKED) FROM TEs GROUP BY species, Rclass"), n=-1)
genome.len.total <- fetch(dbSendQuery(TEdb, "SELECT species, genomesize FROM species"), n=-1)
TE.genprop.class <- merge(TE.genlen.class, genome.len.total, by="species", all=TRUE)
names(TE.genprop.class) <- c("Species", "Class", "Total.class", "Total.genome")
#TE.genprop.class$Species <- factorViri(TE.genprop.class$Species)
p2 <- ggplot(data = TE.genprop.class, aes(x = Species, y = (100*(Total.class/Total.genome)), width=1, fill = Class)) +  labs( x = "", y = "TEs as Percentage of Genome") + nogrid + palette + vertText + nogap + border + theme(plot.margin = unit(c(1,1,1,1), "cm"))
p2 

g1 <- ggplotGrob(p2)
#g1 <- gtable_add_cols(g1, unit(0,"mm")) # add a column for missing legend
g2 <- ggplotGrob(p1)

ncol(g2)
ncol(g1)
g <- rbind(g1, g2, size="first") # stack the two plots
g$widths <- unit.pmax(g1$widths, g2$widths) # use the largest widths
# center the legend vertically
#g$layout[grepl("guide", g$layout$name),c("t","b")] <- c(1,nrow(g))
grid.newpage()
grid.draw(g)

#gridExtra::grid.arrange(p2, p1, ncol = 1, height = 1:1)


####By proportion of TE copy number
#By class
TE.count.class <- fetch(dbSendQuery(TEdb, "SELECT species, Rclass, SUM(FRG_NB) FROM TEs GROUP BY species, Rclass"), n=-1)
TE.count.total <- fetch(dbSendQuery(TEdb, "SELECT species, SUM(FRG_NB) FROM TEs GROUP BY species"), n=-1)
TE.count.prop.class <- merge(TE.count.class, TE.count.total, by="species", all=TRUE)
names(TE.count.prop.class) <- c("Species", "Class", "Total.class", "Total.TEs")
ggplot(data = TE.count.prop.class, aes(x = Species, y = (100*Total.class/Total.TEs), fill = Class)) +  geom_bar(stat = "identity") + coord_flip() + nogrid  + palette


###Genome length vs TE length
#All
#TE.length.comp.class <- fetch(dbSendQuery(TEdb, "SELECT TEs.species, species.genomesize, SUM(TEs.LEN_MASKED) FROM TEs, species WHERE species.species = TEs.species GROUP BY species"), n=-1)
TE.length.comp.class <- fetch(dbSendQuery(TEdb, "SELECT species, genomesize, TEsize from species"), n=-1)
names(TE.length.comp.class) <- c("Species", "Genome.Size", "TEs.Length")
#Non-log
ggplot(TE.length.comp.class, aes(x=(100*TEs.Length/Genome.Size), y=(Genome.Size))) + geom_point() + geom_smooth(method=lm) + geom_text(aes(label=substring(Species,1,3),hjust=-0.3)) + nogrid + palette
#Log
ggplot(TE.length.comp.class, aes(x=(log(100*TEs.Length/Genome.Size)), y=log(Genome.Size/1000000))) + geom_point(size = 1) + geom_smooth(method=lm) + geom_text_repel(aes(label=gsub('\\.', '', substring(Species,1,2))),size=3,point.padding = 0.1) + nogrid + palette + labs( x = "Log( Total TE Percentage )", y = "Log( Genome Size in Mbp )")

ggplot(TE.length.comp.class, aes(x=(log(100*TEs.Length/Genome.Size)), y=log(Genome.Size/1000000))) + geom_point(size = 1) + geom_smooth(method=lm) + geom_text_repel(aes(label=Species)),size=3,point.padding = 0.1) + nogrid + palette + labs( x = "Log( Total TE Percentage )", y = "Log( Genome Size in Mbp )")
# + geom_point(size=2, shape=23)


###Genome length vs TE length
class.type <- "DNA"
#class.type <- "LTR"
#class.type <- "LINE"
#class.type <- "SINE"
TE.length.comp.class <- fetch(dbSendQuery(TEdb, sprintf("SELECT TEs.species, species.genomesize, SUM(TEs.LEN_MASKED) FROM TEs, species WHERE species.species = TEs.species AND Rclass='%s' GROUP BY species",class.type)), n=-1)
#TE.length.comp.class <- fetch(dbSendQuery(TEdb, "SELECT TEs.species, species.genomesize, species.TEsize, TEs.Rclass from TEs, species WHERE species.species = TEs.species AND TEs.Rclass='%s'",class.type)), n=-1)
names(TE.length.comp.class) <- c("Species", "Genome.Size", "TEs.Length")
#Non-log
ggplot(TE.length.comp.class, aes(y=(100*TEs.Length/Genome.Size), x=(Genome.Size))) + geom_point() + geom_smooth(method=lm) + geom_text_repel(aes(label=substring(Species,1,2))) + nogrid + palette
#Log
ggplot(TE.length.comp.class, aes(x=(log(100*TEs.Length/Genome.Size)), y=log(Genome.Size))) + geom_point() + geom_smooth(method=lm) + geom_text_repel(aes(label=substring(Species,1,2))) + nogrid + palette
# + geom_point(size=2, shape=23)

###Genome length vs TE lengths for all classes

TE.length.comp.class <- fetch(dbSendQuery(TEdb, "SELECT TEs.species, species.genomesize, SUM(TEs.LEN_MASKED), species.TEsize, TEs.Rclass FROM TEs, species WHERE species.species = TEs.species GROUP BY species,Rclass"), n=-1)
#TE.length.comp.class <- fetch(dbSendQuery(TEdb, "SELECT TEs.species, species.genomesize, species.TEsize, TEs.Rclass from TEs, species WHERE species.species = TEs.species AND TEs.Rclass='%s'",class.type)), n=-1)
names(TE.length.comp.class) <- c("Species", "Genome.Size", "TEs.Length", "TEs.totallength", "Class")
TE.length.comp.class 


#Non-log
ggplot(TE.length.comp.class, aes(y=(100*TEs.Length/TEs.totallength), x=(Genome.Size),colour=Class)) + geom_smooth(method=lm) + nogrid + palette  + geom_jitter(aes(colour = Class),size=1.5,width = 0.3)
#Log
ggplot(TE.length.comp.class, aes(y=(log(100*TEs.Length/TEs.totallength)), x=log(Genome.Size/1000000),colour=Class)) + geom_smooth(method=lm)  + nogrid + palette +  labs( x = "Log ( Genome Size (Mbp) )", y = "Log ( Percentage of of Genome )") + geom_jitter(aes(colour = Class),size=1.5,width = 0.3) + theme(legend.position="top")

##By family for class
#class.type <- "DNA"
class.type <- "LTR"
#class.type <- "LINE"
#class.type <- "SINE"
TE.len.family <- fetch(dbSendQuery(TEdb, sprintf("SELECT species, CONCAT(`Rclass`,'.', `Rfam`), SUM(LEN_MASKED) FROM TEs WHERE Rclass='%s' GROUP BY species, CONCAT(`Rclass`,'.', `Rfam`)", class.type)), n=-1)
TE.len.class <- fetch(dbSendQuery(TEdb, sprintf("SELECT species, SUM(LEN_MASKED) FROM TEs WHERE Rclass='%s' GROUP BY species, Rclass", class.type)), n=-1)
TE.prop.class <- merge(TE.len.family, TE.len.class, by="species", all=TRUE)
names(TE.prop.class) <- c("Species", "Subfamily", "Total.class", "Total.TEs")
ggplot(data = TE.prop.class, aes(x = Species, y = (100*Total.class/Total.TEs), fill = Subfamily)) +  geom_bar(stat = "identity") + coord_flip() + grid + palette

###By subfamily for class
#class.type <- "DNA"
class.type <- "LTR"
#class.type <- "LINE"
#class.type <- "SINE"
TE.len.subfamily <- fetch(dbSendQuery(TEdb, sprintf("SELECT species, CONCAT(`Rclass`,'.', `Rfam`), SUM(LEN_MASKED) FROM TEs WHERE Rclass='%s' GROUP BY species, CONCAT(`Rclass`,'.', `Rfam`)", class.type)), n=-1)
TE.len.class <- fetch(dbSendQuery(TEdb, sprintf("SELECT species, SUM(LEN_MASKED) FROM TEs WHERE Rclass='%s' GROUP BY species", class.type)), n=-1)
TE.prop.class <- merge(TE.len.subfamily, TE.len.class, by="species", all=TRUE)
names(TE.prop.class) <- c("Species", "Subfamily", "Total.class", "Total.TEs")
ggplot(data = TE.prop.class, aes(x = Species, y = (100*Total.class/Total.TEs), fill = Subfamily)) +  geom_bar(stat = "identity") + coord_flip() + grid + palette

##TE content grid
#TE.len.family <- fetch(dbSendQuery(TEdb, sprintf("SELECT species, CONCAT(`Rclass`,'.', `Rfam`), SUM(LEN_MASKED) FROM TEs WHERE Rclass='%s' GROUP BY species, CONCAT(`Rclass`,'.', `Rfam`)", class.type)), n=-1)
#TE.len.class <- fetch(dbSendQuery(TEdb, sprintf("SELECT species, SUM(LEN_MASKED) FROM TEs WHERE Rclass='%s' GROUP BY species, Rclass", class.type)), n=-1)
TE.len.family <- fetch(dbSendQuery(TEdb, "SELECT species, CONCAT(`Rclass`,'.', `Rfam`), SUM(LEN_MASKED) FROM TEs GROUP BY species, CONCAT(`Rclass`,'.', `Rfam`)"), n=-1)
TE.len.class <- fetch(dbSendQuery(TEdb, "SELECT species, SUM(LEN_MASKED) FROM TEs GROUP BY species"), n=-1)
TE.prop.class <- merge(TE.len.family, TE.len.class, by="species", all=TRUE)
names(TE.prop.class) <- c("Species", "Subfamily", "Total.class", "Total.TEs")
TE.prop.class$Prop.TE <- (TE.prop.class$Total.class/TE.prop.class$Total.TEs)*100

dendro.TE.prop.class <- TE.prop.class
dendro.TE.prop.class$Total.class <- NULL
dendro.TE.prop.class$Total.TEs <- NULL
t.dendro.TE.prop.class <- spread(dendro.TE.prop.class, Subfamily, Prop.TE)
t.dendro.TE.prop.class[is.na(t.dendro.TE.prop.class)] <- 0

t.dendro.TE.prop.class  <- data.frame(t.dendro.TE.prop.class [,-1], row.names=t.dendro.TE.prop.class [,1])

x <- as.matrix(t.dendro.TE.prop.class)

dd.col <- as.dendrogram(hclust(dist(x),method="single"))
dd.row <- as.dendrogram(hclust(dist(t(x)),method="single"))
dx <- dendro_data(dd.row)
dy <- dendro_data(dd.col)

TE.dendro <- ggplot(dy$label) +
  geom_segment(data = dy$segments, aes(x=x, y=y, xend=xend, yend=yend)) +
  labs(x = "", y = "") + theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank()) + coord_flip() + scale_y_reverse(expand=c(0.2, 0)) +
  #Remove below if labels are unwanted
  geom_text(aes(x = x, y = y, label = label, hjust = 0), data= label(dy))

TE.dendro

heatmap.TE.prop.class <- TE.prop.class
heatmap.TE.prop.class$Prop.TE.bin <- findInterval(heatmap.TE.prop.class$Prop.TE, c(0,0.001,0.1,1,5,50))

heatmap.TE.prop.class$Prop.TE <- NULL
heatmap.TE.prop.class$Total.class <- NULL
heatmap.TE.prop.class$Total.TEs <- NULL
t.heatmap.TE.prop.class <- spread(heatmap.TE.prop.class, Subfamily, Prop.TE.bin)
t.heatmap.TE.prop.class[is.na(t.heatmap.TE.prop.class)] <- 0
t.heatmap.TE.prop.class  <- data.frame(t.heatmap.TE.prop.class [,-1], row.names=t.heatmap.TE.prop.class [,1])
heatmapdata <- t.heatmap.TE.prop.class 


x <- as.matrix(t.heatmap.TE.prop.class)

col.ord <- order.dendrogram(dd.col)
row.ord <- order.dendrogram(dd.row)
xx <- x[col.ord, row.ord]
xx_names <- attr(xx, "dimnames")
df <- as.data.frame(xx)
colnames(df) <- xx_names[[2]]
df$Species <- xx_names[[1]]
df$Species <- with(df, factor(Species, levels=Species, ordered=TRUE))
mdf <- reshape2::melt(df, id.vars="Species")
#mdf <- mdf[order(as.character(mdf$variable)),]

mdf$variable <- factor(mdf$variable, levels = mdf$variable[order(as.character(mdf$variable))])

TE.heatmap <- ggplot(mdf, aes(x = reorder(variable, as.character(variable),ordered=TRUE), y = Species)) + geom_tile(aes(fill = value), colour = "white", size = 2.5) + scale_fill_gradient(low = "white" , high = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.0, vjust=0), axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +scale_x_discrete(position = "top") +
  coord_equal()
TE.heatmap

TE.dendro

grid.newpage()
print(TE.heatmap, vp=viewport(0.8, 0.8, x=0.55, y=0.4))
print(TE.dendro, vp=viewport(0.2, 0.8, x=0.1, y=0.4))

TE.plot <- subplot(TE.dendro, TE.heatmap)
TE.plot <- layout(TE.plot,
                margin = list(l = 150,
                              r = 0,
                              b = 50,
                              t = 0
                )
)

TE.plot


# Extract dendrograms for rows and columns from 'heat'
row.dendro <- dd.row
col.dendro <- dd.col

# Convert dendrograms to nwk (via .hcclust and .phylo formats!)
as.hclust (row.dendro)  ->row.hcclust
as.phylo  (row.hcclust) ->row.phylo
write.tree(row.phylo)   ->row.nwk

as.hclust (col.dendro)  ->col.hcclust
as.phylo  (col.hcclust) ->col.phylo
write.tree(col.phylo)   ->col.nwk

row.nwk
col.nwk

ggplot() + labs(x="0%     0 < % < 0.01   0.01 < % < 0.1       0.1 < % < 1     1 < % < 5     5 < %    ")
ggplot() + labs(x="Percentage of TEs")

## RepeatMasker functions ### Generate tree from Newick ############################
setwd("D:\\Documents\\School\\Current Courses\\BIOL 440\\Data")
species.tree <- read.newick("Newick.tre")
species.tree<-collapse.singles(species.tree)
plotTree(species.tree)

plot(species.tree)

species.tree <- read.nexus("Nexus.txt")
species.tree


read.phylip("Newick.tre")

x <- as.dendrogram("Newick.tre")

x <- read.dendrogram("Newick.tre")
plot(x)
ggdendrogram(as.dendrogram(species.tree))
ggdendrogram(x, rotate = TRUE, theme_dendro = TRUE)

########################3Box plot for subfamilies###########3

data <- fetch(dbSendQuery(TEdb, sprintf("SELECT * FROM landscapes WHERE (Rclass='LTR' OR Rclass='SINE' OR Rclass='LINE' OR Rclass='DNA')", species)), n=-1)
data$row_names <- NULL
data$Rclass <- NULL
long_data <- gather(data, bin, value, 1:40)
box_data <- long_data %>%
  group_by(species, bin) %>%
  summarise(val=sum(value))

TE.len.family <- fetch(dbSendQuery(TEdb, "SELECT species, CONCAT(`Rclass`,'.', `Rfam`), SUM(LEN_MASKED) FROM TEs GROUP BY species, CONCAT(`Rclass`,'|', `Rfam`)"), n=-1)
names(TE.len.family) <- c("Species", "Family", "Length")
ggplot(TE.len.family, aes(Family, log(Length)))+ geom_boxplot() + nogrid +vertText  + labs( x = "Kimura Distance", y = "Log(Proportion of Genome)") + geom_boxplot(fill = "white")

#By family

ggplot(data = TE.len.family, aes(x = Species, y = Total.Length, fill = Class.Family)) +  geom_bar(stat = "identity") + coord_flip() + nogrid + palette

## RepeatMasker functions ### Landscape plots ############################
#setwd("D:\\Documents\\School\\Current Courses\\BIOL 440\\Data\\Viridiplantae")
setwd("D:\\Documents\\School\\Current Courses\\BIOL 440\\Data\\Arabidopsis\\Landscapes")
TEdb = dbConnect(MySQL(), user='drew', password='drew', host='localhost', dbname="TEdb")
dbSendQuery(TEdb, "drop table if exists landscapes, landscapeSpecies")
rm(landscapes.df)
rm(landscapeSpecies.df)

input.filelist <- list.files()

for (input.file in input.filelist){

  filename <- toString(input.file)
  filename <- strsplit(filename,".",fixed = TRUE)[[1]]
  
  species <- filename[1]
  genomesize <- as.numeric(as.character(filename[2]))
  
  # if the merged dataset doesn't exist, create it
  if (!exists("landscapes.df")){
    landscapes.df <- read.table(input.file, header=TRUE, sep="\t", comment.char = "", fill = TRUE)
    landscapes.df.Rclass <- landscapes.df$Rclass
    landscapes.df$Rclass <- NULL
    landscapes.df <- (landscapes.df*100)/genomesize
    landscapes.df$Rclass <- landscapes.df.Rclass
    landscapes.df$species <- species
  } else {
    # if the merged dataset does exist, append to it
    landscapes.temp <-read.table(input.file, header=TRUE, sep="\t", comment.char = "", fill = TRUE)
    landscapes.temp.Rclass <- landscapes.temp$Rclass
    landscapes.temp$Rclass <- NULL
    landscapes.temp <- (landscapes.temp*100)/genomesize
    landscapes.temp$Rclass <- landscapes.temp.Rclass
    landscapes.temp$species <- species
    landscapes.df<-rbind(landscapes.df, landscapes.temp)
    rm(landscapes.temp)
  }
  
  if (!exists("landscapeSpecies.df")){
    landscapeSpecies.df <- data.frame(species, genomesize)
  } else {
    landscapeSpecies.temp <- data.frame(species,genomesize)
    landscapeSpecies.df<-rbind(landscapeSpecies.df, landscapeSpecies.temp)
    rm(landscapeSpecies.temp)
  }
}

names(landscapes.df) <- gsub("X", "", names(landscapes.df), fixed = TRUE)

dbWriteTable(TEdb,name='landscapes',value=landscapes.df)
dbWriteTable(TEdb,name='landscapeSpecies',value=landscapeSpecies.df)

##Generate grid of landscape charts
landscape.graphs <- list()
i <- 1
setSessionTimeLimit(cpu = 99999, elapsed = 99999)

for (species in landscapeSpecies.df$species) {
  data <- fetch(dbSendQuery(TEdb, sprintf("SELECT * FROM landscapes WHERE species='%s' AND (Rclass='LTR' OR Rclass='SINE' OR Rclass='LINE' OR Rclass='DNA')", species)), n=-1)
  #rownames(data) <- data$Rclass
  data$row_names <- NULL
  data$species <- NULL
  #data$Rclass <- NULL

  long_data <-  gather(data, bin, value, 1:50)
  
  long_data <- long_data[!(apply(long_data[c("value")], 1, function(y) any(y == 0))),]
  long_data <- long_data[!(apply(long_data[c("bin")], 1, function(y) any(as.numeric(y) > 44))),]
  
  
  
  
  graph <- ggplot(long_data, aes(x=as.numeric(bin),y=value,fill=Rclass)) +
    geom_bar(stat='identity', width=1) + nogap +
    labs(title = species, x = "Kimura Distance", y = "Proportion of Genome")+ guides(fill=guide_legend(title=""))+
  theme(plot.title = element_text(size=8))

#ggtitle(sprintf("Species: %s", species)) +
  graph  
  break
  landscape.graphs[[i]] = graph

  i <- i+1

}


multiplot(plotlist = landscape.graphs, cols = 3)

################################ Box plot
data <- fetch(dbSendQuery(TEdb, sprintf("SELECT * FROM landscapes WHERE (Rclass='LTR' OR Rclass='SINE' OR Rclass='LINE' OR Rclass='DNA')", species)), n=-1)
data$row_names <- NULL
data$Rclass <- NULL
long_data <- gather(data, bin, value, 1:40)
box_data <- long_data %>%
  group_by(species, bin) %>%
  summarise(val=sum(value))


ggplot(box_data, aes(reorder(bin, as.numeric(bin)), val))+ geom_boxplot() + nogrid +vertText  + labs( x = "Kimura Distance", y = "Log(Proportion of Genome)") + geom_boxplot(fill = "white")


################################ Jitter plot
data <- fetch(dbSendQuery(TEdb, sprintf("SELECT * FROM landscapes WHERE (Rclass='LTR' OR Rclass='SINE' OR Rclass='LINE' OR Rclass='DNA')", species)), n=-1)
data$row_names <- NULL
data$species <- NULL
long_data <- gather(data, bin, value, 1:50)

long_data <- long_data[!(apply(long_data[c("value")], 1, function(y) any(as.numeric(y) == 0))),]
long_data <- long_data[!(apply(long_data[c("bin")], 1, function(y) any(as.numeric(y) > 40))),]

colnames(long_data)[colnames(long_data)=="Rclass"] <- "Class"

ggplot(long_data, aes(reorder(bin, as.numeric(bin)), log(value)))+  nogrid +vertText + labs( x = "Kimura Distance", y = "Log(Proportion of Genome") + geom_jitter(aes(colour = Class),size=1.5,width = 0.3) + geom_smooth(method=lm) + scale_y_continuous(labels = comma) + palette +
  geom_smooth(aes(x=as.integer(bin),y=log(value),color=Class,fill=Class),method=loess)



#+ 
#################### Stats tests
genomesize <- TE.length.comp.class <- fetch(dbSendQuery(TEdb, "SELECT genomesize FROM species"), n=-1)
genomesize <- genomesize$genomesize

t.test (genomesize, mu=mean(genomesize))


shapiro.test(genomesize)
qqnorm(genomesize)



### GEnome vs TE size
TE.comp <- fetch(dbSendQuery(TEdb, "SELECT species, genomesize, TEsize from species"), n=-1)
TE.comp <- melt(TE.length.comp.class, id.vars='species')

shapiro.test(TE.comp$genomesize)
shapiro.test(TE.comp$TEsize)

cor.test(TE.comp$genomesize, TE.comp$TEsize, 
         method = "pearson")

TE.comp <- fetch(dbSendQuery(TEdb, "SELECT species, genomesize, TEsize from species WHERE (species='Garboreum' OR species='Ghirsutum' OR species='Graimondii')"), n=-1)
shapiro.test(TE.comp$genomesize)
shapiro.test(TE.comp$TEsize)

cor.test(TE.comp$genomesize, TE.comp$TEsize, 
         method = "pearson")


TE.comp <- fetch(dbSendQuery(TEdb, "SELECT species, genomesize, TEsize from species WHERE (species='Ppersica' OR species='Fvesca' OR species='Tpratense' OR species='Gmax' OR species='Lusitatissimum' OR species='Rcommunis')"), n=-1)
shapiro.test(TE.comp$genomesize)
shapiro.test(TE.comp$TEsize)

cor.test(TE.comp$genomesize, TE.comp$TEsize, 
         method = "pearson")


###### genome vs te size for families

TE.class.df <- fetch(dbSendQuery(TEdb, "SELECT TEs.species, species.genomesize, SUM(TEs.LEN_MASKED), species.TEsize, TEs.Rclass FROM TEs, species WHERE species.species = TEs.species GROUP BY species,Rclass"), n=-1)
#TE.length.comp.class <- fetch(dbSendQuery(TEdb, "SELECT TEs.species, species.genomesize, species.TEsize, TEs.Rclass from TEs, species WHERE species.species = TEs.species AND TEs.Rclass='%s'",class.type)), n=-1)
names(TE.class.df) <- c("Species", "Genome", "TEs", "TEstotal", "Rclass")

TE.class <- TE.class.df

TE.class <- subset(TE.class, Rclass != "LTR  ")
#TE.class <- subset(TE.class, Rclass != "SINE  ")
TE.class <- subset(TE.class, Rclass != "LINE  ")

TE.class <- subset(TE.class, Rclass != "DNA  ")

TE.class$prop <- (TE.class$TEs/TE.class$Genome)*100
TE.class$propofTES <- (TE.class$TEs/TE.class$TEstotal)*100

shapiro.test(TE.class$Genome)
shapiro.test(TE.class$prop)

#TE.class <- melt(TE.class, id.vars='Species')
#t.test(value ~ variable, TE.class)

cor.test(log(TE.class$Genome), log(TE.class$prop), 
         method = "pearson")



################################################################################################################### Legacy stuff

#Parse the fasta file into a dataframe
dbSendQuery(TEdb, "drop table if exists te, species")
TEs.df <- data.frame(species=character(),name=character(),family=character(),status=character(),sequence=character(),stringsAsFactors=FALSE)
species.df <- data.frame(species=character(),genomesize=character())

parseFasta <- function(fspecies,fstatus) {
  
  fastaFile <- readDNAStringSet(paste(fspecies,".",fstatus,".fasta",sep=""))
  header = names(fastaFile)
  sequence = paste(fastaFile)
  species <- rep(fspecies, times = length(sequence))
  status <- rep(fstatus, times = length(sequence))
  headers <- strsplit(header,"#")
  name <- unlist(lapply(headers, `[[`, 1))
  family <- unlist(lapply(headers, `[[`, 2))
  data <- data.frame(species,name,family,status,sequence)
  return (data)
}

setwd("D:\\Documents\\School\\Current Courses\\BIOL 440\\Data\\Arabidopsis")

for (i in 1:length(speciesList)) {
  species = speciesList[i]
  
  genomesize = genomesizeList[i]
  species.df <- rbind(species.df,data.frame(species,genomesize))
  
  if(file.exists(paste(species,".","novel",".fasta",sep=""))) 
    TEs.df <- rbind(TEs.df,parseFasta(species,"novel"))
  
  if(file.exists(paste(species,".","known",".fasta",sep=""))) 
    TEs.df <- rbind(TEs.df,parseFasta(species,"known"))
}

dbWriteTable(TEdb,name='te',value=TEs.df)
dbWriteTable(TEdb,name='species',value=species.df)

#species <- "Osativa" ---- Species name for database
#genome.size <- 4500000000 ---- Total length of species genome
#mobilome.size <- 3100000 ---- Total length of mobilome





#parseRM 
#Show graph with the relative percentage of each class
TE.dist.class <- fetch(dbSendQuery(TEdb, "SELECT species, Rclass, SUM(LEN_MASKED_NR) FROM TEs GROUP BY species, Rclass"), n=-1)
names(TE.dist.class) <- c("Species", "Class", "Total.Length")
#Non-log
ggplot(data = TE.dist.class, aes(x = Species, y = Total.Length, fill = Class)) +  geom_bar(stat = "identity") + coord_flip()
#Log
ggplot(data = TE.dist.class, aes(x = Species, y = log(Total.Length), fill = Class)) +  geom_bar(stat = "identity") + coord_flip()

#Show graph with the relative percentage of each family
TE.dist.family <- fetch(dbSendQuery(TEdb, "SELECT species, CONCAT(`Rclass`,'.', `Rfam`), SUM(LEN_MASKED_NR) FROM TEs GROUP BY species, CONCAT(`Rclass`,'|', `Rfam`)"), n=-1)
names(TE.dist.family) <- c("Species", "Class.Family", "Total.Length")
#Non-log
ggplot(data = TE.dist.family, aes(x = Species, y = Total.Length, fill = Class.Family)) +  geom_bar(stat = "identity") + coord_flip()
#Log
ggplot(data = TE.dist.family, aes(x = Species, y = log(Total.Length), fill = Class.Family)) +  geom_bar(stat = "identity") + coord_flip()

#Show graph with the total length of TEs
TE.dist.total <- fetch(dbSendQuery(TEdb, "SELECT species, SUM(LEN_MASKED_NR) FROM TEs GROUP BY species"), n=-1)
names(TE.dist.class) <- c("Species", "Total.Length")
#Non-log
ggplot(data = TE.dist.class, aes(x = Species, y = Total.Length)) +  geom_bar(stat = "identity") + coord_flip()
#Log
ggplot(data = TE.dist.class, aes(x = Species, y = log(Total.Length))) +  geom_bar(stat = "identity") + coord_flip()


