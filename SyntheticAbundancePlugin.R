###########################################################################
# Function: Abundance table simulation and sequencing reads simulation
# Call: Rscript simulation_abd.R type[gut/oral/skin/buliding/ocean/vaginal] [sample#] [min_species#] [max_species#] [distribution_type] [enable_sequencing_reads_generation][T/F] [reads#]
# Example: Rscript simulation_abd.R gut 100 100 200 L T 100000
# Explanation: Generating an abundance table and sequencing data from zero-inflated log-normal distribution with 100 samples, species number is randomly selected range from 100 to 200 in each sample.
# Genome_ref and Wgsim are necessary for the script
# Species number included in database:
# ----------------------------------------------------------------------------------
#	Gut	Oral	Skin	Vaginal	Ocean	human	building
# ----------------------------------------------------------------------------------
# 	409	445	272	158	109	1287	349
# ----------------------------------------------------------------------------------
# For more information about GCF accesion number, please refers to GCF_annotation.
# L: lognormal distribution, R: random distribution, N: one-side normal distribution, E:even distribution.
# Last update: 2020-09-02, Zheng Sun 
###########################################################################

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")
source("RIO.R")



input <- function(inputfile) {
   parameters <<- readParameters(inputfile)
}

run <- function() {}

output <- function(outputfile) {
pfix <- prefix()
	args1 <- paste(pfix, parameters["input", 2], sep="/")
args2 <- parameters["numsamples", 2]
args3 <- parameters["minspecies", 2]
args4 <- parameters["maxsample", 2]
args5 <- parameters["distributiontype", 2]


#Args assignment
genome_length <- read.table(args1,header = T,sep="\t")
sample_num <- args2
min_spp <- args3 #assigning the minimum nubmer of spp in simulation data
max_spp <- args4 #assigning the maximum number of spp in simulation data
max_species_db <- dim(genome_length)[1] #calculate the maximum number of spp in database
smo_taxonabd <- cbind() # original simulated matrix (taxonomic abd)
smo_seqabd <- cbind() # original simulated matrix (sequence abd)

if (args5=="L") { #log-normal distribution
for(i in 1:sample_num) { # indicating how many samples should be simulated 
	species_num<- sample(min_spp:max_spp,1) #generating 100-199 random integer
	Sample_i <- c(rlnorm(species_num),rep(0, times=(max_species_db-species_num))) #generaing data from zero-inflated lognormal distribution
	Sample_r <- Sample_i/sum(Sample_i) #normalization
	Sample_R <- sample(Sample_r,length(Sample_r))#mess up the order
#sequence abundance
	Sample_i_seq <- Sample_R*genome_length[,2] #calculating sequencing abd
	Sample_r_seq <- Sample_i_seq/sum(Sample_i_seq) #normalization
#comnination abundance
	smo_taxonabd <- cbind(smo_taxonabd,Sample_R) #combine the taxonomic abd
	smo_seqabd <- cbind(smo_seqabd,Sample_r_seq) #combine the sequence abd
}
} else if (args5=="N") { #one-side normal distribution
for(i in 1:sample_num) { # indicating how many samples should be simulated 
	species_num<- sample(min_spp:max_spp,1) #generating 100-199 random integer
	Sample_i <- c(abs(rnorm(species_num)),rep(0, times=(max_species_db-species_num))) #generaing data from zero-inflated lognormal distribution
	Sample_r <- Sample_i/sum(Sample_i) #normalization
	Sample_R <- sample(Sample_r,length(Sample_r))#mess up the order
#sequence abundance
	Sample_i_seq <- Sample_R*genome_length[,2] #calculating sequencing abd
	Sample_r_seq <- Sample_i_seq/sum(Sample_i_seq) #normalization
#comnination abundance
	smo_taxonabd <- cbind(smo_taxonabd,Sample_R) #combine the taxonomic abd
	smo_seqabd <- cbind(smo_seqabd,Sample_r_seq) #combine the sequence abd
}
} else if (args5=="R") { #random distribution
for(i in 1:sample_num) { # indicating how many samples should be simulated 
	species_num<- sample(min_spp:max_spp,1) #generating 100-199 random integer
	Sample_i <- c(runif(species_num),rep(0, times=(max_species_db-species_num))) #generaing data from zero-inflated lognormal distribution
	Sample_r <- Sample_i/sum(Sample_i) #normalization
	Sample_R <- sample(Sample_r,length(Sample_r))#mess up the order
#sequence abundance
	Sample_i_seq <- Sample_R*genome_length[,2] #calculating sequencing abd
	Sample_r_seq <- Sample_i_seq/sum(Sample_i_seq) #normalization
#comnination abundance
	smo_taxonabd <- cbind(smo_taxonabd,Sample_R) #combine the taxonomic abd
	smo_seqabd <- cbind(smo_seqabd,Sample_r_seq) #combine the sequence abd
}
} else if (args5=="E") { #even distribution
for(i in 1:sample_num) { # indicating how many samples should be simulated 
	species_num<- sample(min_spp:max_spp,1) #generating 100-199 random integer
	Sample_i <- c(rep(1, times=(species_num)),rep(0, times=(max_species_db-species_num))) #generaing data from zero-inflated lognormal distribution
	Sample_r <- Sample_i/sum(Sample_i) #normalization
	Sample_R <- sample(Sample_r,length(Sample_r))#mess up the order
#sequence abundance
	Sample_i_seq <- Sample_R*genome_length[,2] #calculating sequencing abd
	Sample_r_seq <- Sample_i_seq/sum(Sample_i_seq) #normalization
#comnination abundance
	smo_taxonabd <- cbind(smo_taxonabd,Sample_R) #combine the taxonomic abd
	smo_seqabd <- cbind(smo_seqabd,Sample_r_seq) #combine the sequence abd
}
} else {
stop("Please assign a type of distribution for abundance simulation.")
}


#generating tables
species_ID <- c("SpeciesID",as.character(genome_length[,3])) #factor to character
smo_taxonabd <- rbind(paste("sample",seq(from=1,to=args2),sep="_"),smo_taxonabd) #header added
smo_taxonabd <- cbind(species_ID,smo_taxonabd) #rownames added
smo_seqabd <- rbind(paste("sample",seq(from=1,to=args2),sep="_"),smo_seqabd ) #header added
smo_seqabd <- cbind(species_ID,smo_seqabd ) #rownames added

write.table(smo_taxonabd,file=paste(outputfile,"_taxonomic_abd.txt",sep=""),quote=F,,row.names=F,col.names = F,sep = "\t")
write.table(smo_seqabd,file=paste(outputfile,"_sequence_abd.txt",sep=""),quote=F,,row.names=F,col.names = F,sep = "\t")


}
