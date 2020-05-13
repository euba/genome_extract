library(seqinr) #load in library

#load in assembly in R
assembly <- read.fasta("Downloads/Ace_trimmed_contigs_spades.FASTA/Ace_trimmed_contigs_spades.fa")

#calculate length, gc, and coverage statistics for assembled contigs
lengths <- unlist(lapply(assembly,length)) #calculate the length for each contig
gcs <- unlist(lapply (assembly,GC)) #calculate the GC for each contig
covs <- lapply(strsplit(names(assembly),split="_"),function(x){x[length(x)]}) #get coverage from the name of each contig
covs <- as.numeric(unlist(covs))
names(covs) = names(assembly)

#plot the gc and coverage against each other
plot(gcs,log10(covs),type="n")
points(gcs,log10(covs),pch=16,cex=0.5)
#plot the gc and length against each other
plot(gcs,lengths,type="n")
points(gcs,lengths,pch=16,cex=0.5)

#extract specific contigs from assembly based on gc
select <- assembly[names(which(gcs<0.25))] #select contigs with gcs less than 25%
#extract specific contigs from assembly based on gc and length
select <- assembly[intersect(names(which(gcs<0.25)),names(which(lengths>1000)))] #select contigs with gcs less than 25% and contigs longer than 1000bp
#extract specific contigs from assembly based on gc and coverage
select <- assembly[intersect(names(which(gcs<0.25)),names(which(covs>200)))] #select contigs with gcs less than 25% and a coverage higher than 200

#write the selected contigs as a separate file
write.fasta(select,file.out="selected_contigs.")
