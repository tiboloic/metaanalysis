# LT 11/06/2015

# library for meta analyses
library(metafor)

# regionalization method used. one of:
# c("sample","fn100id", "fn500id", "ECOREGION", "PROVINCE", "REALM", "EEZ") 

regionalization = "EEZ"

# file name for coverage standardized genetic diversity data
stats_file = paste("DIPnet_stats_061015", regionalization, ".Rdata", sep="")

# file name for spatial data (species richess per regionalization)
spatial_file = "ipdb_sp.tsv"

# read spatial file and add a column with "sample" population ID (gsl_lat_log)
spatial = read.table(spatial_file, header=T, sep="\t",stringsAsFactors = F, na.strings=c("NA"," ",""), quote="")

# define population as same species_locus, same locality, same rounded lat and long
spatial$sample = paste(spatial$locality,round(spatial$decimalLat, digits=0),round(spatial$decimalLon, digits=0),sep="_")

# some locations have hard quotes in them that we need to remove to match Cynthia's names
spatial$sample = gsub("\"", "", spatial$sample)

# same thing applies to PROVINCE regionalization
spatial$PROVINCE = gsub("\"", "", spatial$PROVINCE)

# EEZ Réunion is misspelled
spatial$EEZ[spatial$EEZ == "R?®union Exclusive Economic Zone"] = "Reunion Exclusive Economic Zone"

# EEZ Curaçaoan EEZ
spatial$EEZ[spatial$EEZ == "Cura?ºaoan Exclusive Economic Zone"] = "Curacaoan Exclusive Economic Zone"


# aggregate spatial data to get one measure of species diversity per location
species_richness = aggregate(spatial$total_rich, by=list(spatial[,regionalization]), FUN=mean)
colnames(species_richness) = c("location","total_rich")

# load genetic diversity data
load(stats_file)

# data structure to store the slopes and associated standard errors
meta = c()

# loop on gsl

#for (gsl in divstats) {
for (igsl in 1:length(divstats)){

	gsl = divstats[[igsl]]

	# set location, stored in rownames, as a column 
	gsl$location = rownames(gsl)
	
	# clean up hard quotes in location names
	gsl$location = gsub("\"", "", gsl$location)

	# some extra spaces in Cynthia's location names break the matching, remove them
	# example: "Caohagan, Lapu-Lapu, Philippines _10_124" vs "Caohagan, Lapu-Lapu, Philippines_10_124"
	gsl$location = gsub(" _", "_", gsl$location)
	
	# get standard error of coverage standardized diversity
	gsl$se = (gsl$qD - gsl$qD.95.LCL) /qnorm(.975)

	# keep only populations for which standard error is positive 
	gsl = subset(gsl, gsl$se>0)

	# merge species richness from spatial file
	new_gsl = merge(gsl, species_richness)
	
	if (nrow(new_gsl) < nrow(gsl)) cat(paste("Some locations missing ", igsl, "\n"))
	gsl=new_gsl

	# remove populations with species richness = 0
	new_gsl = subset(gsl, total_rich >0)
	if (nrow(new_gsl) < nrow(gsl)) cat(paste("Some locations have species richness 0", igsl, "\n"))
	gsl=new_gsl

	# test if we have variance in species richness
	if (length(unique(gsl$total_rich)) == 1) {cat(paste("Same spp richness for all pops ", igsl, "\n")); next;}

	# perform analysis only when we have more than 3 populations
	if (nrow(gsl)<3) {cat(paste("Less than 3 populations ", igsl, "\n"));next;}

	# perform meta-analysis regression of genetic diversity vs species richness
	fit = rma.uni(qD, sei=se, data=gsl, mods=~total_rich)

	# save slope ans associated standard error 
	meta = rbind(meta, c(fit$b[2,1], se=fit$se[2]))
}

# perform final meta analysis on slopes
fit.final = rma.uni(total_rich, sei=se, data=meta)
summary(fit.final)

save(meta, fit.final, file=paste("Meta.", regionalization, ".Rdata", sep=""))