# LT 15/06/2015

# Quick figure for restitution of preliminary meta analysis results
# genetic diversity over total richness

regions = c("sample","fn100d", "fn500d", "ECOREGION", "PROVINCE", "REALM", "EEZ")
effects=c()

for (region in regions) {

	# load meta analysis results
	load(paste("Meta.", region, ".Rdata", sep=""))

	# store effect size and standard error
	effects = rbind (effects,  c(slope=fit.final$b, se=fit.final$se))
}
rownames(effects) = regions

# plot 

plot(effects[,"slope"], ylim=c(-0.01,0.05), xaxt="n", xlab="regionalization", ylab="slope genetic diversity over total richness")
axis(1, at=1:7, labels=regions)
arrows(x0=1:7, y0=effects[,"slope"], y1=effects[,"slope"]+ effects[,"se"], length=0.1, angle=90)
arrows(x0=1:7, y0=effects[,"slope"], y1=effects[,"slope"]- effects[,"se"], length=0.1, angle=90)
abline(h=0, lty="dashed")