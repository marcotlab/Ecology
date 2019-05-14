require(compiler)
# require(MASS)

getSpecimenMatFromMeasurements <- function(filename="~/Dropbox/code/R/amandaTeeth/dat/specimens.csv") {
	require(abind)
	# datMeasures<-read.csv("http://dl.dropbox.com/s/6fwrim5z29lcv58/specimens.csv", strip.white=TRUE)
	datMeasures <- read.csv(filename, strip.white=TRUE)
	# datMeasures$species<-synonymize(as.character(datMeasures$species))
	
	m <- c("P2_L","P2_W","P3_L","P3_W","P4_L","P4_W","M1_L","M1_W","M2_L","M2_W","M3_L","M3_W","p2_l","p2_w","p3_l","p3_w","p4_l","p4_w","m1_l","m1_w","m2_l","m2_w","m3_l","m3_w","upP","upM","loP","loM")
	cube <- abind::abind(datMeasures[seq(from=1, to=nrow(datMeasures), by=3),m], datMeasures[seq(from=2, to=nrow(datMeasures), by=3),m], along=3)
	cube <- abind::abind(cube, datMeasures[seq(from=3, to=nrow(datMeasures), by=3),m], along=3)
	replicateMeasureMat <- apply(cube, c(1,2), mean, na.rm=TRUE)
	rownames(replicateMeasureMat) <- datMeasures$species[seq(from=1, to=nrow(datMeasures), by=3)]
	replicateMeasureMat[!is.finite(replicateMeasureMat)] <- NA
	# replicateMeasureMat<-data.frame(species=datMeasures$species[seq(from=1, to=nrow(datMeasures), by=3)], specimen=datMeasures$Specimen.no.[seq(from=1, to=nrow(datMeasures), by=3)], replicateMeasureMat, stringsAsFactors=FALSE, row.names=NULL)
	
	theseSpecimenNos <- datMeasures$Specimen.no.[seq(from=1, to=nrow(datMeasures), by=3)]
	specimenNos <- unique(theseSpecimenNos)
	specimenMat <- matrix(nrow=0, ncol=ncol(replicateMeasureMat))
	for (i in seq_len(length(specimenNos))) {
		index <- which(theseSpecimenNos==specimenNos[i])
		thisSpecimen <- replicateMeasureMat[index,]
		if (!is.null(dim(thisSpecimen))) thisSpecimen <- apply(thisSpecimen, 2, mean, na.rm=TRUE)
		specimenMat <- rbind(specimenMat, thisSpecimen)
		rownames(specimenMat)[nrow(specimenMat)] <- rownames(replicateMeasureMat)[index[1]]
	}
	specimenMat <- data.frame(species=rownames(specimenMat), specimen=specimenNos, specimenMat, stringsAsFactors=FALSE, row.names=NULL)
	specimenMat
}

getBlastoSpecimenMat <- function() {
	lengths<-read.csv("~/Dropbox/code/R/blasto/dat/lengths.csv", strip.white=TRUE)
	widths<-read.csv("~/Dropbox/code/R/blasto/dat/widths.csv", strip.white=TRUE)
	l_mean<-data.frame("specimen"=lengths[,1], "p2_l"=apply(lengths[,grep("p2", colnames(lengths))], 1, mean, na.rm=TRUE), "p3_l"=apply(lengths[,grep("p3", colnames(lengths))], 1, mean, na.rm=TRUE), "p4_l"=apply(lengths[,grep("p4", colnames(lengths))], 1, mean, na.rm=TRUE), "m1_l"=apply(lengths[,grep("m1", colnames(lengths))], 1, mean, na.rm=TRUE), "m2_l"=apply(lengths[,grep("m2", colnames(lengths))], 1, mean, na.rm=TRUE), "m3_l"=apply(lengths[,grep("m3", colnames(lengths))], 1, mean, na.rm=TRUE))
	w_mean<-data.frame("specimen"=widths[,1], "p2_w"=apply(widths[,grep("p2", colnames(widths))], 1, mean, na.rm=TRUE), "p3_w"=apply(widths[,grep("p3", colnames(widths))], 1, mean, na.rm=TRUE), "p4_w"=apply(widths[,grep("p4", colnames(widths))], 1, mean, na.rm=TRUE), "m1_w"=apply(widths[,grep("m1", colnames(widths))], 1, mean, na.rm=TRUE), "m2_w"=apply(widths[,grep("m2", colnames(widths))], 1, mean, na.rm=TRUE), "m3_w"=apply(widths[,grep("m3", colnames(widths))], 1, mean, na.rm=TRUE))
	dat <- merge(l_mean, w_mean, by="specimen")
	dat[!sapply(dat, is.finite)]<-NA
	dat <- merge(read.csv("~/Dropbox/code/R/blasto/dat/info.csv", strip.white=TRUE), dat, all=TRUE)
	specimenMat <- data.frame(species=dat$sp_current, dat[,c(1,13:ncol(dat))])
	specimenMat
}

getLiteratureSpecimenMat <- function(filename="~/Dropbox/code/R/amandaTeeth/dat/literature.csv") {
	dat <- read.csv(filename, strip.white=TRUE)
	# dat <- dat[apply(is.finite(data.matrix(dat[,3:ncol(dat)])), 1, any),] # removes taxa that with no measurements
	dat
}

makeOneSpeciesMatFromSpecimenMat <- function (specimenMat) {
	# species<-unique(specimenMat$species)
	oneSpeciesMat <- aggregate(specimenMat, by=list(specimenMat$species), mean, na.rm=TRUE)
	# oneSpeciesMat <- aggregate(specimenMat, by=list(specimenMat$species), median, na.rm=TRUE)
	oneSpeciesMat <- data.frame(oneSpeciesMat, row.names=oneSpeciesMat[,1])
	oneSpeciesMat[sapply(oneSpeciesMat,is.nan)] <- NA
	colnames(oneSpeciesMat)[1] <- "species"
	oneSpeciesMat[,-(c(2:3, which(colnames(oneSpeciesMat) %in% c("Published.name",  "n", "Reference", "PaleoDB.ref", "Notes",  "X", "X.1", "X.2", "X.3", "X.4"))))]
}

getToothRowLengths<-function(species) {
	thisRow<-vector(mode="numeric", length=0)
	if (!is.nan(species$upP)) { thisRow<-c(thisRow, upP=species$upP)
	} else {
		thisList<-c(species$P2_L, species$P3_L, species$P4_L)
		if (!any(is.nan(thisList))) thisRow<-c(thisRow, upP=sum(thisList)) else thisRow<-c(thisRow, upP=NA)
	}

	if (!is.nan(species$upM)) { thisRow<-c(thisRow, upM=species$upM)
	} else {
		thisList<-c(species$M1_L, species$M2_L, species$M3_L)
		if (!any(is.nan(thisList))) thisRow<-c(thisRow, upM=sum(thisList)) else thisRow<-c(thisRow, upM=NA)
	}

	if (!is.nan(species$loP)) { thisRow<-c(thisRow, loP=species$loP)
	} else {
		thisList<-c(species$p2_l, species$p3_l, species$p4_l)
		if (!any(is.nan(thisList))) thisRow<-c(thisRow, loP=sum(thisList)) else thisRow<-c(thisRow, loP=NA)
	}

	if (!is.nan(species$loM)) { thisRow<-c(thisRow, loM=species$loM)
	} else {
		thisList<-c(species$m1_l, species$m2_l, species$m3_l)
		if (!any(is.nan(thisList))) thisRow<-c(thisRow, loM=sum(thisList)) else thisRow<-c(thisRow, loM=NA)
	}
	return(thisRow)
}

getBMEstimateFromMeasurementAndRegression <- function(m, regRow) {
	10^((log10(m) * regRow$slope) + regRow$intercept)
}
getBMEstimateFromMeasurementAndRegression_compiled <- cmpfun(getBMEstimateFromMeasurementAndRegression)

getBodyMassFromMeasurement <- function(thisCol, regRow) {
	sapply(thisCol, getBMEstimateFromMeasurementAndRegression_compiled, regRow)
}

# meanUsingReg<-function(species, reg) {
	# valueList<-vector(mode="numeric", length=0)
	# weightList<-vector(mode="numeric", length=0)
	# for (i in seq_len(length(species))) {
		# if (!is.na(species[i])) {
			# r2 <- reg$r2[which(reg$m==names(species[i]))]
			# valueList <- c(valueList, species[i]*r2)
			# weightList <- c(weightList, r2)
		# }
	# }
	# sum(valueList)/sum(weightList)
# }
# cmur<-cmpfun(meanUsingReg)

getMendozaBodyMasses <- function (this) {
	sevenOne <- (2.543 * this$m1_w) + (1.827 * this$loM) + (1.402 * this$m2_w) + 1.280 #7.1
	sevenTwo <- (1.850 * this$m1_w) + (1.883 * this$loM) + (1.073 * this$m2_w) + (0.308 * this$p4_w) + 1.216 #7.2
	sevenThree <- (1.460 * this$m1_w) + (1.363 * this$loM) + (1.182 * this$m2_w) + (0.387 * this$p4_w) + (0.955 * this$m2_l) + 1.578 #7.3
	sevenFour <- (1.226 * this$m1_w) + (1.313 * this$loM) + (1.090 * this$m2_w) + (0.320 * this$p4_w) + (1.095 * this$m2_l) + (0.150 * this$loP) + 1.361 #7.4
	sevenFive <- (1.355 * this$m1_w) + (1.427 * this$loM) + (1.322 * this$m2_w) + (0.457 * this$p4_w) + (1.177 * this$m2_l) + (0.234 * this$loP) + (0.340 * this$p4_l) + 1.168 #7.5
	sevenSix <- (2.045 * this$m1_l) + (1.073 * this$loM) + 1.507	#7.6*
	sevenSeven <- (1.213 * this$m1_l) + (1.421 * this$loM) + (0.422 * this$p4_w) + 1.380	#7.7*
	c(sevenOne, sevenTwo, sevenThree, sevenFour, sevenFive, sevenSix, sevenSeven)
}

appendStDevToReg <- function(thisReg) {
	cbind(thisReg, array((log10(100 + thisReg["see"]) - 2), dimnames=list(rownames(thisReg), "stdev")))
	# cbind(thisReg, array((log10(100 + thisReg["see"]) - 2) * sqrt(thisReg["n"]), dimnames=list(rownames(thisReg), "stdev")))
	# thisReg["see"]=(2+see)	
}

getLikelihoodOneBodyMassFromMeasurementAndReg <- function(bm, this, thisReg) {
	lnl <- vector()
	for (i in seq_along(this)) lnl[i] <- dnorm(bm, mean = (log10(this[i]) * thisReg$slope[thisReg$m == colnames(this)[i]]) + thisReg$intercept[thisReg$m == colnames(this)[i]], sd=thisReg$stdev[thisReg$m == colnames(this)[i]], log=TRUE)
	sum(-lnl)
}

getMLbodyMassForOneSpecies <- function(this, thisReg, best.only=FALSE) {
	if (best.only) thisReg <- thisReg[thisReg$j %in% c("FLML","FLMA","SLML","SLMA","SUML","SUMA","TLMA","LMRL") ,]		#Janis's "best" variables Table 13.3, p.281
	this <- this[sapply(this, is.finite) & names(this) %in% thisReg$m]
	if (ncol(this)==0) return (NA)
	# bmVec <- vector()
	# for (i in seq_along(this)) bmVec <- c(bmVec, getBMEstimateFromMeasurementAndRegression_compiled(this[i], thisReg[match(colnames(this)[i], thisReg$m),]))
	# # sapply(this, getBMEstimateFromMeasurementAndRegression, thisReg[match(names(this[i]), as.character(thisReg$m)),])
	# mean(bmVec, na.rm=TRUE)
	optim(2, fn=getLikelihoodOneBodyMassFromMeasurementAndReg, this=data.matrix(this), thisReg=thisReg, method="L-BFGS-B", lower=-Inf, upper=Inf)$par
}
getMLbodyMassForOneSpecies_compiled <- cmpfun(getMLbodyMassForOneSpecies)

getMLBodyMasses <- function(thisMat, regList, best.only=FALSE) {
	bmVec <- array(NA, dim=c(nrow(thisMat), 1), dimnames=list(rownames(thisMat), "bodyMass"))
	for (i in seq_len(nrow(thisMat))) {
		if (!is.na(thisMat$reg[i])) {
			bmVec[i] <- getMLbodyMassForOneSpecies_compiled(this=thisMat[i,], thisReg=regList[[which(names(regList) == thisMat$reg[i])]])
		} else bmVec[i] <- NA
	} 
	bmVec
} 
getMLBodyMasses_compiled <- cmpfun(getMLBodyMasses)

# getBodyMassVectorFromThisMat <- function(thisMat) {
	# famList <- unique(read.csv("~/Dropbox/code/common_dat/taxonomy.csv"), strip.white=TRUE)
	# regList <- list(ruminantia=read.csv("~/Dropbox/code/R/amandaTeeth/dat/regRuminantia.csv"), perissodactyla=read.csv("~/Dropbox/code/R/amandaTeeth/dat/regRuminantia.csv"), ungulate=read.csv("~/Dropbox/code/R/amandaTeeth/dat/regAllUngulates.csv"))
	# regList <- lapply(regList, appendStDevToReg)
	# thisMat[,"reg"] <- famList$reg[match(gsub(pattern="\"", replacement="", x=thisMat$species), famList$taxon)]  #clears any quotes from species names
	# thisMat$reg[thisMat$reg==""] <- NA
	# getMLBodyMasses_compiled(thisMat, regList, best.only=FALSE)
# }

appendRegTypeToThisMat <- function(thisMat) {
		famList <- unique(read.csv("~/Dropbox/code/common_dat/taxonomy.csv"), strip.white=TRUE)
		thisMat[,"reg"] <- famList$reg[match(x=thisMat$species, famList$taxon)]	# now assuming reg will already be a part of thisMat
		thisMat
}

getBodyMassVectorFromThisMatAllMeasures <- function(thisMat, linked.files=FALSE) {
	#######################################################################################################################################
	##### read Janis 1990 regression parameters from file, and append standard deviations
	#######################################################################################################################################
	if (linked.files) {
		# famList <- unique(read.csv("https://dl.dropbox.com/s/wrzovc89oqispyk/taxonomy.csv"), strip.white=TRUE)
		regList <- list(ruminantia=read.csv("https://dl.dropbox.com/s/dcd0bs1x5v9e7lh/regRuminantia.csv"), perissodactyla=read.csv("https://dl.dropbox.com/s/04k387q7yh4wp9u/regPerissodactyla.csv"), ungulate=read.csv("https://dl.dropbox.com/s/310ayur1s1dc8sl/regAllUngulates.csv"))
	} else {
		regList <- list(ruminantia=read.csv("~/Dropbox/code/R/amandaTeeth/dat/regRuminantia.csv"), perissodactyla=read.csv("~/Dropbox/code/R/amandaTeeth/dat/regPerissodactyla.csv"), ungulate=read.csv("~/Dropbox/code/R/amandaTeeth/dat/regAllUngulates.csv"))
	}
	regList <- lapply(regList, appendStDevToReg)

	#######################################################################################################################################
	##### get body mass for only those taxa that have all (i.e., are not missing any) of the published measurements (about 325 species)
	#######################################################################################################################################
	thisMat$reg[thisMat$reg==""] <- NA
	theseColumns <- c(as.character(regList[[1]]$m)[-which(as.character(regList[[1]]$m) %in% c("loP", "loM"))], "reg") # theseColumns is the names of columns that are also in reg - published measuremnts
	bm <- getMLBodyMasses_compiled(thisMat[complete.cases(thisMat[,theseColumns]), theseColumns], regList, best.only=FALSE)
	
	#######################################################################################################################################
	##### get regression parameters for measurements not in published regression
	#######################################################################################################################################
	other_m <- c("P2_L", "P2_W", "P3_L", "P3_W", "P4_L", "P4_W", "M1_L", "M1_W", "M3_L", "M3_W")
# 	other_m <- colnames(thisMat[,sapply(thisMat, is.numeric)])[!colnames(thisMat[,sapply(thisMat, is.numeric)]) %in% theseColumns]

	otherReg <- lapply(X=names(regList), FUN=function(this.group) {
		shortMat <- thisMat[rownames(thisMat) %in% rownames(bm) & thisMat$reg==this.group, other_m]
		short.bm <- bm[rownames(bm) %in% rownames(shortMat),]
		apply(log10(shortMat), MARGIN=2, FUN=function(x, bm) {
			lm(bm ~ x) } , bm=short.bm ) 
		} )

	names(otherReg) <- names(regList)

	otherList <- lapply(otherReg, function(thisReg) t(sapply(thisReg, function(x) c(coef(x), summary(x)$sigma))))
	
	#######################################################################################################################################
	##### merge regression parameters for unpublished measurements (otherList) with those from published (regList)
	#######################################################################################################################################
	for (i in seq_along(regList)) {
		colnames(otherList[[i]]) <- c("intercept", "slope", "stdev")
		otherList[[i]] <- data.frame(m=rownames(otherList[[i]]), otherList[[i]], stringsAsFactors=FALSE)
		regList[[i]] <- merge(regList[[i]], otherList[[i]], all=TRUE, sort=FALSE)
	}

	#######################################################################################################################################
	##### recalculate body masses of all taxa with all (merged) measurements
	#######################################################################################################################################
	bm <-  getMLBodyMasses_compiled(thisMat, regList, best.only=FALSE)
	bm[match(thisMat$species, rownames(bm))]	
}

appendMissingPaleoDBSpecies <- function(thisMat, tax.vec) {
	# this adds taxa that are in PaleoDB (i.e., occurrence data), but not in the measurement files
	tax.vec <- tax.vec[!tax.vec %in% thisMat$species]
	tax.frame <- data.frame(array(NA, dim=c(length(tax.vec), ncol(thisMat)), dimnames=list(tax.vec, colnames(thisMat))))
	tax.frame$species <- tax.vec
	if (any(tax.frame$species %in% thisMat$species)) { 
		# thisMat[match(rangesNotThisMat$species, thisMat$species),c("FO", "LO")] <- ranges[,c("FO", "LO")]
	} else thisMat <- merge(thisMat, tax.frame, all=TRUE, sort=FALSE)
	rownames(thisMat) <- thisMat$species
	thisMat
}

fillMissingBodyMasses <- function(thisMat) {
	require(stringr)
	noMass <- rownames(thisMat)[!is.finite(thisMat$bodyMass)]
	noMass <- data.frame(species=noMass, bodyMass=sapply(X=noMass, FUN=function(x) mean(thisMat[grep(str_split(x, pattern=" ")[[1]][1], rownames(thisMat)),"bodyMass"], na.rm=TRUE)))
	thisMat$bodyMass[match(rownames(noMass), rownames(thisMat))] <- noMass[,"bodyMass"]
	array(thisMat$bodyMass, dimnames=list(rownames(thisMat)))
}

