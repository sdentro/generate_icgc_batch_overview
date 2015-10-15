d = read.table("input/santa_cruz_pilot.v2.2015_0504.tsv", header=T, stringsAsFactors=F, sep="\t")
b6 = read.table("input/batch06.txt", stringsAsFactors=F)
b7 = read.table("input/batch07.txt", stringsAsFactors=F)
b8 = read.table("input/batch08.txt", stringsAsFactors=F)

getColumns = function(d) {
	# Split sample ids into two lists for single and multiples
	samplelist = d$Tumour.WGS.aliquot.ID.s.
	singles = samplelist[-grep(",", samplelist)]
	multiples = samplelist[grep(",", samplelist)]
	multiples = unlist(lapply(multiples, function(x) { unlist(strsplit(x, ",")) }))
	allsamples = c(singles, multiples)

	# Fetch the donor ids
	multiples_rowids = unlist(lapply(multiples, function(x) { grep(x, d$Tumour.WGS.aliquot.ID.s.) }))
	donorids = c(d$Submitter.donor.ID[-grep(",", samplelist)], d$Submitter.donor.ID[multiples_rowids])

	# Fetch the project code
	projectcodes = c(d$Project.code[-grep(",", samplelist)], d$Project.code[multiples_rowids])

	return(list(samplelist=samplelist, singles=singles, multiples=multiples, allsamples=allsamples, multiples_rowids=multiples_rowids, donorids=donorids, projectcodes=projectcodes))
}

singleCols = getColumns(d)
samplelist = singleCols$samplelist
allsamples = singleCols$allsamples
multiples_rowids = singleCols$multiples_rowids
donorids = singleCols$donorids
projectcodes = singleCols$projectcodes

output = data.frame(projectcode=projectcodes, 
			donorid=donorids, 
			sampleid=allsamples, 
			batch06=allsamples %in% b6[,1], 
			batch07=allsamples %in% b7[,1], 
			batch08=allsamples %in% b8[,1], 
			qc_failed=!((allsamples %in% b6[,1]) | (allsamples %in% b7[,1]) | (allsamples %in% b8[,1])))

# Samples that were not sent out, yet marked as pass
not_sent_out = read.table("input/not_sent_out.lst", header=F)
for (samplename in not_sent_out[,1]) {
	output[output$sampleid==samplename,4:7] = c(F,F,F,T)
}

# Samples that were shared but are marked as fail - all in batch07
shared_eroneously = read.table("input/shared_eroneously.lst", header=F)
for (samplename in shared_eroneously[,1]) {
	output[output$sampleid==samplename,4:7] = c(F,T,F,F)
}

# Samples that were shared more than once are assigned to the first batch in which they were shared
shared_more_than_once = read.table("input/shared_more_than_once.lst", header=F)
for (samplename in shared_more_than_once[,1]) {
	first = head(which(head(output[output$sampleid==samplename,4:7], 1)==TRUE), 1)+3
	output[output$sampleid==samplename,4:7] = c(F,F,F,F)
	output[output$sampleid==samplename,first] = TRUE
}

write.table(output, file="santa_cruz_pilot_battenberg_batches_wide.tsv", quote=F, row.names=F, sep="\t")


meltData = function(output) {
	library(reshape2)
	dat.m = melt(output, id.vars=c("projectcode","donorid","sampleid"))
	dat.m = dat.m[dat.m$value,]
	dat.m = dat.m[,1:4]
	colnames(dat.m)[4] = "release"
	return(dat.m)
}

dat.m = meltData(output)
write.table(dat.m, file="santa_cruz_pilot_battenberg_batches.tsv", quote=F, row.names=F, sep="\t")

# Add in the August release new samples
august = read.table("input/release_aug2015.v1.new.tsv", header=T, stringsAsFactors=F, sep="\t")
# Setting the column names to what they were before in the Santa Cruz pilot
# Note that this works but columns 3 and 4 have changed usage between the two.
# 3 was Data.train and is now a boolean representing in santa_cruz_pilot
# 4 was a boolean representing Train2.pilot and is now validation_by_deep_seq
colnames(august) = colnames(d)
singleCols = getColumns(august)
samplelist = singleCols$samplelist
allsamples = singleCols$allsamples
donorids = singleCols$donorids
projectcodes = singleCols$projectcodes

output$in_qc = rep(F, nrow(output))
output = rbind(output, data.frame(projectcode=projectcodes,
                                  donorid=donorids,
                                  sampleid=allsamples,
                                  batch06=rep(F, length(projectcodes)),
                                  batch07=rep(F, length(projectcodes)),
                                  batch08=rep(F, length(projectcodes)),
                                  qc_failed=rep(F, length(projectcodes)),
                                  in_qc=rep(T, length(projectcodes))))
dat.m = meltData(output)
write.table(output, file="icgc_battenberg_batches_wide.tsv", quote=F, row.names=F, sep="\t")
write.table(dat.m, file="icgc_battenberg_batches.tsv", quote=F, row.names=F, sep="\t")
