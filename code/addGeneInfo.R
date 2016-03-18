
##dir of differentially expressed genes
degDir="/home/zxh/Projects/jkCheng_shanghaiJiaoda/genomicSeq/results/Cuffdiff/Unorm"
##file types
ft=c("gene_exp.diff", "isoform_exp.diff")
##kgXref table path
kgXrefPath="/home/zxh/Projects/jkCheng_shanghaiJiaoda/genomicSeq/kgXref_hg19.txt"
##read in kgXref table
#kr <- read.table(kgXrefPath, header=T, sep="\t", comment.char="")
kr=read.delim(kgXrefPath)

##cycling through files
##subDirs under degDir
sds=dir(degDir)
for(sd in sds){
	##cycle through file types
	for(f in ft){
		##full file path
		ff=file.path(degDir, sd, f)
		message(ff, " :: ", date())
		##read in deg table
    deg=read.delim(ff)
    ##attach gene info to deg
    dk=merge(deg, kr, by.x="gene_id", by.y="mRNA", all.x=T)
    ##write out
    of <- file.path(degDir, sd, f)
    of <- paste(of,".geneInfo", sep="")
    write.table(dk, of, row.names=F, col.names=T, sep="\t", quote=F)
	}

}



