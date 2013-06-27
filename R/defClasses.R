##################################################
# class def biomvRCNS using GRanges
##################################################	
setClass("biomvRCNS",
	representation(
		x = "GRanges",
		res = "GRanges",
		param = "list"
	),
	prototype = list(
		x = GRanges(),
		res = GRanges(),
		param = list()
	)
)
## Show method
setMethod("show", "biomvRCNS", 
	function(object) {
		res = c(sprintf("Object is of class: 'biomvRCNS'\n"),
					sprintf("List of parameters used in the model:\n"))
		cat(res, sep="")
		cat(names(object@param), fill=T, sep=', ')
		cat("\nThe segmented ranges:\n")
		show(object@res)
	})
	
## Plot method
setMethod("plot", "biomvRCNS",  
	function(x, sampleInOne=TRUE,...) {
		for(seq in as.character(unique(seqnames(x@x)))){
			if(sampleInOne){
				biomvRGviz(exprgr=x@x[seqnames(x@x)==seq, unique(mcols(x@res)[,'SAMPLE'])], seggr=x@res,...) 
			} else {
				for(s in unique(values(x@res)[,'SAMPLE'])){
					biomvRGviz(exprgr=x@x[seqnames(x@x)==seq, s], seggr=x@res[values(x@res)[,'SAMPLE']==s],...) 
				}
			}
		}
	})	
	
	

