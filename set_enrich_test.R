#!/usr/bin/Rscript
#####################################################################################################
#Functional gene-set enrichment test, for non-genetic data only: Fisher's Exact test, plus refinement
#Andrew Pocklington
#V0.2 03/09/2019

#Functional gene-set enrichment test: Fisher's Exact test, plus refinement

#NOTE: Do not use to test for enrichment in genetic association

#Required parameters:

#--ann_path: Path to file containing annotation (eg. GO, MP) data

#--set_path: Path to file containing gene-sets to be tested

#--dest_dir: Path to output directory for results

#Optional parameters:

#--dest_prefix: Prefix for log & output file names [default = none]

#--ann_format: Format of annotation file [default = standard; other = magma. NOTE: this doesn't recognise all files that can be read into MAGMA, just a specific MAGMA-compatible format: 2 tab-delimited columns (no header, 1 line per gene-set) with field 1 = annotation name & field 2 = space-separated list of entrez gene ids]

#--set_names: Names of gene-sets to be tested [default is to test all gene-sets in set_path file, apart from the background set (if present)]

#--set_format: Format of gene set file (see ann_format) [default = ann_format]

#--bkgd_name: Name of background gene-set [default is to compare each gene-set to all genes for which annotation data is available; other = name of background gene-set]

#--bkgd_path: Path to file containing background gene-set. Only used if --bkgd_name specified [default = set_path]

#--bkgd_format: Format of background file (see ann_format) [default = ann_format]

#--alt: Alternative hypothesis [default = two.sided; other = greater or less]

#--p_thr: Bonferroni-corrected p-value threshold for significance, used to identify gene-sets taken forward for refinement [default = 0.05]

#--p_refine_thr: Raw p-value threshold used to determine gene-sets whose enrichment signal is captured by one or more gene-sets of larger-effect [default = 0.01]

#--core: Number of processing cores to use [default = 1]

#--ann_min: Minimum number of genes per annotation: those with less (after restricting to genes in background) are removed [default = 20]

#--collate: Whether output for all gene_sets bis to be collated into single file [default = FALSE]

####################################################################################
#required libraries
library(optparse)
suppressMessages(library(data.table))
library(Matrix)
library(foreach)
suppressMessages(library(doMC))
#######################################################################
# options - default
#######################################################################

#Path to file containing annotation (eg. GO, MP) data
default_ann_path = NULL

#Path to file containing gene-sets to be tested
default_set_path = NULL

#Path to output directory for results
default_dest_dir = NULL

#--dest_prefix: Prefix for log & output file names 
default_dest_prefix = NULL

#Format of annotation file
default_ann_format = 'standard'

#Format of gene-set file
default_set_format = NULL

#Names of gene-sets to be tested
default_set_names = NULL

#Name of background gene-set
default_bkgd_name = NULL

#Path to file containing background gene-set
default_bkgd_path = NULL

#Format of background file
default_bkgd_format = NULL

#Alternative hypothesis
default_alt = 'two.sided'

#Bonferroni-corrected p-value threshold for significance
default_p_thr = 0.05

#Raw p-value threshold used for refinement
default_p_refine_thr = 0.05

#Number of processing cores to use
default_core = 1

#Minimum number of genes per annotation: those with less (after restricting to genes in background) are removed
default_ann_min = 20

#Maximum number of genes per annotation: those with more (after restricting to genes in background) are removed
default_ann_max = 20000

#Whether output for all gene_sets bis to be collated into single file
default_collate = FALSE

###########################################################################
# option settings
###########################################################################
option_list = list(

make_option("--ann_path", action="store", default=default_ann_path, type='character',help="Path to file containing annotation (eg. GO, MP) data [required]"),
make_option("--set_path", action="store", default=default_set_path, type='character',help="Path to file containing gene-sets to be tested [required]"),
make_option("--dest_dir", action="store", default=default_dest_dir, type='character',help="Path to output directory for results [required]"),
make_option("--dest_prefix", action="store", default=default_dest_prefix, type='character',help="Prefix for log & output file names [default = none]"),
make_option("--ann_format", action="store", default=default_ann_format, type='character',help="Format of annotation file [default = standard; other = magma. NOTE: this doesn't recognise all files that can be read into MAGMA, just a specific MAGMA-compatible format: 2 tab-delimited columns (no header, 1 line per gene-set) with field 1 = annotation name & field 2 = space-separated list of entrez gene ids]"),
make_option("--set_format", action="store", default=default_set_format, type='character',help="Format of gene set file (see ann_format) [default = ann_format]"),
make_option("--set_names", action="store", default=default_set_names, type='character', help="Names of gene-sets to be tested [default is to test all gene-sets in set_path file, apart from the background set (if present)]"),
make_option("--bkgd_name", action="store", default=default_bkgd_name, type='character', help="Name of background gene-set [default is to compare each gene-set to all genes for which annotation data is available]"),
make_option("--bkgd_path", action="store", default=default_bkgd_path, type='character',
	help="Path to file containing background gene-set. Only used if --bkgd_name specified. [default = set_path]"),
make_option("--bkgd_format", action="store", default=default_bkgd_format, type='character',help="Format of background file (see ann_format) [default = ann_format]"),
make_option("--alt", action="store", default=default_alt, type='character', help="Alternative hypothesis [default = two.sided; other = greater or less]"),
make_option("--p_thr", action="store", default=default_p_thr, type='double',help="Bonferroni-corrected p-value threshold for significance, used to identify gene-sets taken forward for refinement [default = 0.05]"),
make_option("--p_refine_thr", action="store", default=default_p_refine_thr, type='double',help="Bonferroni-corrected p-value threshold used to determine gene-sets whose enrichment signal is captured by one or more gene-sets of larger-effect [default = 0.05]"),
make_option("--core", action="store", default=default_core, type='integer',help="Number of processing cores to use [default = 1]"),
make_option("--ann_min", action="store", default=default_ann_min, type='integer',help="Minimum number of genes per annotation: those with less (after restricting to genes in background) are removed [default = 20]"),
make_option("--ann_max", action="store", default=default_ann_max, type='integer',help="Maximum number of genes per annotation: those with more (after restricting to genes in background) are removed [default = 20,000]"),
make_option("--collate", action="store", default=default_collate, type='logical',help="Whether output for all gene_sets bis to be collated into single file [default = FALSE]")
)
#######################################################################
# functions
#######################################################################

###########################################################################
#read gene annotations from 'standard' format file: tab-delimited, 2 required fields ('ann','ncbi'), 1 line per gene-set
#                                                   may have additional fields (e.g., 'id','symbol')
#'ann' field = annotation name
#'ncbi' field = pipe ('|') separated list of ncbi/entrez gene ids
#'symbol' field = pipe ('|') separated list of gene symbols (symbol field not used here)
#order of genes in ncbi & symbol field (if present) must be the same
# e.g.
#ann	ncbi	symbol
#ann_name1	entrez1|entrez2|...	symbol1|symbol2|...
#ann_name2...
read_standard <- function(next_path) {
	
	tmp <- fread(next_path,sep = '\t')
	
	tmp_list = lapply(tmp$ann,function(k) unlist(strsplit(tmp[ann == k,ncbi],'|',fixed = TRUE)))
	names(tmp_list) = tmp$ann
	
	return(tmp_list)
}
###########################################################################
#read gene annotations from MAGMA format file: 2 tab-delimited fields (no field names/header), 1 line per gene-set
#field 1 = annotation name
#field 2 = space (' ') separated list of entrez gene ids
# e.g.
#ann_name1	entrez1 entrez2 entrez3 ...
#ann_name2...
read_MAGMA <- function(next_path) {
	
	tmp = fread(next_path,header = FALSE,col.names = c('ann','ncbi'),sep = '\t')
	
	tmp_list = lapply(tmp$ann,function(k) unlist(strsplit(tmp[ann == k,ncbi],' ',fixed = TRUE)))
	names(tmp_list) = tmp$ann
	
	return(tmp_list)
}

#read file
read_file <- function(path,format) {
	
	switch(format,
	
		'standard' = read_standard(path),
		'magma' = read_MAGMA(path)
	)
	
}

######################
#read annotation sets
read_ann_data <- function(path,format) {
	
	cat('Reading annotations from file:',path,'\n')
	
	start.time = Sys.time()
	data = read_file(path,format)
	time.taken = Sys.time() - start.time
	
	cat(length(data),'annotations read, time taken =',time.taken,attr(time.taken, 'units'),'\n\n')
	
	data
}

#read gene-sets for analysis
read_set_data <- function(path,format) {
	
	cat('Reading gene-sets from file:',path,'\n')
	
	start.time = Sys.time()
	if (is.null(format)) {
		
		data = read_file(path,opt$ann_format)
				
	} else {
				
		data = read_file(path,format)
	}	
	time.taken = Sys.time() - start.time
	
	cat(length(data),'gene-sets read, time taken =',time.taken,attr(time.taken, 'units'),'\n\n')
	
	data
}

############################
#read/extract background set
get_background <- function(path,name,format,ann_data,set_data) {
	
	#if path not specified
	if(is.null(path)) {
	
	#if name not specified
	if(is.null(name)) {
		
		#set comparator to all genes for which annotation data is available
		bkgd_genes = as.character(unique(unlist(ann_data)))
		
		cat('All genes for which annotation data is available used to define background set\n')
	}
	else {#otherwise, extract & remove comparator set from set_data
		
		bkgd_genes = as.character(unique(set_data[[name]]))
		set_data <<- set_data[setdiff(names(set_data),c(name))]
		
		cat('Background set:',name,'\nExtracted from set file\n')
	}	
	} else {#if path specified
	
		#if name not specified, raise an error
		if(is.null(name)) {stop('--bkgd_name must be specified if --bkgd_path specified')
		}
		else{#otherwise, read file and extract comparator set
			
			if (is.null(format)) {
				
				bkgd_genes = as.character(unique(read_ann_data(path,opt$ann_format)[[name]]))
				
			} else {
				
				bkgd_genes = as.character(unique(read_ann_data(path,format)[[name]]))
			}
			
			cat('Background set:', name,'\nRead from file:',path,'\n')
		}
	
	}
	cat(length(bkgd_genes),'genes in background\n\n')
	
	bkgd_genes
}

######################################
#remove any gene-sets not for analysis
process_set_data <- function(set_data,opt_names) {
	
	#if --set_names not specified by user, test all sets in set_data
	if (is.null(opt_names)) {
		
		set_names = names(set_data)
		
		cat('Testing all gene-sets:\n')
		
	} else {#if --set_names specified by user, only test these sets & remove all others from set_data
		
		set_names = unlist(strsplit(opt_names,',',fixed = TRUE))
		set_data <<- set_data[set_names]
		
		cat('Testing specified gene-sets:\n')
	}
	for(next_name in set_names) {cat(length(set_data[[next_name]]),'\t',next_name,'\n')}
	cat('\n')
	
	set_names
}

###################################################################################################################
#create binary matrix recording background genes present in each annotation & gene-set (row = gene, col = gene-set)
#removing any annotations with < ann_min or > ann_max genes after restricting to background: also remove any gene sets
#containing all background genes
make_gene_ann_mtx <- function(bkgd_genes,ann_data,set_data,set_names,ann_min,ann_max) {
	
	cat('Constructing gene-annotation matrix for',length(bkgd_genes),'background genes...\n\n')
	start.time = Sys.time()
	
	#collate all annotation names
	ann_names = names(ann_data)
	all_names = c(ann_names,set_names)
	
	#number of genes in background set
	bkgd_n = length(bkgd_genes) 
	
	gene_ann_mtx  <- Matrix(0, nrow = bkgd_n, ncol = length(all_names), dimnames = list(bkgd_genes,all_names),sparse = TRUE)
	
	for (next_ann in ann_names) {
		
		next_genes = intersect(bkgd_genes,ann_data[[next_ann]])
		if (length(next_genes) > 0) { gene_ann_mtx[next_genes,next_ann] = 1 }
	}
	
	for (next_ann in set_names) {
		
		next_genes = intersect(bkgd_genes,set_data[[next_ann]])
		if (length(next_genes) > 0) { gene_ann_mtx[next_genes,next_ann] = 1 }
	}
	
	#identify annotations with >= ann_min genes and <= ann_max (or number of background genes-1, whichever is smaller)
	ann_n = colSums(gene_ann_mtx)
	effective_max = min(ann_max,bkgd_n - 1)
	test_ann = intersect(ann_names,dimnames(gene_ann_mtx)[[2]][which(ann_n >= ann_min & ann_n <= effective_max)])
	
	#restrict annotation data to these
	ann_data <<- ann_data[test_ann]
	
	#write summary to log
	time.taken = Sys.time() - start.time
	cat(length(test_ann),'annotations with >=',ann_min,'genes and <=', effective_max,'genes, time taken =',time.taken,attr(time.taken, 'units'),'\n\n')
	
	gene_ann_mtx
}

###################################
#perform annotation enrichment test
test_set <- function(set_name,ann_names,gene_ann_mtx,alt) {
	
	cat(set_name,'\nPerforming enrichment test...')
	start.time = Sys.time()
	
	#binary vector denoting which background genes are in gene-set
	set_col = gene_ann_mtx[,set_name]
	
	results <- foreach(ann = ann_names,.combine=rbind) %dopar% {
	
		fisher_test(set_col,ann,gene_ann_mtx[,ann],alt)
	
	} #foreach(i = 1:GO_test_N, .combine=cbind) %dopar% 
	
	results$P_bonf = p.adjust(results$P,method = 'bonferroni')
	
	#write summary to log
	time.taken = Sys.time() - start.time
	cat('time taken =',time.taken,attr(time.taken, 'units'),'\n\n')
	
	#order results & return
	order_indx = order(results$P)
	
	results[order_indx,]
}

###################################
fisher_test <- function(set_col,ann,ann_col,alt) {
	
	t = fisher.test(set_col,ann_col,alternative = alt)
	
	data.table(ann = ann,ann_N = sum(ann_col),overlap_N = sum(set_col*ann_col),OR = t$estimate,P = t$p.value)

}

#############
#save results
save_test_results <- function(dest_path_root,set_name,results) {
	
	path = paste(dest_path_root,set_name,'_results.txt',sep='')
	
	write.table(results,file = path,row.names = FALSE,col.names = TRUE,sep = '\t',quote = FALSE)
	
	cat('Test results written to file: ',path,'\n\n')
}
###########################################################################
#refine results, prioritising gene-sets by odds ratio (largest to smallest)
refine_results <- function(test_results,alt,p_thr,p_refine_thr,next_set,gene_ann_mtx) {

	#identify gene-sets with Bonferroni-corrected P < p_thr
	sig_ann = subset(test_results,P_bonf < p_thr)
	
	cat(dim(sig_ann)[1],'gene-sets with Bonferroni-corrected P <',p_thr,'\n\n')
	cat('Refining enrichment...')
	
	#if there are no significant annotations to refine, return empty table
	if (dim(sig_ann)[1] == 0) {
		
		sig_ann$refined_ann = sig_ann$ann
		
		return(sig_ann)
	}
	
	#names of refined set of annotations
	refined = NULL
	
	#binary vector denoting which background genes are in gene-set
	set_col = gene_ann_mtx[,next_set]
	
	#sort annotations by odds ratio
	order_indx = order(sig_ann$OR,decreasing = TRUE)
	ann_table = data.table(ann = sig_ann$ann[order_indx])
	
	#binary vector denoting which background genes are to be removed
	remove_col = rep(0,dim(gene_ann_mtx)[1])
	
	#while there is more than one annotation left...
	while (dim(ann_table)[1] > 1) {
		
		#take annotation with largest odds ratio...
		next_ann = ann_table$ann[1]
		
		#...add it's genes to those to be removed
		next_ann_col = gene_ann_mtx[ ,next_ann]
		
		remove_col = remove_col + next_ann_col - (remove_col*next_ann_col)
		
		#...remove refined set genes from annotations and re-test their enrichment
		ann_table = cbind(ann_table,remove_col_test(set_col,remove_col,ann_table$ann,alt,gene_ann_mtx))
		
		#calculate Bonferroni-corrected p-value threshold
		bonf_thr = p_refine_thr/dim(ann_table)[1]
		
		#...add annotation to list of refined sets, storing data for it & other annotations
		#whose remaining enrichment it captures (ie. those with P >= p_refine_thr)
		next_non_sig = subset(ann_table,P >= bonf_thr)$ann
		
		next_summary = cbind(data.table('refined_ann' = rep.int(next_ann,length(next_non_sig))),subset(sig_ann,ann %in% next_non_sig))
		
		refined = rbind(refined,next_summary)
		
		#extract all remaining annotations with P < p_refine_thr & sort by residual odds ratio
		tmp = subset(ann_table,P < bonf_thr)
		order_indx = order(tmp$OR,decreasing = TRUE)
		ann_table = data.table(ann = tmp$ann[order_indx])
		
	}#while (dim(ann_table)[1] > 1)
	
	#if there is one annotation remaining, add it to refined list of annotations
	if (dim(ann_table)[1] == 1) {
		
		next_ann = ann_table$ann[[1]]
		
		next_summary = cbind(data.table('refined_ann' = next_ann),subset(sig_ann,ann == next_ann))
		
		refined = rbind(refined,next_summary)
	}
	
	#return refined results
	refined
}

#remove refined genes from annotations and re-test their enrichment
remove_col_test <- function(set_col,remove_col,ann,alt,gene_ann_mtx) {
	
	test_results <- foreach(k = ann, .combine=rbind) %dopar% {
		
		#binary vector denoting which background genes are in annotation
		next_col = gene_ann_mtx[,k]
		
		#remove refined genes from annotation
		residual_col = next_col - (next_col*remove_col)
			
		#perform Fisher's Exact test on residual annotation
		if (max(residual_col) == 0) {
			
			next_results  = data.table(P = 1,OR = -9999)
				
		} else {
				
			t = fisher.test(set_col,residual_col,alternative = alt)
	
			next_results = data.table(P = t$p.value,OR = t$estimate)

		}
	
		next_results
	
	} #foreach(k = ann, .combine=cbind) %dopar%
	
	test_results
}

#####################
#save refined results
save_refined <- function(dest_path_root,set_name,refined) {
	
	path = paste(dest_path_root,set_name,'_refined.txt',sep='')
	
	write.table(refined,file = path,row.names = FALSE,col.names = TRUE,sep = '\t',quote = FALSE)
	
	cat('Refined results written to file: ',path,'\n\n')
}
####################################################################################
#save refined groups (annotations whose enrichment is captured by each refined term)
save_refined_groups <- function(dest_path_root,set_name,refined) {
	
	path = paste(dest_path_root,set_name,'_refined_groups.txt',sep='')
	
	write.table(refined,file = path,row.names = FALSE,col.names = TRUE,sep = '\t',quote = FALSE)
	
	cat('Refined groups written to file: ',path,'\n\n')
}
###########################################################################
#initialisation
###########################################################################
#start time
start.time = Sys.time()

#variables
opt = parse_args(OptionParser(option_list=option_list))

#set up parallel processing
registerDoMC(opt$core)

#set up root path for output files
if (is.null(opt$dest_prefix)) {
	
	dest_path_root = opt$dest_dir
} else {
	
	dest_path_root = paste(opt$dest_dir,opt$dest_prefix,'_',sep='')
}


#log file
log_path = paste(dest_path_root,'set_enrich_test.log',sep='')
sink(file = log_path, append = F)
cat('#######################################################################################################
##Functional gene-set enrichment test, for non-genetic data only: Fisher\'s Exact test, plus refinement
#Andrew Pocklington
#V0.11 29/07/2019
#######################################################################################################\n')
print(opt)
cat('Started at',as.character(start.time),'\n\n')
#########################################################################################################################################
# main - start
#########################################################################################################################################

#read annotation sets
ann_data = read_ann_data(opt$ann_path,opt$ann_format)

#read gene-sets for analysis (file may contain extra gene-sets)
set_data = read_set_data(opt$set_path,opt$set_format)

#read/extract background set
bkgd_genes = get_background(opt$bkgd_path,opt$bkgd_name,opt$bkgd_format,ann_data,set_data)

#remove any gene-sets not for analysis
set_names = process_set_data(set_data,opt$set_names)

#create binary matrix recording background genes present in each annotation & gene-set (row = gene, col = gene-set)
#removing any annotations with < ann_min or > ann_max genes after restricting to background
gene_ann_mtx = make_gene_ann_mtx(bkgd_genes,ann_data,set_data,set_names,opt$ann_min,opt$ann_max)

#names of annotations to be tested
ann_names = names(ann_data)

#collated output (if necessary)
collated_results = NULL
collated_refined = NULL
collated_refined_groups = NULL

#for each gene-set
for (next_set in set_names) {
	
	#perform annotation enrichment test
	test_results = test_set(next_set,ann_names,gene_ann_mtx,opt$alt)

	#save raw results
	if (opt$collate) {
		
		tmp_N = dim(test_results)[1]
		set_num = sum(gene_ann_mtx[,next_set])
		
		tmp = data.table(set = rep(next_set, tmp_N),set_N = rep(set_num,tmp_N))
		
		collated_results = rbind(collated_results,cbind(tmp,test_results))
		
		
	} else {
		
		save_test_results(dest_path_root,next_set,test_results)
	}
	
	#refine results
	refined_groups = refine_results(test_results,opt$alt,opt$p_thr,opt$p_refine_thr,next_set,gene_ann_mtx)
		
	refined = subset(refined_groups,refined_ann == ann,select = c('ann','ann_N','overlap_N','OR','P','P_bonf'))

	#write summary to log
	cat(dim(refined)[1],'refined annotations\n\n')
	
	#save refined results
	if (opt$collate) {
		
		set_num = sum(gene_ann_mtx[,next_set])
		
		tmp_N = dim(refined)[1]
		tmp = data.table(set = rep(next_set,tmp_N),set_N = rep(set_num,tmp_N))
		collated_refined = rbind(collated_refined,cbind(tmp,refined))
		
		tmp_N = dim(refined_groups)[1]
		tmp = data.table(set = rep(next_set,tmp_N),set_N = rep(set_num,tmp_N))
		collated_refined_groups = rbind(collated_refined_groups,cbind(tmp,refined_groups))
		
	} else {
		
		save_refined(dest_path_root,next_set,refined)
		save_refined_groups(dest_path_root,next_set,refined_groups)
	}

	#save collated data if necessary
	if (opt$collate) {
		
		save_test_results(dest_path_root,'collated',collated_results)
		save_refined(dest_path_root,'collated',collated_refined)
		save_refined_groups(dest_path_root,'collated',collated_refined_groups)
	}
	
}
#########################################################################################################################################
# main - end
#########################################################################################################################################
time.taken = Sys.time() - start.time
cat('Total time taken =',time.taken,attr(time.taken, 'units'),'\n\n')

##############################
#stop sending data to log file
sink()