#!/usr/bin/env Rscript

options(warn=-1)
suppressPackageStartupMessages({
        require(stringr)
})

# Parsing arguments
args = commandArgs(trailingOnly=TRUE)

# reading command arguments
input = args[1]
output = args[2]
region = args[3]

# defining functions
GetExonID = function(x){
    full_name = x[5]
    first_split = str_split(full_name, pattern = "=", simplify = T)
    if(ncol(first_split) == 1){ 
    	geneID = first_split[,1] 
    }
    else{
    	first_split = first_split[,2]
        geneID = str_split(first_split,pattern = ";", simplify = T)[1,1]
    }
}

GetGeneIDForExons = function(x){
    full_gene_info = x[5]
    first_split = str_split(full_gene_info, pattern = "gene=", simplify = T)
    first_split = first_split[,2]
    geneID = str_split(first_split,pattern = ";", simplify = T)[1,1]
    paste0("gene-", geneID)
}

GetUCENames = function(x){
	UCE = x[1]
	paste0(str_split(UCE, pattern="_", simplify=T)[,1], ".nexus")
}

GetGeneIDForIntrons = function(x){
    full_gene_info = x[5]
    first_split = str_split(full_gene_info, pattern = "ID=", simplify = T)
    # make another split by ";" to get rid of additional information
    # and get only the gene ID
    first_split = first_split[,2]
    str_split(first_split,pattern = ";", simplify = T)[1,1]
}

# read dataframe
df = read.table(input, header=F, sep="\t")

if(region == "exon"){
	df$UCE = apply(df, MARGIN = 1, FUN = GetUCENames)
	df$Exon = apply(df, MARGIN = 1, FUN = GetExonID)
	df$Gene = apply(df, MARGIN = 1, FUN = GetGeneIDForExons)
	df = df[,c(6,7,8)]
}
if(region == "intron"){
	df$UCE = apply(df, MARGIN = 1, FUN = GetUCENames)
	df$Intron = apply(df, MARGIN=1, FUN= GetGeneIDForIntrons)
	df = df[,c(6,7)]
}
if(region == "intergenic"){
	df$UCE = apply(df, MARGIN = 1, FUN = GetUCENames)
	df = df$UCE
}
write.table(df, output, row.names = F, quote=F, col.names = F, sep="\t")