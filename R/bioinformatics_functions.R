# Eskew et al. wood frog/American bullfrog transcriptomics study

# Bioinformatics Functions

#==============================================================================


# Custom functions used for bioinformatics analyses


# Lookup annotation data by contig name

anno_lookup_contig <- function(x, BP = FALSE) {
  
  x <- paste("Lithobates.clamitans_contig_", as.character(x), sep = "")
  
  if (BP == TRUE) {
    print(data$genes[data$genes$Contig_Name == x, ])
    print(anno.table[anno.table$contig_identifier == x & 
                       anno.table$namespace_1003 == "biological_process", ])
  }
  
  else {
    print(data$genes[data$genes$Contig_Name == x, ])
    print(anno.table[anno.table$contig_identifier == x, ])
  }
}


# Lookup annotation data by contig name in the Trinotate database

anno_lookup_contig_trinotate <- function(x, bd = FALSE) {
  
  if (bd == TRUE) {
    x <- paste("BDET_", as.character(x), sep = "")
    print(x)
    
    print(data$genes[data$genes$Contig_Name == x, ])
    print(trinotate.table[trinotate.table$"#gene_id" == x, c(3, 7, 13)])
  }
  
  else {
    x <- paste("Lithobates.clamitans_contig_", as.character(x), sep = "")
    
    print(data$genes[data$genes$Contig_Name == x, ])
    print(trinotate.table[trinotate.table$"#gene_id" == x, c(3, 7, 13)])
  }
}


# Lookup annotation data by Entrez ID
# (There may be cases where a contig doesn't appear in the anno.table
# despite the Entrez ID being present if it's also assigned to another contig)

anno_lookup_entrez <- function(x, BP = FALSE) {
  
  if (BP == TRUE) {
    print(data$genes[match(x, data$genes$Entrez_ID)[1], ])
    print(anno.table[anno.table$entrezgene == x &
                       anno.table$namespace_1003 == "biological_process", ])
  }
  
  else {
    print(data$genes[match(x, data$genes$Entrez_ID)[1], ])
    print(anno.table[anno.table$entrezgene == x, ])
  }
}


# Lookup annotation data by GO ID

anno_lookup_go <- function(x, head = F) {
  
  x <- paste("GO:", x, sep = "")
  if(head == T) print(head(anno.table[anno.table$go_id == x, ]))
  else print(anno.table[anno.table$go_id == x, ])
}


# Lookup all annotation data in a list of genes

anno_lookup_list <- function(x) {
  
  for(gene in x) {
    print(anno.table[anno.table$contig_identifier == gene, ])
  }
}


# Calculate number of differentially expressed genes in the set

get_num_de <- function(x) {
  
  length(which(decideTestsDGE(x) != 0))
}


# Get differentially expressed genes

get_de <- function(x) {
  
  rownames(data$counts)[decideTestsDGE(x) != 0]
}


# Get upregulated genes

get_de_up <- function(x) {
  
  rownames(data$counts)[decideTestsDGE(x) == 1]
}


# Get downregulated genes

get_de_down <- function(x) {
  
  rownames(data$counts)[decideTestsDGE(x) == -1]
}


# Generate GOstats parameters for a vector of contigs

get.go.params <- function(x) {
  
  GSEAGOHyperGParams(name = "My custom annotation",
                     geneSetCollection = gsc, geneIds = x, 
                     universeGeneIds = universe, ontology = "BP", 
                     pvalueCutoff = 0.05, conditional = TRUE, 
                     testDirection = "over")
}
