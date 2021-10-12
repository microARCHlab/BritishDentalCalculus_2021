
# British Genus and Species Levels Co-Occurrence Plots
#
# Gancz 5/27/2021
#
# Usage: Analyze CCLasso objects from British datasets and generate visualizations

############################################################################################################################

# Package Installation

############################################################################################################################

set.seed(5272021)

#setting up libraries
library(testit)
#install.packages("igraph")
library(igraph)
#install.packages("dplyr")
library(dplyr)
#install.packages("RColorBrewer")
library(RColorBrewer)
#install.packages("flashClust")
library(flashClust)
#install.packages("dynamicTreeCut")
library(dynamicTreeCut)
#install.packages("reshape2")
library(reshape2)
#install.packages("parallel")
library(parallel)
#install.packages("svglite")
library(svglite)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("cowplot")
library(cowplot)

############################################################################################################################

# Import Data

############################################################################################################################

#Genus Level CCLasso object: Genus_res_ccl_count
Genus_pvals  <- data.frame(Genus_res_ccl_count$p_vals)
  colnames(Genus_pvals) <- names(Genus)
  row.names(Genus_pvals) <- names(Genus)
Genus_cors <- Genus_res_ccl_count$cor_w
  colnames(Genus_cors) <- names(Genus)
  row.names(Genus_cors) <- names(Genus)
  
#Genus Strep Level CCLasso object: Genus_strep_res_ccl_count
Genus_strep_pvals  <- data.frame(Genus_strep_res_ccl_count$p_vals)
  colnames(Genus_strep_pvals) <- names(Genus_strep)
  row.names(Genus_strep_pvals) <- names(Genus_strep)
Genus_strep_cors <- Genus_strep_res_ccl_count$cor_w
  colnames(Genus_strep_cors) <- names(Genus_strep)
  row.names(Genus_strep_cors) <- names(Genus_strep)
  
#Genus meth Level CCLasso object: Genus_strep_res_ccl_count
  Genus_meth_pvals  <- data.frame(Genus_meth_res_ccl_count$p_vals)
    colnames(Genus_meth_pvals) <- names(Genus_meth)
    row.names(Genus_meth_pvals) <- names(Genus_meth)
  Genus_meth_cors <- Genus_meth_res_ccl_count$cor_w
    colnames(Genus_meth_cors) <- names(Genus_meth)
    row.names(Genus_meth_cors) <- names(Genus_meth)
    
#Genus Actino Level CCLasso object: Genus_strep_res_ccl_count
  Genus_Actino_pvals  <- data.frame(Genus_Actino_res_ccl_count$p_vals)
    colnames(Genus_Actino_pvals) <- names(Genus_Actino)
    row.names(Genus_Actino_pvals) <- names(Genus_Actino)
  Genus_Actino_cors <- Genus_Actino_res_ccl_count$cor_w
    colnames(Genus_Actino_cors) <- names(Genus_Actino)
    row.names(Genus_Actino_cors) <- names(Genus_Actino)
  
#Species Level CCLasso objet: Species_res_ccl_count
Species_pvals  <- Species_res_ccl_count$p_vals
  colnames(Species_pvals) <- names(Species)
  row.names(Species_pvals) <- names(Species)
Species_cors <- Species_res_ccl_count$cor_w
  colnames(Species_cors) <- names(Species)
  row.names(Species_cors) <- names(Species)
  
  
#Genus Level CCLasso object: genus_samp70_res_ccl_count
genus_samp70_pvals  <- data.frame(genus_samp70_res_ccl_count$p_vals)
  colnames(genus_samp70_pvals) <- names(genus_samp70_trans)
  row.names(genus_samp70_pvals) <- names(genus_samp70_trans)
genus_samp70_cors <- genus_samp70_res_ccl_count$cor_w
  colnames(genus_samp70_cors) <- names(genus_samp70_trans)
  row.names(genus_samp70_cors) <- names(genus_samp70_trans)
  
  
#General parameters
cor_threshold <- 0.3
pval_threshold <- 0.1


############################################################################################################################

# Network Metrics

############################################################################################################################

#Setting up functions for Network Analysis 
# Function to calculate network stats:
statme <- function(x) {
  # Centrality
  V(x)$degree      <- degree(x, mode = "total")
  V(x)$betweenness <- betweenness(x)
  V(x)$evcent      <- evcent(x)$vector
  V(x)$closeness   <- closeness(x)
  E(x)$ebetweenness <- edge.betweenness(x)
  V(x)$transitivity <- transitivity(x, type = "local")
  
  # Local position
  V(x)$effsize     <- effective.size(x, mode = "all")
  V(x)$constraint  <- constraint(x)
  
  # Clustering
  com <- edge.betweenness.community(x)
  V(x)$memb        <- com$membership
  
  # Whole network
  set.graph.attribute(x, "density", graph.density(x))
  set.graph.attribute(x, "avgpathlength", average.path.length(x))
  set.graph.attribute(x, "modularity", modularity(com))
  set.graph.attribute(x, "betcentralization", centralization.betweenness(x)$centralization)
  set.graph.attribute(x, "degcentralization", centralization.degree(x, mode = "total")$centralization)
  set.graph.attribute(x, "size", vcount(x))
  set.graph.attribute(x, "edgecount", ecount(x))
  
  return(x)
}

#Genus Network Stats
        #reformat tables
          Genus_pvals_reformat <- Genus_pvals[, 1:ncol(Genus_pvals)]
          Genus_cors_pruned <- as.matrix(Genus_cors)
        #prune insig edges  
          Genus_cors_pruned[which(Genus_pvals_reformat > 0.0002)] <- 0 #check what the pvalue thereshold should be based on CCLASSO lit
        #prune out absolute value of correlations less than a certain value
          Genus_cors_pruned[abs(Genus_cors_pruned) < cor_threshold] <- 0
        #turn into igraph object
          g_Genus <- graph.adjacency(Genus_cors_pruned, mode = "undirected", weighted = TRUE)
          g_Genus <- simplify(g_Genus)
          plot(g_Genus)
          
#Genus strep Network Stats
        #reformat tables
          Genus_strep_pvals_reformat <- Genus_strep_pvals[, 1:ncol(Genus_strep_pvals)]
          Genus_strep_cors_pruned <- as.matrix(Genus_strep_cors)
        #prune insig edges  
          Genus_strep_cors_pruned[which(Genus_strep_pvals_reformat > 0.0002)] <- 0 #check what the pvalue thereshold should be based on CCLASSO lit
        #prune out absolute value of correlations less than a certain value
          Genus_strep_cors_pruned[abs(Genus_strep_cors_pruned) < cor_threshold] <- 0
        #turn into igraph object
          g_Genus_strep <- graph.adjacency(Genus_strep_cors_pruned, mode = "undirected", weighted = TRUE)
          g_Genus_strep <- simplify(g_Genus_strep)
          plot(g_Genus_strep)
          
#Genus meth Network Stats
        #reformat tables
          Genus_meth_pvals_reformat <- Genus_meth_pvals[, 1:ncol(Genus_meth_pvals)]
          Genus_meth_cors_pruned <- as.matrix(Genus_meth_cors)
        #prune insig edges  
          Genus_meth_cors_pruned[which(Genus_meth_pvals_reformat > 0.0002)] <- 0 #check what the pvalue thereshold should be based on CCLASSO lit
        #prune out absolute value of correlations less than a certain value
          Genus_meth_cors_pruned[abs(Genus_meth_cors_pruned) < cor_threshold] <- 0
        #turn into igraph object
          g_Genus_meth <- graph.adjacency(Genus_meth_cors_pruned, mode = "undirected", weighted = TRUE)
          g_Genus_meth <- simplify(g_Genus_meth)
          plot(g_Genus_meth)
          
          
#Genus actino Network Stats
        #reformat tables
          Genus_Actino_pvals_reformat <- Genus_Actino_pvals[, 1:ncol(Genus_Actino_pvals)]
          Genus_Actino_cors_pruned <- as.matrix(Genus_Actino_cors)
        #prune insig edges  
          Genus_Actino_cors_pruned[which(Genus_Actino_pvals_reformat > 0.0002)] <- 0 #check what the pvalue thereshold should be based on CCLASSO lit
        #prune out absolute value of correlations less than a certain value
          Genus_Actino_cors_pruned[abs(Genus_Actino_cors_pruned) < cor_threshold] <- 0
        #turn into igraph object
          g_Genus_Actino <- graph.adjacency(Genus_Actino_cors_pruned, mode = "undirected", weighted = TRUE)
          g_Genus_Actino <- simplify(g_Genus_Actino)
          plot(g_Genus_Actino)

#Species Network Stats
          #reformat tables
          Species_pvals_reformat <- Species_pvals[, 1:ncol(Species_pvals)]
          Species_cors_pruned <- as.matrix(Species_cors)
          #prune insig edges  
          Species_cors_pruned[which(Species_pvals_reformat > 0.0002)] <- 0 #check what the pvalue thereshold should be based on CCLASSO lit
          #prune out absolute value of correlations less than a certain value
          Species_cors_pruned[abs(Species_cors_pruned) < cor_threshold] <- 0
          #turn into igraph object
          g_Species <- graph.adjacency(Species_cors_pruned, mode = "undirected", weighted = TRUE)
          g_Species <- simplify(g_Species)
          plot(g_Species)
          
#Genus Network Stats Sammplee 70, min freq 50000
          #reformat tables
          genus_samp70_pvals_reformat <- genus_samp70_pvals[, 1:ncol(genus_samp70_pvals)]
          genus_samp70_cors_pruned <- as.matrix(genus_samp70_cors)
          #prune insig edges  
          genus_samp70_cors_pruned[which(genus_samp70_pvals_reformat > 0.0002)] <- 0 #check what the pvalue thereshold should be based on CCLASSO lit
          #prune out absolute value of correlations less than a certain value
          genus_samp70_cors_pruned[abs(genus_samp70_cors_pruned) < cor_threshold] <- 0
          #turn into igraph object
          g_genus_samp70 <- graph.adjacency(genus_samp70_cors_pruned, mode = "undirected", weighted = TRUE)
          g_genus_samp70s <- simplify(g_genus_samp70)
          plot(g_genus_samp70)
          
          
############################################################################################################################

# Visualize Correlation Matrix: Genus

############################################################################################################################

#Genus
    #Graph of all edges transformed into positive 
    Genus_cors_pos <- abs(Genus_cors_pruned)
    g_Genus_pos <- graph.adjacency(Genus_cors_pos, mode = "undirected", weighted= TRUE)
    g_Genus_pos <- simplify(g_Genus_pos)     
    plot(g_Genus_pos)

    #Graph of only positive co-occurrences
    Genus_pos_cors <- Genus_cors_pruned
    Genus_pos_cors[Genus_pos_cors< 0] <- 0
    Genus_pos_graph <- graph.adjacency(Genus_pos_cors, mode ='undirected', weighted = TRUE)
    Genus_pos_graph <- simplify(Genus_pos_graph)
    plot(Genus_pos_graph)
    plot(Genus_pos_graph, rescale = TRUE, edge.color = ifelse(E(Genus_pos_graph)$weight> 0, "pink", "blue"), vertex.label.font = 3, vertex.size = deg*.5, vertex.label.cex = 1)
    
    deg <- degree(Genus_pos_graph, mode = "all")
    
    #export pos cooccur plot    
    svglite(filename = "all_pos_cooccur_clusters.svg", width =10, height =10)
    plot(Genus_pos_graph, rescale = TRUE, edge.color = ifelse(E(Genus_pos_graph)$weight> 0, "pink", "blue"), vertex.label.font = 3, vertex.size = deg*.5, vertex.label.cex = 1)
    dev.off()   

    #plot of only negative co-occurrences
    Genus_neg_cors <- Genus_cors_pruned
    Genus_neg_cors[Genus_neg_cors > 0] <- 0
    Genus_neg_graph <- graph.adjacency(Genus_neg_cors, mode ='undirected', weighted = TRUE)
    Genus_pos_graph <- simplify(Genus_neg_graph)
    plot(Genus_neg_graph)

    #Visualizations
    g <- plot(g_Genus, rescale = TRUE, edge.width = abs(E(g_Genus)$weight)*1, edge.color = ifelse(E(g_Genus)$weight >0, "lightpink", "lightblue"), vertex.label.font = 3)
    
    A <- layout_with_fr(g_Genus)
    plot(g_Genus, rescale = TRUE, edge.width = abs(E(g_Genus)$weight)*.5, edge.color = ifelse(E(g_Genus)$weight >0, "lightpink", "lightblue"), vertex.label.font = 3, layout = A, vertex.size = deg*.1, vertex.label.cex = deg*.01)
    
    deg <- degree(g_Genus, mode = "all")
    plot(g_Genus, rescale = TRUE, edge.width = abs(E(g_Genus)$weight)*.5, edge.color = ifelse(E(g_Genus)$weight >0, "pink", "blue"), vertex.label.font = 3, vertex.label.cex = deg*.01, vertex.size = deg*.5, main ="Authorities")
    
    A = layout_with_fr(g_Genus)
    plot(g_Genus, rescale = TRUE, edge.color = ifelse(E(g_Genus)$weight >0, "pink", "lightblue"), vertex.label.font = 3, layout = A, vertex.size = deg*.1, vertex.label.cex= 1)
    
    #export cooccur plot    
    svglite(filename = "all_cooccur.svg", width =10, height =10)
    plot(g_Genus, rescale = TRUE, edge.color = ifelse(E(g_Genus)$weight >0, "pink", "lightblue"), vertex.label.font = 3, layout = A, vertex.size = deg*.1, vertex.label.cex= 1)
    dev.off()  
    
############################################################################################################################
    
# Visualize Correlation Matrix: Genus: Strep
    
############################################################################################################################       
    
    deg <- degree(g_Genus_strep, mode = "all")
    
    #Graph of all edges transformed into positive 
    Genus_strep_cors_pos <- abs(Genus_strep_cors_pruned)
    g_Genus_strep_pos <- graph.adjacency(Genus_strep_cors_pos, mode = "undirected", weighted= TRUE)
    g_Genus_strep_pos <- simplify(g_Genus_strep_pos)     
    plot(g_Genus_strep_pos)
    
    #Graph of only positive co-occurrences
    Genus_strep_pos_cors <- Genus_strep_cors_pruned
    Genus_strep_pos_cors[Genus_strep_pos_cors< 0] <- 0
    Genus_strep_pos_graph <- graph.adjacency(Genus_strep_pos_cors, mode ='undirected', weighted = TRUE)
    Genus_strep_pos_graph <- simplify(Genus_strep_pos_graph)
    plot(Genus_strep_pos_graph)
    plot(Genus_strep_pos_graph, rescale = TRUE, edge.color = ifelse(E(Genus_strep_pos_graph)$weight> 0, "pink", "blue"), vertex.label.font = 3, vertex.size = deg*.2, vertex.label.cex = 1)
    
    #export pos cooccur plot    
    svglite(filename = "strep_pos_cooccur.svg", width =10, height =10)
    plot(Genus_strep_pos_graph, rescale = TRUE, edge.color = ifelse(E(Genus_strep_pos_graph)$weight> 0, "pink", "blue"), vertex.label.font = 3, vertex.size = deg*.2, vertex.label.cex = 1)
    dev.off()   

    A <- layout_with_fr(g_Genus_strep)
    plot(g_Genus_strep, rescale = TRUE, edge.color = ifelse(E(g_Genus_strep)$weight >0, "lightpink", "lightblue"), vertex.label.font = 3, layout = A, vertex.size = deg*.1, vertex.label.cex = deg*.01)
    
    
############################################################################################################################
    
# Visualize Correlation Matrix: Genus: Meth
    
############################################################################################################################       
  
    deg <- degree(g_Genus_meth, mode = "all")      
    
  #Graph of all edges transformed into positive 
    Genus_meth_cors_pos <- abs(Genus_meth_cors_pruned)
    g_Genus_meth_pos <- graph.adjacency(Genus_meth_cors_pos, mode = "undirected", weighted= TRUE)
    g_Genus_meth_pos <- simplify(g_Genus_meth_pos)     
    plot(g_Genus_meth_pos)
    
  #Graph of only positive co-occurrences
    Genus_meth_pos_cors <- Genus_meth_cors_pruned
    Genus_meth_pos_cors[Genus_meth_pos_cors< 0] <- 0
    Genus_meth_pos_graph <- graph.adjacency(Genus_meth_pos_cors, mode ='undirected', weighted = TRUE)
    Genus_meth_pos_graph <- simplify(Genus_meth_pos_graph)
    plot(Genus_meth_pos_graph)
    plot(Genus_meth_pos_graph, rescale = TRUE, edge.color = ifelse(E(Genus_meth_pos_graph)$weight> 0, "pink", "blue"), vertex.label.font = 3, vertex.size = deg*.8, vertex.label.cex = 1)
    
    #plot of only negative co-occurrences
    Genus_meth_neg_cors <- Genus_meth_cors_pruned
    Genus_meth_neg_cors[Genus_meth_neg_cors > 0] <- 0
    Genus_meth_neg_graph <- graph.adjacency(Genus_meth_neg_cors, mode ='undirected', weighted = TRUE)
    Genus_meth_pos_graph <- simplify(Genus_meth_neg_graph)
    plot(Genus_meth_neg_graph)
    plot(Genus_meth_neg_graph, rescale = TRUE, edge.color = ifelse(E(Genus_meth_pos_graph)$weight> 0, "pink", "blue"), vertex.label.font = 3, vertex.size = deg*.8, vertex.label.cex = 1)
    
  #export pos cooccur plot    
    svglite(filename = "meth_pos_cooccur.svg", width =10, height =10)
    plot(Genus_meth_pos_graph, rescale = TRUE, edge.color = ifelse(E(Genus_meth_pos_graph)$weight> 0, "pink", "blue"), vertex.label.font = 3, vertex.size = deg*.8, vertex.label.cex = 1)
    dev.off()    
    
  #regular plot
    A <- layout_with_fr(g_Genus_meth)
    plot(g_Genus_meth, rescale = TRUE, edge.color = ifelse(E(g_Genus_meth)$weight >0, "lightpink", "lightblue"), vertex.label.font = 3, layout = A, vertex.size = deg*.1, vertex.label.cex = deg*.01)
    
############################################################################################################################
    
# Visualize Correlation Matrix: Genus: Meth
    
############################################################################################################################       
    
    deg <- degree(g_Genus_Actino, mode = "all")      
    
    #Graph of all edges transformed into positive 
    Genus_Actino_cors_pos <- abs(Genus_Actino_cors_pruned)
    g_Genus_Actino_pos <- graph.adjacency(Genus_Actino_cors_pos, mode = "undirected", weighted= TRUE)
    g_Genus_Actino_pos <- simplify(g_Genus_Actino_pos)     
    plot(g_Genus_Actino_pos)
    
    #Graph of only positive co-occurrences
    Genus_Actino_pos_cors <- Genus_Actino_cors_pruned
    Genus_Actino_pos_cors[Genus_Actino_pos_cors< 0] <- 0
    Genus_Actino_pos_graph <- graph.adjacency(Genus_Actino_pos_cors, mode ='undirected', weighted = TRUE)
    Genus_Actino_pos_graph <- simplify(Genus_Actino_pos_graph)
    plot(Genus_Actino_pos_graph)
    plot(Genus_Actino_pos_graph, rescale = TRUE, edge.color = ifelse(E(Genus_Actino_pos_graph)$weight> 0, "pink", "blue"), vertex.label.font = 3, vertex.size = deg*.4, vertex.label.cex = 1)
    
    #plot of only negative co-occurrences
    Genus_Actino_neg_cors <- Genus_Actino_cors_pruned
    Genus_Actino_neg_cors[Genus_Actino_neg_cors > 0] <- 0
    Genus_Actino_neg_graph <- graph.adjacency(Genus_Actino_neg_cors, mode ='undirected', weighted = TRUE)
    Genus_Actino_pos_graph <- simplify(Genus_Actino_neg_graph)
    plot(Genus_Actino_neg_graph)
    plot(Genus_Actino_pos_graph, rescale = TRUE, edge.color = ifelse(E(Genus_Actino_pos_graph)$weight> 0, "pink", "blue"), vertex.label.font = 3, vertex.size = deg*.4, vertex.label.cex = 1)
    
    #export pos cooccur plot    
    svglite(filename = "Actino_pos_cooccur.svg", width =10, height =10)
    plot(Genus_Actino_pos_graph, rescale = TRUE, edge.color = ifelse(E(Genus_Actino_pos_graph)$weight> 0, "pink", "blue"), vertex.label.font = 3, vertex.size = deg*.4, vertex.label.cex = 1)
    dev.off()    
    
    #regular plot
    A <- layout_with_fr(g_Genus_Actino)
    plot(g_Genus_Actino, rescale = TRUE, edge.color = ifelse(E(g_Genus_Actino)$weight >0, "lightpink", "lightblue"), vertex.label.font = 3, layout = A, vertex.size = deg*.4, vertex.label.cex = 1)
    
    #export pos cooccur plot    
    svglite(filename = "Actino_cooccur.svg", width =10, height =10)
    plot(g_Genus_Actino, rescale = TRUE, edge.color = ifelse(E(g_Genus_Actino)$weight >0, "lightpink", "lightblue"), vertex.label.font = 3, layout = A, vertex.size = deg*.4, vertex.label.cex = 1)
    dev.off()    
    
    
############################################################################################################################
    
# Visualize Correlation Matrix: Species
    
############################################################################################################################    
    
#Species
    #Graph of all edges transformed into positive 
    Species_cors_pos <- abs(Species_cors_pruned)
    g_Species_pos <- graph.adjacency(Species_cors_pos, mode = "undirected", weighted= TRUE)
    g_Species_pos <- simplify(g_Species_pos)     
    plot(g_Species_pos)
    
    #Graph of only positive co-occurrences
    Species_pos_cors <- Species_cors_pruned
    Species_pos_cors[Species_pos_cors< 0] <- 0
    Species_pos_graph <- graph.adjacency(Species_pos_cors, mode ='undirected', weighted = TRUE)
    Species_pos_graph <- simplify(Species_pos_graph)
    plot(Species_pos_graph)
    plot(Species_pos_graph, rescale = TRUE, edge.width = abs(E(g_Species)$weight)*.5, edge.color = ifelse(E(g_Species)$weight> 0, "lightpink", "lightblue"), vertex.label.font = 3, vertex.size = deg*.1, vertex.label.cex = deg*.01)
    
############################################################################################################################
    
    # Visualize Correlation Matrix: Genus: 70 and minfreq 5000
    
############################################################################################################################       
    
    deg <- degree(genus_samp70_Actino, mode = "all")      
    
    #Graph of all edges transformed into positive 
    g_genus_samp70_cors_pos <- abs(genus_samp70_cors_pruned)
    g_genus_samp70_pos <- graph.adjacency(genus_samp70_cors_pos, mode = "undirected", weighted= TRUE)
    g_genus_samp70_pos <- simplify(g_genus_samp70_pos)     
    plot(g_genus_samp70_pos)
    
    #Graph of only positive co-occurrences
    genus_samp70_pos_cors <- genus_samp70_cors_pruned
    genus_samp70_pos_cors[genus_samp70_pos_cors< 0] <- 0
    genus_samp70_pos_graph <- graph.adjacency(genus_samp70_pos_cors, mode ='undirected', weighted = TRUE)
    genus_samp70_pos_graph <- simplify(genus_samp70_pos_graph)
    plot(genus_samp70_pos_graph)
    plot(genus_samp70_pos_graph, rescale = TRUE, edge.color = ifelse(E(genus_samp70_pos_graph)$weight> 0, "pink", "blue"), vertex.label.font = 3, vertex.size = deg*.1, vertex.label.cex = .75)
    
    #plot of only negative co-occurrences
    genus_samp70_neg_cors <- genus_samp70_cors_pruned
    genus_samp70_neg_cors[genus_samp70_neg_cors > 0] <- 0
    genus_samp70_neg_graph <- graph.adjacency(genus_samp70_neg_cors, mode ='undirected', weighted = TRUE)
    genus_samp70_pos_graph <- simplify(genus_samp70_neg_graph)
    plot(genus_samp70_neg_graph)
    plot(genus_samp70_neg_graph, rescale = TRUE, edge.color = ifelse(E(genus_samp70_neg_graph)$weight> 0, "pink", "blue"), vertex.label.font = 3, vertex.size = deg*.4, vertex.label.cex = 1)
    
    #export pos cooccur plot    
    svglite(filename = "All_pos_cooccur_samp70.svg", width =10, height =10)
    plot(genus_samp70_pos_graph, rescale = TRUE, edge.color = ifelse(E(genus_samp70_pos_graph)$weight> 0, "pink", "blue"), vertex.label.font = 3, vertex.size = deg*.1, vertex.label.cex = .75)
    dev.off()    
    
    #regular plot
    A <- layout_with_fr(g_genus_samp70)
    plot(g_genus_samp70, rescale = TRUE, edge.color = ifelse(E(g_genus_samp70)$weight >0, "lightpink", "lightblue"), vertex.label.font = 3, layout = A, vertex.size = deg*.1, vertex.label.cex = .75)
    
    #export pos cooccur plot    
    svglite(filename = "All_cooccur_samp70.svg", width =10, height =10)
    plot(g_genus_samp70, rescale = TRUE, edge.color = ifelse(E(g_genus_samp70)$weight >0, "lightpink", "lightblue"), vertex.label.font = 3, layout = A, vertex.size = deg*.1, vertex.label.cex = .75)
    dev.off()    
    
############################################################################################################################

    
    
    
    
    
        
############################################################################################################################

# Exporting from R

############################################################################################################################

#Export of regular plot
svglite(filename = "plot.svg", width =10, height =10)
plot(g_Genus, rescale = TRUE, edge.width = abs(E(g_Genus)$weight)*1, edge.color = ifelse(E(g_Genus)$weight >0, "lightpink", "lightblue"), vertex.label.font = 3)
dev.off()

svglite(filename = "plot.svg", width =10, height =10)
plot(g_Genus, rescale = TRUE, edge.width = abs(E(g_Genus)$weight)*.5, edge.color = ifelse(E(g_Genus)$weight >0, "lightpink", "lightblue"), vertex.label.font = 3, layout = A, vertex.size = deg*.1, vertex.label.cex= ifelse(deg>55, .01*deg, .001))
dev.off()

svglite(filename = "plot.svg", width =10, height =10)
plot(Genus_pos_graph, rescale = TRUE, edge.width = abs(E(g_Genus)$weight)*.5, edge.color = ifelse(E(g_Genus)$weight> 0, "lightpink", "lightblue"), vertex.label.font = 3, vertex.size = deg*.1, vertex.label.cex = deg*.01)
dev.off()
