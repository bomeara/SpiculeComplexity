

getGenus <- function(x) {
  genera <- rep(NA, length(x)) 
  for(i in seq_along(x)) {
    genera[i] <- strsplit(x[i], "_")[[1]][1]
  }
  
  return(genera)
}

badLabels <- function(x) {
  x[is.na(x)] <- ""
  badones <- which(nchar(x)==0)
  return(badones)
}

getComplexity <- function() {
	complexity <- read.csv("../Matrix/Complexity.csv")
	complexity$Species <- gsub(" ", "_", complexity$Species)
	return(complexity)
}

getRaxmlTree <- function() {
	phy <- ape::read.tree("raxml.tre")[[1]]
	root_node <- ape::getMRCA(phy, c("Corticium_vaceleti", "Pheronemoides_crustiformis"))
	phy <- ape::root(phy, node=root_node)
	return(phy)
}

dateRaxmlTree <- function(phy) {
	phy <- compute.brlen(phy, 1)
	phy$edge.length <- 600*phy$edge.length/max(phytools::nodeHeights(phy)[,2])
	return(ape::chronos(phy, calibration=data.frame(node=ape::Ntip(phy)+1, age.min=600, age.max=600, soft.bounds=FALSE)))
}


getPaleoDBTree <- function() {
	paleodb_phy <- paleotree::makePBDBtaxonTree(paleotree::getCladeTaxaPBDB("Porifera"), rankTaxon = "species")
	return(paleodb_phy)
}


getCharSpecies <- function(complexity) {
	charSpecies <- complexity$Species
	return(charSpecies)
}

getCharGenus <- function(complexity) {
	charSpecies <- getCharSpecies(complexity)
	charGenus <- getGenus(charSpecies)

	return(charGenus)
}

getOTTTrees <- function(charSpecies, charGenus) {
	ott_taxa <- tnrs_match_names(charSpecies, context = "Animals")
	ott_tree <- tol_induced_subtree(ott_id(ott_taxa)[is_in_tree(ott_id(ott_taxa))], label_format="name")
	missing_in_ott_tree <- charSpecies[!charSpecies %in% ott_tree$tip.label]
	genera_missing_in_ott_tree <- charGenus[!charGenus %in% getGenus(ott_tree$tip.label)]

	ott_taxa2 <- tnrs_match_names(c(charSpecies, getGenus(missing_in_ott_tree)), context = "Animals")

	ott_tree2 <- tol_induced_subtree(ott_id(ott_taxa2)[is_in_tree(ott_id(ott_taxa2))], label_format="name")

	ott_taxonomy <- taxonomy_subtree(ott_id=ott_id(tnrs_match_names("Porifera")), label_format="name", output_format="newick")
	ott_taxonomy_tree <- ape::read.tree(text=ott_taxonomy)
	
	
	ott_taxonomy_tree_no_duplicates <- ott_taxonomy_tree
	ott_taxonomy_tree_no_duplicates <- ape::drop.tip(ott_taxonomy_tree_no_duplicates, tip=sequence(length(ott_taxonomy_tree_no_duplicates$tip.label))[duplicated(ott_taxonomy_tree_no_duplicates$tip.label)])
	badones <- badLabels(ott_taxonomy_tree_no_duplicates$tip.label)
	if(length(badones)>0) {
	ott_taxonomy_tree_no_duplicates <- ape::drop.tip(ott_taxonomy_tree_no_duplicates, tip=badones)
	}
	
	return(list(missing_in_ott_tree=missing_in_ott_tree, genera_missing_in_ott_tree=genera_missing_in_ott_tree, ott_tree=ott_tree, ott_tree2=ott_tree2, ott_taxonomy_tree=ott_taxonomy_tree, ott_taxonomy_tree_no_duplicates=ott_taxonomy_tree_no_duplicates))
}

getCleanedRaxmlTreeData <- function(charSpecies, raxml_phy, charGenus, missing_in_ott_tree, complexity, genera_missing_in_ott_tree) {
	phy <- raxml_phy
	missing_in_raxml_tree <- charSpecies[!charSpecies %in% phy$tip.label]
	genera_missing_in_raxml_tree <- charGenus[!charGenus %in% getGenus(phy$tip.label)]
	missing_in_both_trees <- intersect(missing_in_raxml_tree, missing_in_ott_tree)
	genera_missing_in_both_trees <- intersect(genera_missing_in_raxml_tree, genera_missing_in_ott_tree)



	matching_species <- complexity$Species[complexity$Species %in% phy$tip.label]
	#unmatching_species <- complexity$Species[!complexity$Species %in% phy$tip.label]
	#genera_of_unmatching_species <- getGenus(unmatching_species)

	for (i in seq_along(complexity$Species)) {
	if(!(complexity$Species[i] %in% matching_species)) {
		complexity$Species[i] <- getGenus(complexity$Species[i])
	}
	}

	complexity_vector <- complexity$SpiculeTypes
	names(complexity_vector) <- complexity$Species

	for(i in seq_along(phy$tip.label)) {
		if(!(phy$tip.label[i] %in% matching_species)) {
			if(!any(phy$tip.label==getGenus(phy$tip.label[i]))) {
			phy$tip.label[i] <- getGenus(phy$tip.label[i])
			}
		}
	}

	cleaned <- geiger::treedata(phy, complexity_vector)
	cleaned$data <- cleaned$data[!duplicated(rownames(cleaned$data)),]
	return(cleaned)
}

plotComplexity <- function(cleaned_raxml) {
	cleaned <- cleaned_raxml
	phytools::contMap(ape::compute.brlen(cleaned$phy), cleaned$data)

}


# phy_no_duplicates <- phy
# phy_no_duplicates <- ape::drop.tip(phy_no_duplicates, tip=sequence(length(phy_no_duplicates$tip.label))[duplicated(ott_taxonomy_tree_no_duplicates$tip.label)])
# badones <- badLabels(phy_no_duplicates$tip.label)
# if(length(badones)>0) {
#   phy_no_duplicates <- ape::drop.tip(phy_no_duplicates, tip=badones)
# }

#super <- phangorn::superTree(c(phy_no_duplicates, ott_taxonomy_tree_no_duplicates), trace=4)


#cleaned_super <- geiger::treedata(super, complexity_vector)

#phytools::contMap(ape::compute.brlen(cleaned_super$phy), cleaned_super$data)

MakeConstraints <- function(fossil_phy) {
	constraints <- data.frame()
	unique_nodes <- unique(fossil_phy$edge[,1])
	if(is.null(fossil_phy$node.label)) {
		fossil_phy$node.label <- fossil_phy$edge[,1]
	}
	for (i in sequence(length(unique_nodes))) {
		descendant_numbers <- phangorn::Descendants(fossil_phy, unique_nodes[i], type="tips")[[1]]
		descendant_names <- fossil_phy$tip.label[descendant_numbers]
		non_descendant_names <- fossil_phy$tip.label[!(sequence(ape::Ntip(fossil_phy)) %in% descendant_numbers)]
		constraints <- rbind(constraints, data.frame(node=unique_nodes[1], name=fossil_phy$node.label[i], descendants=paste(descendant_names, collapse=" "), nondescendants=paste(non_descendant_names, collapse=" ")))
	}	
	constraints <- constraints[nchar(constraints$nondescendants)>0,]
	return(constraints)
}

PrintConstraints <- function(constraints, output_file="constraints.nex") {
	cat("#NEXUS\n\nbegin mrbayes;\n", file=output_file, append=FALSE)

	for (i in sequence(nrow(constraints))) {
		cat("\nconstraint backbone partial = ", constraints$descendants[i], " : ", constraints$nondescendants[i], file=output_file, append=TRUE)	
	}
	cat("\nend;\n", file=output_file, append=TRUE)
}

GenerateTipAges <- function(paleodb_ages) {
	species <- subset(paleodb_ages, taxon_rank=="species")
	return(species[, c("taxon_name", "lastapp_min_ma")])
}

MakePaleoDBTreeDated <- function(paleodb_ages, paleodb_tree) {
	species <- subset(paleodb_ages, taxon_rank=="species")
	species$taxon_name <- gsub(" ", "_", species$taxon_name)
	species <- species[!is.na(species$lastapp_min_ma),]
	calibrations <- data.frame()
	unmatched_taxa <- paleodb_tree$tip.label[!(paleodb_tree$tip.label %in% species$taxon_name)]
	paleodb_tree <- ape::drop.tip(paleodb_tree, tip=unmatched_taxa)
	for (i in sequence(length(paleodb_tree$tip.label))) {
		focal_taxon <- paleodb_tree$tip.label[i]
		row_number <- which(species$taxon_name==focal_taxon)[1]
		if(length(row_number)==1) {
			calibrations <- rbind(calibrations, data.frame(node=i, age.min=species$lastapp_min_ma[row_number], age.max=species$firstapp_max_ma[row_number]))
		} else {
			print(length(row_number))
			print(paste0("No calibration for ", focal_taxon))
		}
	}
	#calibrations <- calibrations[!is.na(calibrations$age.min),]
	#calibrations <- calibrations[calibrations$age.min>1e-8,]

	#calibrations <- rbind(calibrations, data.frame(node=Ntip(paleodb_tree)+1, age.min=50+max(calibrations$age.max), age.max=100+max(calibrations$age.max)))
	paleodb_tree <- compute.brlen(paleodb_tree, method=1)
	paleodb_tree$edge.length <- paleodb_tree$edge.length * max(calibrations$age.max ) / max(phytools::nodeHeights(paleodb_tree)[,2])
	
	
	node.date <- estimate.dates(paleodb_tree, calibrations$age.min)

     ## To rescale the tree over time
     paleodb_tree$edge.length <- node.date[paleodb_tree$edge[, 2]] - node.date[paleodb_tree$edge[, 1]]
	 paleodb_tree$edge.length <- 600*paleodb_tree$edge.length / max(phytools::nodeHeights(paleodb_tree)[,2])
	#chronogram <- chronos(paleodb_tree, calibration=calibrations)
	return(paleodb_tree)

}

PrintTipAges <- function(species_ages) {
	cat("#NEXUS\n\nbegin mrbayes;\ncalibrate ", file="tip_ages.nex", append=FALSE)
	for (i in sequence(nrow(species_ages))) {
		cat(gsub(" ", "_", species_ages$taxon_name[i]), " = Fixed(", species_ages$lastapp_min_ma[i], ")\n", file="tip_ages.nex", append=TRUE)	
	}
	cat("\n;	prset brlenspr=clock:uniform;
	prset clockvarpr=igr;	
	prset igrvarpr=exp(37.12);
	prset clockratepr = lognorm(-7.08069,2.458582);
	calibrate root=offsetexp(315,0.01234568);
	calibrate holometabola_with_fossils=offsetexp(302,0.0106383);
	prset topologypr=constraints(root, holometabola_with_fossils);
	prset nodeagepr = calibrated;
	end;", file="tip_ages.nex", append=TRUE)

}

PrintMrBayesBatchFile <- function() {
	cat("#NEXUS\n\nbegin mrbayes;
	execute combined_only_relevant_aligned_pruned_for_mb.nex;
	execute constraints.nex;
	execute tip_ages.nex;
	mcmc ngen=20000000 temp=0.1 nchain=4 samplefreq=1000 printfr=100 nruns=2 filename=spongeDating;	
	end;", file="spongeDating.nex", append=FALSE)
}

RunMrBayes <- function(...) {
	system("mb spongeDating.nex")
}

GetTreeDataTree <- function(x) { return(x$phy) }

GetTreeDataData <- function(x) { return(x$data) }

PrintContMap <- function(phy, data, output_file="contmap.pdf", complexity=complexity) {
	pdf(file=output_file, width=16, height=16)
	phytools::contMap(phy, data, lims=range(complexity$SpiculeTypes))
	dev.off()	
}

ComputeRates <- function(phy, data) {
	BM <- geiger::fitContinuous(phy=ape::multi2di(phy), dat=data, model="BM", control=list(niter = 100, FAIL = 1e+200, hessian=TRUE))
	return(BM$opt)
}

SaveBMOutput <- function(raxml_BM, paleodb_BM) {
	write.csv(cbind(raxml_BM, paleodb_BM), file="BM_output.csv")	
}

PruneCombinedToRelevantTaxa <- function(charGenus, charSpecies) {
	unique_genera <- unique(charGenus)
	cat("\n", file="combined_only_relevant.fasta", append=FALSE)
	for (i in sequence(length(unique_genera))) {
		print(paste0("Trying ", unique_genera[i]))
		system(paste0("grep -A1 ", unique_genera[i],  " combined_not_interleaved.fasta >> combined_only_relevant.fasta"))
	}
	system("mafft-einsi combined_only_relevant.fasta > combined_only_relevant_aligned.fasta")
	dna <- as.character(ape::read.FASTA(file="combined_only_relevant_aligned.fasta"))
	dna_matrix <- t(simplify2array(dna))
	duplicated_taxa <- rownames(dna_matrix)[duplicated(rownames(dna_matrix))]
	rownames(dna_matrix)[which(duplicated(rownames(dna_matrix)))] <- paste0(duplicated_taxa, " 2")
	missing_taxa <- setdiff(gsub("_", " ", charSpecies), rownames(dna_matrix))
	dna_with_missing <- rbind(dna_matrix, matrix("-", nrow=length(missing_taxa), ncol=length(dna_matrix[1,])))
	rownames(dna_with_missing) <- c(rownames(dna_matrix), missing_taxa)
	number_missing_states <- rep(NA, ncol(dna_with_missing))
	for (i in sequence(length(number_missing_states))) {
		number_missing_states[i] <- sum(dna_with_missing[,i]=="-")
	}
	number_having_states <- nrow(dna_with_missing) - number_missing_states
	dna_with_missing <- dna_with_missing[,number_having_states>=4]
	rownames(dna_with_missing) <- gsub(" ", "_", rownames(dna_with_missing))
	dna_DNAbin <- ape::as.DNAbin(dna_with_missing)
	ape::write.nexus.data(dna_DNAbin, file="combined_only_relevant_aligned_pruned.nex")
	print("now open this in mesquite and save as combined_only_relevant_aligned_pruned_for_mb.nex")
}