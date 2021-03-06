library(targets)
source("functions.R")


tar_option_set(packages=c("phytools", "geiger", "rotl", "phangorn", "paleotree"))



list(
 tar_target(complexity, getComplexity()),
 tar_target(raxml_tree, getRaxmlTree()),
 tar_target(paleodb_tree, getPaleoDBTree()),
 tar_target(charSpecies, getCharSpecies(complexity)),
 tar_target(charGenus, getCharGenus(complexity)),
 tar_target(ott_trees, getOTTTrees(charSpecies, charGenus)),
 tar_target(paleodb_ages, paleotree::getCladeTaxaPBDB("Porifera")),
 tar_target(raxml_treedata, getCleanedRaxmlTreeData(charSpecies=charSpecies, raxml_phy=raxml_tree, charGenus=charGenus, missing_in_ott_tree=ott_trees$missing_in_ott_tree, complexity=complexity, genera_missing_in_ott_tree=ott_trees$genera_missing_in_ott_tree)),
 tar_target(constraints, MakeConstraints(paleodb_tree)),
 tar_target(print_constraints, PrintConstraints(constraints)),
 tar_target(print_tip_constraints, PrintTipAges(GenerateTipAges(paleodb_ages) )),
 tar_target(print_mrbayes_batch_file, PrintMrBayesBatchFile()),
 tar_target(raxml_matched_tree, GetTreeDataTree(raxml_treedata)),
 tar_target(paleodb_matched_tree, GetTreeDataTree(paleodb_treedata)),
 tar_target(raxml_matched_data, GetTreeDataData(raxml_treedata)),
 tar_target(paleodb_matched_data, GetTreeDataData(paleodb_treedata)),
 tar_target(dated_raxml_tree, dateRaxmlTree(raxml_matched_tree)),
 tar_target(dated_paleodb_tree, MakePaleoDBTreeDated(paleodb_ages, paleodb_matched_tree)),
 tar_target(paleodb_treedata2, getCleanedRaxmlTreeData(charSpecies=charSpecies, raxml_phy=dated_paleodb_tree, charGenus=charGenus, missing_in_ott_tree=ott_trees$missing_in_ott_tree, complexity=complexity, genera_missing_in_ott_tree=ott_trees$genera_missing_in_ott_tree)),
 tar_target(dated_paleodb_tree2, GetTreeDataTree(paleodb_treedata2)),
 tar_target(paleodb_matched_data2, GetTreeDataData(paleodb_treedata2)),

 tar_target(paleodb_treedata, getCleanedRaxmlTreeData(charSpecies=charSpecies, raxml_phy=paleodb_tree, charGenus=charGenus, missing_in_ott_tree=ott_trees$missing_in_ott_tree, complexity=complexity, genera_missing_in_ott_tree=ott_trees$genera_missing_in_ott_tree)),
 tar_target(paleodb_contmap, PrintContMap(dated_paleodb_tree2, paleodb_matched_data2, "paleodb_contmap.pdf", complexity=complexity)),
 tar_target(raxml_contmap, PrintContMap(dated_raxml_tree, raxml_matched_data, "raxml_contmap.pdf", complexity=complexity)),
 tar_target(paleodb_BM, ComputeRates(dated_paleodb_tree2, paleodb_matched_data2)),
tar_target(raxml_BM, ComputeRates(dated_raxml_tree, raxml_matched_data)),
tar_target(bm_output, SaveBMOutput(raxml_BM, paleodb_BM)),
tar_target(clean_fasta, PruneCombinedToRelevantTaxa(charGenus))

 #tar_target(run_mb, RunMrBayes(print_constraints, print_tip_constraints, print_mrbayes_batch_file))
)