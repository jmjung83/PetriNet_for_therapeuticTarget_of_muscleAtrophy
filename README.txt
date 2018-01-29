=== Read me file ===

[Title]
Petri net-based prediction of therapeutic targets that recovers abnormally phosphorylated proteins in muscle atrophy (Jung et al)

[Environment]
 python-2.7.13 (Anaconda2-5.0.1)
 gseapy package

[Exclustion of the simulation results that need large space]
o Excluded files
  1) ./REFR/Sim_result/*
  2) ./MYST/Sim_result/*
  3) ./GeneInhibit/Sim_result_1/*
  4) ./GeneInhibit/Sim_result_2/*
  5) ./GeneInhibit/Sim_result_3/*

[Order of the program]
1. network construction (i.e. myocyte specific phosphorylation network)
 o code
    ./network/network_construction.py
 o input files
    1) ./network/KEGG/phosDeph_network(KEGG).txt
    2) ./network/PhosphoNet/phos_network(KSI).txt
    3) ./network/DEPOD/deph_network(DEPOD).txt
    4) ./network/HPA_protein_expression (BSML_gene).tsv
 o output files
    1) ./network/transition.txt (network edge information)
    2) ./network/place.txt (network node information)

2. Petrinet simulation for reference state (10 simulations for each of 11 enabling thresholds (TH))
 o code: ./REFR/v1_CODA_simulator_S.py
 o input files
    1) ./network/transition.txt
    2) ./network/place.txt
 o output files
    1) ./REFR/Sim_result/*

3. Petrinet simulation for atrophic state (10 simulations for each of 11 enabling thresholds (TH))
 o code: ./MYST/v1_CODA_simulator_S.py
 o input files
    1) ./network/transition.txt
    2) ./network/place.txt
 o output files
    1) ./MYST/Sim_result/*

4. Compute rank correlation of 10 final marking for each of 11 enabling thresholds in atrophic and referenc states
 o code: 1_get_finalMarking_RankCorr(REFR_MYST).py
 o input files
    1) ./REFR/Sim_result/*
    2) ./MYST/Sim_result/*
 o output files
    1) printed output (rank correlations in atrophic and reference states)
    2) ./REFR/final_marking_bioMarker(TH35).txt (final marking in referece state for the highst correlation threshold(35))
    3) ./MYST/final_marking_bioMarker(TH35).txt (final marking in atrophic state for the highst correlation threshold(35))

5. Compute final marking difference between atrophic and reference state (which is used for validation with Western blot results)
 o code: 2_get_finalMarkingDiff(REFR_MYST).py
 o input files
    1) ./REFR/final_marking_bioMarker(TH35).txt
    2) ./MYST/final_marking_bioMarker(TH35).txt
 o output files (simulation results for various enabling threshold, TH)
    1) printed output (difference, five among six proteins are consistent with Western Blot in terms of directions of difference)

6. Petrinet simulation for gene-inhibited state (10 simulations for each of 331 proteins connected to the five validated markers)
 o code: ./GeneInhibit/v1_CODA_simulator_S.py (three times)
 o input files
    1) ./network/transition.txt
    2) ./network/place.txt
 o output files
    1) ./GeneInhibit/Sim_result_1/*
    2) ./GeneInhibit/Sim_result_2/*
    3) ./GeneInhibit/Sim_result_3/*

7. Compute final markings of the five validated markers in gene-inhibited states
 o code: 3_get_finalMarking(GeneInhibit).py
 o input files
    1) ./GeneInhibit/Sim_result_1/*
    2) ./GeneInhibit/Sim_result_2/*
    3) ./GeneInhibit/Sim_result_3/*
 o output files
    1) ./GeneInhibit/final_marking_bioMarker(TH35)_1time.txt
    2) ./GeneInhibit/final_marking_bioMarker(TH35)_2time.txt
    3) ./GeneInhibit/final_marking_bioMarker(TH35)_3time.txt

8. Compute candidiate therapeutic targets by comparing directions of differecce in the five markers
 o code: 4_get_therapeuticTarget.py
 o input files
    1) ./REFR/final_marking_bioMarker(TH35).txt
    2) ./MYST/final_marking_bioMarker(TH35).txt
    3) ./GeneInhibit/final_marking_bioMarker(TH35)_1time.txt
    4) ./GeneInhibit/final_marking_bioMarker(TH35)_2time.txt
    5) ./GeneInhibit/final_marking_bioMarker(TH35)_3time.txt
 o output files
    1) all_proteins.txt (all 331 considered proteins)
    2) target_candidates.txt (37 therapeutic targets satisfying the predefined conditioin)

9. gene enrichment test for the 37 predicted targets
 o requirement: gseapy package
 o code: 5_geneEnrichTest.py
 o input files
    1) None
 o output files
    1) ./enrichedGOTerms/GO_Biological_Process_2015.BP2015.enrichr.reports.txt
    2) ./enrichedGOTerms/KEGG_2016.KEGG2016.enrichr.reports.txt

10. draw a plot of simulation results for PTPRS and SMAD3
 o code: drawPlot.py
 o input files
    1) ./GeneInhibit/Sim_result_1/5802_PTPRS/*
    2) ./GeneInhibit/Sim_result_2/5802_PTPRS/*
    3) ./GeneInhibit/Sim_result_3/5802_PTPRS/*
    4) ./REFR/Sim_result/TH35/*
    4) ./MYST/Sim_result/TH35/*
 o output files
    1) plot

11. evaluation with muscle atrophy related genes from CTD
 o code: 7_CTD_eval.py
 o input files
    1) ./CTD/CTD_genes_related_to_atrophy.txt
    2) target_candidates.txt
 o output files
    1) printed output (hypergeometric test results)

12. check network inconsistency
 o code: 8_networkInconsistency.py
 o input files
    1) ./network/transition.txt
 o output files
    1) printed output (the results of inconsistency)
=========================================================  End of file ========================================