# Enriched citrus vascular proteomics highlights the role of peroxidases and serine protease during Huanglongbing disease progression

Jessica T. Franco, Shree P. Thapa, Zhiqian Pang, Fatta B. Gurung, Thomas Liebrand, Danielle M. Stevens, Veronica Ancona, Nian Wang, Gitta Coaker.

-----------------------

Purpose: The script in this repository is for processing and ploting z-scores from total proteomic data as well as Go-terms (Gene Ontology). If you have any questions, concerns, or found a bug, please contact me at dmstev@ucdavis.edu.



Due to the redundacy of go-terms called for each protein present, we had to first make a custom script to:
 1) bin the terms by the protein called
 2) plotted to see which terms are most abundant
 3) manually selected based upon the annotation of the protein
 
This list is then writen out to the file. This output (labeled Filtered_GO_terms_to_plot.xlsx) is then inputed, organized and plotted (via ggplot). Specifically, it:
 1) bins the number of proteins based upon their desginated go-terms (molecular function and biological proccess) 
 2) plots weather the term was significantly significantly downregulated or upregulated
 
 Finally, normalized z-scores were clustered and plotted with additonal meta data collected in the study.
 
 Information regarding R and package versions can be found in the session_info.txt document. 
 
 
 ***There is one known problem with the script. One occasion, I have had issues wth the installation and loading of the ComplexHeatmap package. I suspect this may be due to myself working on a windows computer. The current fix I have is reinstalling the package, soemtimes using different methods of installation, which one can see in the script. This should fix the issue and allow you to run the script.
