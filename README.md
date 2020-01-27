# Franco_2019_Proteomics_Paper

The script in this repository are for processing and ploting z-scores from total proteomic data as well as Go-terms (Gene Ontology). If you have any questions, concerns, or found a bug, please contact me at dmstev@ucdavis.edu.


-----------------------

Do to the redundacy of go-terms called for each protein present, we had to first make a custom script to:
 1) bin the terms by the protein called
 2) plotted to see which terms are most abundant
 3) manually selected based upon the annotation of the protein
 
This list is then writen out to the file. Additonally, these proteins were further:
 1) binned based upon their desginated go-terms (molecular function and biological proccess) 
 2) labeled and subsequently plotted upon weather the term was significantly downregulated or upregulated
 
 Finally, normalized z-scores were clustered and plotted with additonal meta data (if the sample was from infected or healthy tissue).
 
 Information regarding package versions can be found in the session_info.txt document.
