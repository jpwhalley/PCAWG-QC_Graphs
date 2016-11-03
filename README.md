# PCAWG-QC_Graphs
Python 2.7.11 code used in producing graphs for the PCAWQ-QC paper. Packages required are numpy version 1.10.4, scipy version 0.17.0 and matplotlib version 1.5.1. Other versions may work, though these are the ones I used.

Imports the QC measures in the form of a tsv (Supplementary Table 1 in the PCAWG-QC paper), calculates which pass for each QC measure and gives a star rating. Various graphs to show various quality measures are also plotted.
		
INPUT: TSV files containing the data, metadata file linking the projects to tissure types. Numerical calculations are done using the numpy and scipy.stats packages. Graphs are plotted using matplotlib.pyplot package.

FUNCTION: star_rating(data)

OUTPUT: A TSV file with the star rating and a series of graphs used to illustrate the different QC measures and the star rating.

Command line to run (from same folder as the supplementary tables):
python -c 'import PCAWG_QC_Star_Rating; PCAWG_QC_Star_Rating.star_rating()'

Time Taken to run: ~ 30 seconds
