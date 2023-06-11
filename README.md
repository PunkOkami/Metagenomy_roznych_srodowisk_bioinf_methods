# Metagenomy_roznych_srodowisk_bioinf_methods
Bioinformatical methods and tools used in 'Microbiomes of deciduous trees prospective for the forests on the example of small-leaved linden (_Tilia cordata_, L.) - bioinformatic-phenomic study of seedlings.' and presented at "Metagenomy różnych środowisk” Symposium in 2023

![Metagenomy23 - Poster w A3]("Graphs/Metagenomy23_poster.gif")


## Overview and description
The project aimed to profile mycobiomes of small-leaved linden (_Tilia cordata_, L.). Methodolody consists of metabarcoding, OTU clustering, taxonomy assignment, ecological role assignment, pairwise and k-means clustering of samples, computational statistical validation, trophic network construction and analysis of constructed networks.

### Sampling
Soil samples of two groups of trees. One from nursery pots and one from trees planted in the forest a year ago. 

### Sequencing
Sequencing was perfomed using IONTorrent platform and ITS1 marker region. 

### Quality filtering and OTU clustering
Quality filtering and OTU clustering was achieved using Qiime2 tool and used commands are listed in file `commands/Qiime_command_used.md`. 

### Taxonomic assigment
Taxonomic assigment was done using local BLASTN searching against fungal UNITE database. OTUs that were not found in fungal databse were filtered by count with cut off 10 occurences and then searched for in UNITE databse of all genomes. That gave us 5 sequences and 4 were found in data base, all being single cell Eukaria species outside of funga. Next we combined abundance of OTUs assigned to the same species into one record. That gave us 820 species of fungi not counting singletons.

### Guilds and ecological metacategories
Using taxonomic information obtained from UNITE database FUNGuild tool was used to assign guilds. To assure confidence in guild information we kept only data labeled as 'Highly Probable". Later we proposed metacategories based on hierarchical order based on ecological effect on plants: Myc (mycorrhizal of all kinds), Pat (plant pathogen), Sap (saprotrophs of all kinds), End (endophytes), Bcn (biocontrol) and Nop (non-plant specific). Except those 6 we named the rest as ENA (ecologically not assigned).

### Occurance filtering
At this stage we used an abundance filter. We kept only the species that togther constitute over 99% of all hits. The filter left us woth 217 taxa that make up 99.26% of all hits and limited ecological metacategories present to %: Myc, Pat, Sap, Nop and ENA. The make up of taxa and groups can be found in Krona charts and heatmaps in folder Graphsas well as in tsv tables in Data folder.

### Clustering
Parallel to occurance filtering our pipeline perfomes pair-wise clustering of samples and produces a clustermap showing the pairwise clustering and biodiversity of samples. K-means clustering is based on two biodiversity indexes: Margalef richness index and Simpson eveness index. Both indexes are also validated using boostrap method.

### Network construction and analysis
To construct trophic networks we use NetCoMi R package. We used SpiecEasi method to construct networks due to relibly high modularity of constructed networks (over 0.36). We analysed networks usign graph parameters such as clustering coefficient, modularity and natural connectivity. We also analysed clusters that were created in said networks to deterine what metacategories make up those clusters.

## Authors
Sampling and sequncing: Mikołaj Charchuta - mikcha1@st.amu.edu.pl\
Bioinformatics work: Maksymilian Chmielewski - makchm@st.amu.edu.pl\
Supervisor: Prof. UAM dr hab. Władysław Polcyn - polcyn@amu.edu.pl

## Affiliations
Bioinformatics Section of Natural Sciences Club\
Adam Mickiwicz University in Poznań
