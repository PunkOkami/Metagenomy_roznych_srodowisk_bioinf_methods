library(NetCoMi)

# Reading data and preparing it for network costruction with filtering out all taxa not found in set metacategories
data = read.delim('../../Data/species_freq_table_filtered.tsv', sep='\t')
guild_we_want = c('Pat')
data = data[data$Guild.category %in% guild_we_want,]
otu_names = data[,1]
data = data[,-c(1, 16:30)]
data = as.matrix(data)
rownames(data) = otu_names
data = apply(t(data), 2, rev)

# Basic single network costruction and analysis to set it ready to plot. Values are set to keep around 100 species
net <- netConstruct(data,
                    filtTax = "totalReads",
                    filtTaxPar = list(totalReads = 200),
                    filtSamp = "totalReads",
                    filtSampPar = list(totalReads = 100),
                    measure = "spieceasi",
		    normMethod = "none",
                    zeroMethod = "none",
                    sparsMethod = "none",
                    dissFunc = "signed",
                    verbose = 3,
                    seed = 123456)

props <- netAnalyze(net_spring,
                    centrLCC = TRUE,
                    clustMethod = "cluster_fast_greedy",
                    hubPar = "eigenvector",
                    weightDeg = FALSE, normDeg = FALSE)

network_title = paste('Species in guilds:', paste(guild_we_want, collapse = ', '), 'over 200 reads')
p <- plot(props_spring,
          nodeColor = "cluster",
          nodeSize = "eigenvector",
          title1 = network_title,
          showTitle = TRUE,
          cexTitle = 1.9,
          repulsion = 0.95,
          labelScale = FALSE,
          cexLabels = 1,
          cexHubLabels = 1.5)

legend(0.4, 1.09, cex = 1.5, title = "estimated association:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"),
       bty = "n")
