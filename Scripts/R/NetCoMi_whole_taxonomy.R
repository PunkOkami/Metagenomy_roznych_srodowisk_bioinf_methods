library(NetCoMi)
library(grDevices)

# Reads data and parses to fir NetCoMi input
data = read.delim('../../Data/species_freq_table.tsv', sep='\t')
otu_names = data[,1]
netcomi_data = data[,-c(1, 16:30)]
netcomi_data = as.matrix(netcomi_data)
rownames(netcomi_data) = otu_names
netcomi_data = apply(t(netcomi_data), 2, rev)

# Sets number of cores, remember to make it fit your machine
ncores = 12

# Calculates network using SpiecEasi using all 452 species that pass filtering
net_spieceasi <- netConstruct(netcomi_data,
                              filtTax = "totalReads",
                              filtTaxPar = list(totalReads = 75),
                              measure = "spieceasi",
                              normMethod = "none",
                              zeroMethod = "none",
                              sparsMethod = "bootstrap",
                              dissFunc = "signed",
                              verbose = 3,
                              cores = ncores)

# Analyses network and prints summary
props_spieceasi <- netAnalyze(net_spieceasi,
                           clustMethod = "cluster_fast_greedy",
                           normDeg = FALSE)
summary(props_spieceasi)

# Plots network and saves the plot obkect as source of colours
p_spec <- plot(props_spieceasi,
          nodeColor = "cluster",
          nodeSize = "eigenvector",
          title1 = "Clustering of T. cordata rizosphere funga, no hidden edges",
          showTitle = TRUE,
          cexTitle = 1.6,
          repulsion = 0.95,
          labelScale = FALSE,
          cexLabels = 0,
          cexHubLabels = 0,
          edgeInvisFilter = "threshold",
          edgeInvisPar = 0.15)

# Pulls out data about what clusters are made up of and splits it into sets of specie names according to cluster
spieceasi_clusters = props_spieceasi$clustering$clust1
spieceasi_clusters = split(spieceasi_clusters, unname(unlist(spieceasi_clusters)))
if (0 %in% unlist(unname(spieceasi_clusters))) {
  spieceasi_clusters = spieceasi_clusters[-1]
}
# Pulls out colours of clusters
colours = p_spec$nodecolor$nodecol1

# Creates a data.frame that have cluster name, guild names, guild proportions and colour of cluster
i = 1
spieceasi_cluster_data = data.frame()
for (cluster in spieceasi_clusters) {
  spiecies = names(cluster)
  colour = colours[spiecies]
  colour = unique((unname(unlist(colour))))
  guilds = data[which(data$Species %in% spiecies),]$Guild.category
  guild_f = factor(guilds)
  freqs = data[which(data$Species %in% spiecies),]$Total.abundance
  summed_freqs = c()
  all_freqs = sum(freqs)
  g_levels = levels(guild_f)
  for (level in g_levels) {
    sum_freq = sum(freqs[(guild_f %in% level)])/all_freqs
    summed_freqs = c(summed_freqs, sum_freq)
  }
  cluster_name = paste0('cluster_', i)
  g_levels = paste(g_levels, collapse = "_")
  summed_freqs = paste(summed_freqs, collapse = '_')
  colour = paste(colour, collapse = '_')
  i = i + 1
  row = c(cluster_name, g_levels, summed_freqs, colour)
  spieceasi_cluster_data = rbind(spieceasi_cluster_data, row)
}
# Saves data.frame  to tsv file
colnames(spieceasi_cluster_data) = c('cluster name', 'guilds in cluster', 'guilds percentiles', 'cluster colour')
write.table(spieceasi_cluster_data, '../../Data/spieceasi_cluster_data.tsv', sep = '\t', quote = FALSE, row.names = FALSE)


# This half of code does the same job as upper, but using SPRING as method of network construction
net_spring <- netConstruct(netcomi_data,
                           filtTax = "totalReads",
                           filtTaxPar = list(totalReads = 75),
                           measurePar = list(nlambda=100,
                                             rep.num=10),
                           measure = "spring",
                           normMethod = "none",
                           zeroMethod = "none",
                           sparsMethod = "none",
                           dissFunc = "signed",
                           verbose = 3,
                           cores = ncores)

props_spring <- netAnalyze(net_spring,
                           clustMethod = "cluster_fast_greedy",
                           normDeg = FALSE)
summary(props_spring)

p_spring <- plot(props_spring,
          nodeColor = "cluster",
          nodeSize = "eigenvector",
          title1 = "Species with over 100 sequences found",
          showTitle = TRUE,
          cexTitle = 1.9,
          repulsion = 0.95,
          labelScale = FALSE,
          cexLabels = 0,
          cexHubLabels = 0)

legend(0.4, 1.09, cex = 1.5, title = "estimated association:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"),
       bty = "n")

spring_clusters = props_spring$clustering$clust1
spring_clusters = split(spring_clusters, unname(unlist(spring_clusters)))
if (0 %in% unlist(unname(spring_clusters))) {
  spring_clusters = spring_clusters[-1]
}
colours = p_spring$nodecolor$nodecol1

i = 1
spring_cluster_data = data.frame()
for (cluster in spring_clusters) {
  spiecies = names(cluster)
  colour = colours[spiecies]
  colour = unique((unname(unlist(colour))))
  guilds = data[which(data$Species %in% spiecies),]$Guild.category
  guild_f = factor(guilds)
  freqs = data[which(data$Species %in% spiecies),]$Total.abundance
  summed_freqs = c()
  all_freqs = sum(freqs)
  g_levels = levels(guild_f)
  for (level in g_levels) {
    sum_freq = sum(freqs[(guild_f %in% level)])/all_freqs
    summed_freqs = c(summed_freqs, sum_freq)
  }
  cluster_name = paste0('cluster_', i)
  g_levels = paste(g_levels, collapse = "_")
  summed_freqs = paste(summed_freqs, collapse = '_')
  colour = paste(colour, collapse = '_')
  i = i + 1
  row = c(cluster_name, g_levels, summed_freqs, colour)
  spring_cluster_data = rbind(spring_cluster_data, row)
}
colnames(spring_cluster_data) =  c('cluster name', 'guilds in cluster', 'guilds percentiles', 'cluster colour')
write.table(spring_cluster_data, '../../Data/spring_cluster_data.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
