import seaborn
import matplotlib.pyplot as plt
import csv
import copy
seaborn.set(font_scale=0.65)

# Open file and reads data of species that are in metacategory specified by guild_wanted
data_file = open('../../Data/species_freq_table_filtered.tsv')
data_reader = csv.reader(data_file, delimiter='\t')

guild_wanted = 'Myc'
most_common_dict = {}
freq_data = {}
sample_names = []
sample_indexes = slice(0, 0)
species_index = 0
guild_index = 0
first_line = True
for row in data_reader:
	if first_line:
		species_index = row.index('Species')
		guild_index = row.index('Guild category')
		sample_indexes = slice(species_index+1, row.index('Kingdom'))
		sample_names = row[sample_indexes]
		first_line = False
		continue
	if row[guild_index] == guild_wanted:
		nums = [float(num) for num in row[sample_indexes]]
		freq_data[row[species_index]] = nums
		

# Changes the data format to match what seaborn expect and converst absolute values to percentile as this method shows to be best suited
heatmap_data = []
species_names = []
for key, value in freq_data.items():
	copyied = copy.deepcopy(value)
	heatmap_data.append(copyied)
	species_names.append(key)
heatmap_data_flattened = []
for row in heatmap_data:
	heatmap_data_flattened.extend(row)
heatmap_data_flattened.sort()
percentile_dict = {num: round((heatmap_data_flattened.index(num) + 1) / len(heatmap_data_flattened) * 100, 2) for num in heatmap_data_flattened}
for i in range(len(heatmap_data)):
	seq = heatmap_data[i]
	for j in range(len(seq)):
		num = seq[j]
		perc = percentile_dict[num]
		seq[j] = perc
	heatmap_data[i] = seq

# Plots clustermap and saves to file
species_guilds_heatmap = seaborn.clustermap(heatmap_data, cmap="crest", yticklabels=False, xticklabels=sample_names,
											col_cluster=True, method='ward')
plt.savefig(f'../../Graphs/Clustering/Pairwise_clustering/{guild_wanted}_clustermap.png', dpi=150)

