import csv
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from PIL import ImageColor
from copy import deepcopy
import math

# Setting several guild order in table as well as colours for each guild metacategory
guilds = ['Bcn', 'ENA', 'End', 'Myc', 'Nop', 'Pat', 'Sap']
colours = ['blue', 'black', 'red', 'green', 'yellow', 'purple', 'brown']
colours = [mcolors.CSS4_COLORS[colour] for colour in colours]
colours = [ImageColor.getcolor(colour, 'RGBA') for colour in colours]
colours = [list(tup) for tup in colours]

# The code is reading from two files produced using two network construction methods: SPRING and SpiecEasi. Both of them show to be good for data we were working with as modularity was stable and usually over 0.25 and for smaller graphs could reach over 0.5 even. SpiecEasi is usually better, always over 0.36 modularity
cluster_data_files = Path('../../Data').glob('*cluster_data.tsv')
for path in cluster_data_files:

	# Opens file and reads data
	file = open(path)
	method = path.stem.split('_')[0]
	reader = csv.reader(file, delimiter='\t')
	all_clusters_data = {}
	cluster_names = []
	cluster_colours = []
	guilds_index = 0
	percentiles_index = 0
	colour_index = 0
	first_line = True
	for line in reader:
		cluster_data = {guild: 0.0 for guild in guilds}
		if first_line:
			guilds_index = line.index('guilds in cluster')
			percentiles_index = line.index('guilds percentiles')
			colour_index = line.index('cluster colour')
			first_line = False
			continue

		# Writes guild proportions from a row into dictionary and also reads and saves what colour was assinged to each cluster. This makes sure that colours used on network plot and in table correspond
		guilds_col = sorted(line[guilds_index].split('_'))
		percentiles_col = [float(val) for val in line[percentiles_index].split('_')]
		for i, guild in enumerate(guilds_col):
			cluster_data[guild] = percentiles_col[i]
		all_clusters_data[line[0]] = cluster_data
		c_colour = line[colour_index].split('_')
		c_colour = list(list(set([ImageColor.getcolor(col, 'RGBA') for col in c_colour]))[0])
		c_colour = [val/255 for val in c_colour]
		cluster_names.append(line.pop(0))
		cluster_colours.append(c_colour)

	# Parses data from a dict to create matrix of values and matrix of colour values that are input for table constructor
	table_data = []
	colours_table = []
	for cluster, cluster_data in all_clusters_data.items():
		num_row = list(cluster_data.values())
		colours_row = deepcopy(colours)
		row = []
		for i, guild_perc in enumerate(num_row):
			colour = colours_row[i]
			colour[3] = round(colour[3]*guild_perc)
			colour = [val/255 for val in colour]
			row.append(round(colour[3], 2))
			colour[3] = math.ceil(colour[3]*10)/10
			colours_row[i] = colour
		colours_table.append(colours_row)
		table_data.append(row)

	# Creates plot, formats it and saves
	ccolours = ['pink' for i in range(len(guilds))]
	fig, ax = plt.subplots(figsize=(6, 5))
	the_table = ax.table(cellText=table_data, rowLabels=cluster_names, rowLoc='center', colLabels=guilds, loc='center',
						colColours=ccolours, rowColours=cluster_colours, colWidths=[0.07 for i in range(len(guilds))],
						cellColours=colours_table, cellLoc='center')
	the_table.scale(1.2, 2)
	the_table.auto_set_font_size(False)
	the_table.set_fontsize(10)
	# plt.legend(handles=handles, loc='right')
	ax = plt.gca()
	ax.get_xaxis().set_visible(False)
	ax.get_yaxis().set_visible(False)
	plt.box(on=None)
	plt.title(f'Clusters by {method}', loc='center', fontdict={'size': 20})
	plt.draw()
	plt.savefig(f'../../Graphs/NetCoMi/{method}_clusters_makeup_table.png', dpi=150)
