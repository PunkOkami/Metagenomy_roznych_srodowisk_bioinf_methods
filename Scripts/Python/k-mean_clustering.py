import csv
import numpy as np
from skbio.diversity import alpha_diversity
from random import sample
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.cluster import KMeans
from scipy.spatial import ConvexHull

# Reads data from two files, one to sample from (that is more random set) to create bootsrap data and real data
freq_data_file = open('species_freq_table_filtered.tsv')
sampling_data_file = open('species_freq_table_expanded.tsv')
freq_data_reader = csv.reader(freq_data_file, delimiter='\t')
sampling_data_reader = csv.reader(sampling_data_file, delimiter='\t')

# Creates matrix of real data to calculate indexes used as clustering parameters
sample_names = []
first_line = True
real_data_matrix = []
samples_indexes = slice(0)
sample_num = 0
species_count = 0
real_data_list = []
for line in freq_data_reader:
	if first_line:
		samples_indexes = slice(line.index('Species')+1, line.index('Kingdom'))
		sample_names = line[samples_indexes]
		sample_num = samples_indexes.stop - samples_indexes.start
		real_data_matrix = [[] for i in range(sample_num)]
		first_line = False
		continue
	line = line[samples_indexes]
	for i in range(sample_num):
		real_data_matrix[i].append(float(line[i]))
	species_count += 1
	for num in line:
		real_data_list.append(float(num))

# Reads data in sampling set, but as a list as there is no need to keep samples separate
first_line = True
sample_data_list = []
samples_indexes = slice(0)
sample_num = 0
for line in sampling_data_reader:
	if first_line:
		samples_indexes = slice(line.index('Species')+1, line.index('Kingdom'))
		first_line = False
		continue
	line = line[samples_indexes]
	for num in line:
		sample_data_list.append(float(num))

# Calculates Margalef index for real data
margalef_of_real_data = list(alpha_diversity('margalef', real_data_matrix))

# Creates 300 random samples and calculates Margalef of that data
rep_num = 300
booted_margalefs = []
for i in range(rep_num):
	random_sample = sample(sample_data_list, species_count)
	random_margalef = alpha_diversity('margalef', random_sample)[0]
	booted_margalefs.append(random_margalef)

# Calculates max value in each sample to later get max value in whole matrix
maxes_in_sample_data = [int(max(sample)) for sample in real_data_matrix]
simpson_of_real_data = list(alpha_diversity('simpson_e', real_data_matrix))

# Here bootstrap method differs slightly from one before as this index is more sensive to randmness of data less to richness of the sample
booted_simpsons = []
for i in range(rep_num):
	random_sample = sample(range(0, int(max(sample_data_list))), species_count)
	random_simpson = alpha_diversity('simpson_e', random_sample)[0]
	booted_simpsons.append(random_simpson)

# Creates simple scatterplot with all points, but no clusters to see how much they differ without that clue
all_margalef = margalef_of_real_data.copy()
all_margalef.extend(booted_margalefs)
all_simpson = simpson_of_real_data.copy()
all_simpson.extend(booted_simpsons)
plt.scatter(all_margalef, all_simpson)
plt.title('Alpha diversity')
plt.xlabel('Margalef index')
plt.ylabel('Simpson evenness index')
plt.show()


all_points = [tup for tup in zip(all_margalef, all_simpson)]
sample_points = np.array([tup for tup in zip(margalef_of_real_data, simpson_of_real_data)])
booted_points = np.array([tup for tup in zip(booted_margalefs, booted_simpsons)])
model = KMeans(n_clusters=3)
model.fit(all_points)
labels = model.predict(all_points)
plt.scatter(all_margalef, all_simpson, c=np.array(['green', 'blue', 'red', 'brown'])[np.array(labels)])
plt.scatter(model.cluster_centers_[:, 0], model.cluster_centers_[:, 1], marker='D', s=200, alpha=0.75)
forest_points = []
pots_points = []
for i in range(len(sample_names)):
	print(f'{sample_names[i]} --- {labels[i]}')
	if labels[i] == 1:
		forest_points.append(list(sample_points[i]))
	elif labels[i] == 2:
		pots_points.append(list(sample_points[i]))

# Creates first clusters of points then convex hulls that show borders of the clusters and plots that
forest_points = np.array(forest_points)
pots_points = np.array(pots_points)
forest_hull = ConvexHull(forest_points)
booted_hull = ConvexHull(booted_points)
pots_hull = ConvexHull(pots_points)
forest_hull_x = np.append(forest_points[forest_hull.vertices, 0], forest_points[forest_hull.vertices, 0][0])
forest_hull_y = np.append(forest_points[forest_hull.vertices, 1], forest_points[forest_hull.vertices, 1][0])
pots_hull_x = np.append(pots_points[pots_hull.vertices, 0], pots_points[pots_hull.vertices, 0][0])
pots_hull_y = np.append(pots_points[pots_hull.vertices, 1], pots_points[pots_hull.vertices, 1][0])
booted_hull_x = np.append(booted_points[booted_hull.vertices, 0], booted_points[booted_hull.vertices, 0][0])
booted_hull_y = np.append(booted_points[booted_hull.vertices, 1], booted_points[booted_hull.vertices, 1][0])
plt.fill(booted_hull_x, booted_hull_y, alpha=0.3, c=np.array([0, 0.501, 0]))
plt.fill(forest_hull_x, forest_hull_y, alpha=0.3, c=np.array([0, 0, 1]))
plt.fill(pots_hull_x, pots_hull_y, alpha=0.3, c=np.array([1, 0, 0]))
legend_elements = [Line2D([0], [0], marker='o', color='w', label='Forest samples', markerfacecolor='blue', markersize=5),
				Line2D([0], [0], marker='o', color='w', label='Pots samples', markerfacecolor='red', markersize=5),
				Line2D([0], [0], marker='o', color='w', label='Bootstrap samples', markerfacecolor='green', markersize=5)]
plt.legend(handles=legend_elements, loc='upper right')
plt.title('Alpha diversity')
plt.xlabel('Margalef index')
plt.ylabel('Simpson evenness index')
plt.show()

# This clustering is to check is there are two clusters that correspond to pots and forest samples
model = KMeans(n_clusters=2)
model.fit(sample_points)
labels = model.predict(sample_points)
plt.scatter(margalef_of_real_data, simpson_of_real_data, c=np.array(['brown', 'black', 'red'])[np.array(labels)])
plt.scatter(model.cluster_centers_[:, 0], model.cluster_centers_[:, 1], marker='D', s=200, alpha=0.75)
plt.title('Alpha diversity of real data\n')
plt.xlabel('Margalef index')
plt.ylabel('Simpson evenness index')
plt.show()
