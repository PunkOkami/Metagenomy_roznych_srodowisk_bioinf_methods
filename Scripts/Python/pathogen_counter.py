import csv

# Reads data and pulls out taxa in metacategory Pathogens
species_table = open('../../Data/species_freq_table_filtered.tsv')
species_reader = csv.reader(species_table, delimiter='\t')

first_line = True
plant_pathogen_dict = {}
pots = []
forest = []
for row in species_reader:
	if first_line:
		sample_slice = slice(row.index('Species')+1, row.index('Kingdom'))
		pots = [row.index(column) for column in row[sample_slice] if 'z' in column or 'c' in column]
		forest = [row.index(column) for column in row[sample_slice] if row.index(column) not in pots]
		first_line = False
		continue
	if 'Pat' in row[-1]:
		plant_pathogen_dict[row[0]] = (sum([float(row[index]) for index in pots]), sum([float(row[index]) for index in forest]))

# Creates output file, sorts Pats by abundance in pots and writes a list of ten species with highest to output file
out_file = open('pathogen_list.txt', mode='w')
out_file.write(f'Number of pathogen species: {len(plant_pathogen_dict.items())}\n\n')
plant_pathogen_dict = sorted(plant_pathogen_dict.items(), key=lambda x: x[1], reverse=True)
out_file.write('10 species that have highest abundance in pots:\n')
for i, pathogen in enumerate(plant_pathogen_dict):
	if i > 9:
		break
	out_file.write(f' - {pathogen[0]}\n')
out_file.write('\n')

# Writes data about all taxa that it found to a file, split into group of pots and forest
out_file.write('Species - abundance in pots samples - abundance in forest samples\n')
for (pathogen, freq) in plant_pathogen_dict:
	out_file.write(f'{pathogen} - {freq[0]} - {freq[1]}\n')
