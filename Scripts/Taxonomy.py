import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def read_table(TABLE, used_index, sample):
	lines = open(TABLE, "r").readlines()
	ID, value = [], []

	for line in lines[1:]:
		ID.append(line.split(",")[0])
		value.append(float(line.split(",")[used_index].strip().replace("\n", "")))

	data = [ID, value]

	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/taxonomy.txt', 'a') as table:
		table.write(f'{sample}: {data}\n')

	return data

with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/taxonomy.txt', 'w') as table:
		table.write(f'The data sets are:\n')

H1_13d = read_table("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/metaphlan_feces_phylum_shared.csv", 1, "H1 13d")
H1_20d = read_table("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/metaphlan_feces_phylum_shared.csv", 2, "H1 20d")
H1_27d = read_table("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/metaphlan_feces_phylum_shared.csv", 3, "H1 27d")
H1_35d = read_table("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/metaphlan_feces_phylum_shared.csv", 4, "H1 35d")
H2_10d = read_table("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/metaphlan_feces_phylum_shared.csv", 5, "H2 10d")
H2_17d = read_table("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/metaphlan_feces_phylum_shared.csv", 6, "H2 17d")
H2_24d = read_table("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/metaphlan_feces_phylum_shared.csv", 7, "H2 24d")
H2_31d = read_table("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/metaphlan_feces_phylum_shared.csv", 8, "H2 31d")
H3_9d = read_table("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/metaphlan_feces_phylum_shared.csv", 9, "H3 9d")
H3_16d = read_table("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/metaphlan_feces_phylum_shared.csv", 10, "H3 16d")
H3_23d = read_table("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/metaphlan_feces_phylum_shared.csv", 11, "H3 23d")
H3_30d = read_table("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/metaphlan_feces_phylum_shared.csv", 12, "H3 30d")

H1_13d_val = H1_13d[1]
H1_20d_val = H1_20d[1]
H1_27d_val = H1_27d[1]
H1_35d_val = H1_35d[1]
H2_10d_val = H2_10d[1]
H2_17d_val = H2_17d[1]
H2_24d_val = H2_24d[1]
H2_31d_val = H2_31d[1]
H3_9d_val = H3_9d[1]
H3_16d_val = H3_16d[1]
H3_23d_val = H3_23d[1]
H3_30d_val = H3_30d[1]

def visualization_heatmap(all_labels, tax_groups, data_array):
	fix, axs = plt.subplots(figsize=(10, 6))
	im = axs.imshow(data_array, cmap="YlGn")
	cbra = axs.figure.colorbar(im, ax=axs, aspect=40, shrink=0.8)
	cbra.ax.set_ylabel("Relative Abundance (%)", rotation=-90, va='bottom')

	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/taxonomy.txt', 'a') as table:
		table.write(f'Data array: {data_array}\n')

	axs.set_xticks(range(len(all_labels)))
	axs.set_xticklabels(all_labels, rotation=45, ha='right', rotation_mode="anchor")
	axs.set_yticks(range(len(tax_groups)))
	axs.set_yticklabels(tax_groups)

	for i in range(len(tax_groups)):
		for j in range(len(all_labels)):
			text = axs.text(j, i, f'{data_array[i, j]:.1f}', ha='center')
	axs.set_title("Taxonomy Assignment of Metagenomic Sequences")
	plt.tight_layout()
	plt.savefig(f"/media/ja/data/Mikra/Jane/plasmidome/new/Results/taxonomy_heatmap.svg")

#Data for heatmap
all_labels = ["H1 13d", "H1 20d", "H1 27d", "H1 35d", "H2 10d", "H2 17d", "H2 24d", "H2 31d", "H3 9d", "H3 16d", "H3 23d", "H3 30d"]
taxonomic_groups = ['Actinobacteria', 'Bacteroidota', 'Candidatus melainabacteria', 'Firmicutes', 'Proteobacteria', 'Verrucomicrobia', 'Apicomplexa']
data_array = np.array([H1_13d_val, H1_20d_val, H1_27d_val, H1_35d_val,
	H2_10d_val, H2_17d_val, H2_24d_val, H2_31d_val,
	H3_9d_val, H3_16d_val, H3_23d_val, H3_30d_val])
data_array = data_array.T

visualization_heatmap(all_labels, taxonomic_groups, data_array)