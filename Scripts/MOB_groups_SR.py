import matplotlib
matplotlib.use('Agg')
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm

def kma_res(res_table):
	#read the kma results, specify where the values are positioned in res file
	lines = open(res_table, "r").readlines()
	
	ID, cov, identity, length = [], [], [], []
	for line in lines[1:]:
		ID.append(line.split("\t")[0])
		cov.append(float(line.split("\t")[5]))
		identity.append(float(line.split("\t")[4]))			
		length.append(float(line.split("\t")[3]))

	return np.array(ID), np.array(cov), np.array(identity), np.array(length)
	# return the values in the form of arrays

def kma_mapstat(mapstat_table, sampel):
	#return the ids and correponding readcounts (fragmentCountAln) from mapstat file 
	lines1 = open(mapstat_table, "r").readlines()
	total_length = float(lines1[3].split()[2].strip().replace("#", ""))

	readCount = []
	IDs = []
	for line in lines1[7:]: 
		IDs.append(line.split()[0])
		readCount.append(float(line.split()[14]))

	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/all_groups_MOB_specific_SR.txt', 'a') as data:
		data.write(f'mapstat readcounts for {sampel} are:\n')
		data.write(f'{IDs}\n')
		data.write(f'{readCount}\n')
		data.write(f'{total_length}\n')

	return [IDs, readCount], total_length

def cov80(ids, covs, identity, lengths, sampel):
	#find hits with coverage greater than 80

	covs_80_find = np.where((covs >= 80) & (identity >= 90))
	ids_80 = ids[covs_80_find]
	covs_80 = covs[covs_80_find]
	identity_90 = identity[covs_80_find]
	lengths_80 = lengths[covs_80_find]

	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/all_groups_MOB_specific_SR.txt', 'w') as data:
		data.write(f'The filtered genes for {sampel} are:\n')
		data.write(f'{ids_80}\t{covs_80}\t{identity_90}\t{lengths_80}\n')

	return ids_80, covs_80, identity_90, lengths_80

def RPKM_count(RPKM_ID, RPKM_numreads, RPKM_IDlength, totalNumreads, sampel):
	#count the RPKM
	rpkm_id = []
	RPKMs = []

	for i in range(len(RPKM_ID)):
		id_for_RPKM = RPKM_ID[i]
		pos_numreads = RPKM_numreads[0].index(id_for_RPKM)
		RPKM_readsCount = RPKM_numreads[1][pos_numreads]
		RPKM_Length = RPKM_IDlength[i]

		RPKM = RPKM_readsCount/(RPKM_Length/100 * totalNumreads/1000000)
		rpkm_id.append(id_for_RPKM)
		RPKMs.append(RPKM)

	total_RPKM = sum(RPKMs)

	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/all_groups_MOB_specific_SR.txt', 'a') as data:
		data.write('RKPM values for %s are:\n' % (sampel))
		data.write(f'{rpkm_id}\n')
		data.write(f'{RPKMs}\n')
		data.write(f'{RPKM_readsCount}\n')
		data.write(f'{RPKM_Length}\n')
		data.write(f'{totalNumreads}\n')
		data.write(f'{total_RPKM}\n')

	return [rpkm_id, RPKMs], total_RPKM

def sort_AMR_groups(RPKM_list, cov80_id_list, total, samp):
    # Create mappings from gene IDs to RPKM values and AMR groups
    gene_id_to_rpkm = dict(zip(RPKM_list[0], RPKM_list[1]))
    gene_id_to_group = {gene_id: gene_id for gene_id in cov80_id_list}
    total_RPKM = float(total)
    
    # Initialize dictionary to hold RPKM values for each group
    group_values = {}
    group_sums = {}
    
    # Aggregate RPKM values under their respective AMR groups
    for gene_id, group in gene_id_to_group.items():
        rpkm = gene_id_to_rpkm.get(gene_id)
        if rpkm is not None:
            group_values.setdefault(group, []).append(float(rpkm))
    
    # Calculate counts and average RPKM values for each group
    group_sums = {group: sum(values) for group, values in group_values.items()}
    group_avg = {
        group: sum(values) / total_RPKM if values else 0.0
        for group, values in group_values.items()
    }
    
    # Write results to the specified file
    with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/all_groups_MOB_specific_SR.txt', 'a') as data:
        data.write(f'The AMR groups, counts, and values for {samp} are:\n')
        data.write(f'Groups: {list(group_values.keys())}\n')
        data.write(f'Group sums: {group_sums}\n')
        data.write(f'Values: {group_values}\n')
        data.write(f'The RPKM averages for groups of {samp} in dictionaries are:\n')
        data.write(f'{group_avg}\n')
    
    # Prepare the result list to return
    keys_list = list(group_avg.keys())
    values_list = list(group_avg.values())
    group_list = [keys_list, values_list]

    for group, value in zip(group_list[0], group_list[1]):
    	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/final_info_groups_MOB_SR.txt', 'a') as table:
    		table.write(f'{samp}\t{group}\t{value}\t{group_sums}\t{total}\n')

    return group_list

def top_groups_all(data_sets):
	group_totals = defaultdict(float)

	for data in data_sets:
		for hall_data in data:
			for group, value in zip(hall_data[0], hall_data[1]):
				group_totals[group] += value

	sorted_groups = sorted(group_totals.items(), key=lambda x: x[1], reverse=True)
	top_groups = {group for group, _ in sorted_groups[:20]}

	return top_groups

##Execute the analysis of RPKM
#Read kma .res files and assign position for ID, coverage and length
H1_13d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/02F_mob.res")
H1_20d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/03F_mob.res")
H1_27d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/04F_mob.res")
H1_35d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/05F_mob.res")
H2_10d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/07F_mob.res")
H2_17d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/08F_mob.res")
H2_24d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/09F_mob.res")
H2_31d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/10F_mob.res")
H3_9d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/12F_mob.res")
H3_16d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/13F_mob.res")
H3_23d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/14F_mob.res")
H3_30d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/15F_mob.res")

H1_13_id, H1_13_cov, H1_13_IDidentity, H1_13_IDlength = H1_13d[0], H1_13d[1], H1_13d[2], H1_13d[3]
H1_20_id, H1_20_cov, H1_20_IDidentity, H1_20_IDlength = H1_20d[0], H1_20d[1], H1_20d[2], H1_20d[3]
H1_27_id, H1_27_cov, H1_27_IDidentity, H1_27_IDlength  = H1_27d[0], H1_27d[1], H1_27d[2], H1_27d[3]
H1_35_id, H1_35_cov, H1_35_IDidentity, H1_35_IDlength = H1_35d[0], H1_35d[1], H1_35d[2], H1_35d[3]
H2_10_id, H2_10_cov, H2_10_IDidentity, H2_10_IDlength = H2_10d[0], H2_10d[1], H2_10d[2], H2_10d[3]
H2_17_id, H2_17_cov, H2_17_IDidentity, H2_17_IDlength = H2_17d[0], H2_17d[1], H2_17d[2], H2_17d[3]
H2_24_id, H2_24_cov, H2_24_IDidentity, H2_24_IDlength = H2_24d[0], H2_24d[1], H2_24d[2], H2_24d[3]
H2_31_id, H2_31_cov, H2_31_IDidentity, H2_31_IDlength = H2_31d[0], H2_31d[1], H2_31d[2], H2_31d[3]
H3_9_id, H3_9_cov, H3_9_IDidentity, H3_9_IDlength = H3_9d[0], H3_9d[1], H3_9d[2], H3_9d[3]
H3_16_id, H3_16_cov, H3_16_IDidentity, H3_16_IDlength = H3_16d[0], H3_16d[1], H3_16d[2], H3_16d[3]
H3_23_id, H3_23_cov, H3_23_IDidentity, H3_23_IDlength = H3_23d[0], H3_23d[1], H3_23d[2], H3_23d[3]
H3_30_id, H3_30_cov, H3_30_IDidentity, H3_30_IDlength = H3_30d[0], H3_30d[1], H3_30d[2], H3_30d[3]

#Filter out the data with coverage =< 80
H1_13_id_80, H1_13_cov_80, H1_13_IDidentity_90, H1_13_IDlength_80 = cov80(H1_13_id, H1_13_cov, H1_13_IDidentity, H1_13_IDlength, 'H1_13')
H1_20_id_80, H1_20_cov_80, H1_20_IDidentity_90, H1_20_IDlength_80 = cov80(H1_20_id, H1_20_cov, H1_20_IDidentity, H1_20_IDlength, 'H1_20')
H1_27_id_80, H1_27_cov_80, H1_27_IDidentity_90, H1_27_IDlength_80 = cov80(H1_27_id, H1_27_cov, H1_27_IDidentity, H1_27_IDlength, 'H1_27')
H1_35_id_80, H1_35_cov_80, H1_35_IDidentity_90, H1_35_IDlength_80 = cov80(H1_35_id, H1_35_cov, H1_35_IDidentity, H1_35_IDlength, 'H1_35')
H2_10_id_80, H2_10_cov_80, H2_10_IDidentity_90, H2_10_IDlength_80 = cov80(H2_10_id, H2_10_cov, H2_10_IDidentity, H2_10_IDlength, 'H2_10')
H2_17_id_80, H2_17_cov_80, H2_17_IDidentity_90, H2_17_IDlength_80 = cov80(H2_17_id, H2_17_cov, H2_17_IDidentity, H2_17_IDlength, 'H2_17')
H2_24_id_80, H2_24_cov_80, H2_24_IDidentity_90, H2_24_IDlength_80 = cov80(H2_24_id, H2_24_cov, H2_24_IDidentity, H2_24_IDlength, 'H2_24')
H2_31_id_80, H2_31_cov_80, H2_31_IDidentity_90, H2_31_IDlength_80 = cov80(H2_31_id, H2_31_cov, H2_31_IDidentity, H2_31_IDlength, 'H2_31')
H3_9_id_80, H3_9_cov_80, H3_9_IDidentity_90, H3_9_IDlength_80 = cov80(H3_9_id, H3_9_cov, H3_9_IDidentity, H3_9_IDlength, 'H3_9')
H3_16_id_80, H3_16_cov_80, H3_16_IDidentity_90, H3_16_IDlength_80 = cov80(H3_16_id, H3_16_cov, H3_16_IDidentity, H3_16_IDlength, 'H3_16')
H3_23_id_80, H3_23_cov_80, H3_23_IDidentity_90, H3_23_IDlength_80 = cov80(H3_23_id, H3_23_cov, H3_23_IDidentity, H3_23_IDlength, 'H3_23')
H3_30_id_80, H3_30_cov_80, H3_30_IDidentity_90, H3_30_IDlength_80 = cov80(H3_30_id, H3_30_cov, H3_30_IDidentity, H3_30_IDlength, 'H3_30')

#Read kma .mapstat files and assing position of ID + read count, and total length
H1_13_numread_pair, H1_13_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/02F_mob.mapstat", 'H1_13')
H1_20_numread_pair, H1_20_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/03F_mob.mapstat", 'H1_20')
H1_27_numread_pair, H1_27_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/04F_mob.mapstat", 'H1_27')
H1_35_numread_pair, H1_35_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/05F_mob.mapstat", 'H1_35')

H2_10_numread_pair, H2_10_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/07F_mob.mapstat", 'H2_10')
H2_17_numread_pair, H2_17_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/08F_mob.mapstat", 'H2_17')
H2_24_numread_pair, H2_24_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/09F_mob.mapstat", 'H2_24')
H2_31_numread_pair, H2_31_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/10F_mob.mapstat", 'H2_31')


H3_9_numread_pair, H3_9_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/12F_mob.mapstat", 'H3_9')
H3_16_numread_pair, H3_16_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/13F_mob.mapstat", 'H3_16')
H3_23_numread_pair, H3_23_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/14F_mob.mapstat", 'H3_23')
H3_30_numread_pair, H3_30_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/15F_mob.mapstat", 'H3_30')

#Count the RPKMs
H1_13_RPKM, H1_13_tot_RPKM = RPKM_count(H1_13_id_80, H1_13_numread_pair, H1_13_IDlength_80, H1_13_totalLength, 'H1_13')
H1_20_RPKM, H1_20_tot_RPKM = RPKM_count(H1_20_id_80, H1_20_numread_pair, H1_20_IDlength_80, H1_20_totalLength, 'H1_20')
H1_27_RPKM, H1_27_tot_RPKM = RPKM_count(H1_27_id_80, H1_27_numread_pair, H1_27_IDlength_80, H1_27_totalLength, 'H1_27')
H1_35_RPKM, H1_35_tot_RPKM = RPKM_count(H1_35_id_80, H1_35_numread_pair, H1_35_IDlength_80, H1_35_totalLength, 'H1_35')

H2_10_RPKM, H2_10_tot_RPKM = RPKM_count(H2_10_id_80, H2_10_numread_pair, H2_10_IDlength_80, H2_10_totalLength, 'H2_10')
H2_17_RPKM, H2_17_tot_RPKM = RPKM_count(H2_17_id_80, H2_17_numread_pair, H2_17_IDlength_80, H2_17_totalLength, 'H2_17')
H2_24_RPKM, H2_24_tot_RPKM = RPKM_count(H2_24_id_80, H2_24_numread_pair, H2_24_IDlength_80, H2_24_totalLength, 'H2_24')
H2_31_RPKM, H2_31_tot_RPKM = RPKM_count(H2_31_id_80, H2_31_numread_pair, H2_31_IDlength_80, H2_31_totalLength, 'H2_31')

H3_9_RPKM, H3_9_tot_RPKM = RPKM_count(H3_9_id_80, H3_9_numread_pair, H3_9_IDlength_80, H3_9_totalLength, 'H3_9')
H3_16_RPKM, H3_16_tot_RPKM = RPKM_count(H3_16_id_80, H3_16_numread_pair, H3_16_IDlength_80, H3_16_totalLength, 'H3_16')
H3_23_RPKM, H3_23_tot_RPKM = RPKM_count(H3_23_id_80, H3_23_numread_pair, H3_23_IDlength_80, H3_23_totalLength, 'H3_23')
H3_30_RPKM, H3_30_tot_RPKM = RPKM_count(H3_30_id_80, H3_30_numread_pair, H3_30_IDlength_80, H3_30_totalLength, 'H3_30')

#Put together all RPKMs and sort them alphabetically
with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/final_info_groups_MOB_SR.txt', 'w') as table:
	table.write(f'Sample ID\tATB group\tGroup AVG (%)\tRPKM sum\ttotal RPKM\n')

H1_13_groups = sort_AMR_groups(H1_13_RPKM, H1_13_id_80, H1_13_tot_RPKM, "H1 13")
H1_20_groups = sort_AMR_groups(H1_20_RPKM, H1_20_id_80, H1_20_tot_RPKM, "H1 20")
H1_27_groups = sort_AMR_groups(H1_27_RPKM, H1_27_id_80, H1_27_tot_RPKM, "H1 27")
H1_35_groups = sort_AMR_groups(H1_35_RPKM, H1_35_id_80, H1_35_tot_RPKM, "H1 35")

H2_10_groups = sort_AMR_groups(H2_10_RPKM, H2_10_id_80, H2_10_tot_RPKM, "H2 10")
H2_17_groups = sort_AMR_groups(H2_17_RPKM, H2_17_id_80, H2_17_tot_RPKM, "H2 17")
H2_24_groups = sort_AMR_groups(H2_24_RPKM, H2_24_id_80, H2_24_tot_RPKM, "H2 24")
H2_31_groups = sort_AMR_groups(H2_31_RPKM, H2_31_id_80, H2_31_tot_RPKM, "H2 31")

H3_9_groups = sort_AMR_groups(H3_9_RPKM, H3_9_id_80, H3_9_tot_RPKM, "H3 9")
H3_16_groups = sort_AMR_groups(H3_16_RPKM, H3_16_id_80, H3_16_tot_RPKM, "H3 16")
H3_23_groups = sort_AMR_groups(H3_23_RPKM, H3_23_id_80, H3_23_tot_RPKM, "H3 23")
H3_30_groups = sort_AMR_groups(H3_30_RPKM, H3_30_id_80, H3_30_tot_RPKM, "H3 30")

color_mapping = {
	'NZ_CP011417|MOBP': '#ffe9f0',
	'NZ_CP013027|MOBF': '#A6C2C0',
	'JN985534|MOB_unknown': '#eec5b7',
	'NC_018997|MOBP': '#ea9999',
	'CP018359|MOB_unknown': '#FFE69F',
	'NC_019078|MOBP': '#F0E188',
	'HM021326|MOBP': '#ffe3c0',
	'CP020736|MOBP': '#FFABAB',
	'MF083142|MOBP': '#d2d7e7',
	'NC_019136|MOBP': '#9edae5',
	'NC_011990|MOBP': '#C78BA5',
	'NZ_CP010316|MOBP': '#cfe2f3',
	'CP022165|MOBP': '#d9ead3',
	'NZ_CP009579|MOBV': '#b4a7d6',
	'NZ_CP008845|MOB_unknown': '#dda298',
	'AB027308|MOBP': '#84D6B3',
	'CP019195|MOBP': '#B8DECC',
	'NC_019986|MOB_unknown': '#d5a6bd',
	'KY014464|MOBP': '#dfbac9',
	'NC_011754|MOBP': '#BCBD22'
}

info_names = ["H1", "H2", "H3"]
time_labels = [
    ["13d", "20d", "27d", "35d"],  # House 1
    ["10d", "17d", "24d", "31d"],  # House 2
    ["9d", "16d", "23d", "30d"]    # House 3
]
data_sets = [
    [H1_13_groups, H2_10_groups, H3_9_groups],
    [H1_20_groups, H2_17_groups, H3_16_groups],
    [H1_27_groups, H2_24_groups, H3_23_groups],
    [H1_35_groups, H2_31_groups, H3_30_groups],
]

def visualisation_combined_bar(info_names, data_sets, time_labels):
	top_groups = top_groups_all(data_sets)
	plt.figure(figsize=(14, 10))
	n_halls = len(info_names)
	n_points = len(data_sets)
	width = 0.7
	group_spacing = 0.1

	x = []  # Position of each bar
	xtick_labels = []

	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/all_groups_MOB_specific_SR.txt', 'a') as data:
		data.write('The color mapping groups are:\n')
		data.write(f'{color_mapping}\n')

	relative_abundance_storage = {key: [] for key in color_mapping.keys()}
	relative_abundance_storage['other'] = []

	# Track which AMR groups have already been added to the legend
	legend_labels = set()

	# Prepare to plot each time point
	for hall_idx, hall_name in enumerate(info_names):
		hall_start = hall_idx * (n_points + group_spacing)

		for time_idx, time_label in enumerate(time_labels[hall_idx]):
			x.append(hall_start + time_idx * (width + 0.05))
			xtick_labels.append(f'{hall_name}_{time_label}')

	for idx, (pos, label) in enumerate(zip(x, xtick_labels)):
		hall_idx = idx // n_points
		time_idx = idx % n_points
		data = data_sets[time_idx][hall_idx]

		all_data = {key: 0 for key in top_groups | {'other'}}  # Initialize with zero for each AMR group
		all_data['other'] = 0

		for info_id, value in zip(data[0], data[1]):
			if info_id in color_mapping:
				all_data[info_id] += value
			else:
				all_data['other'] += value

		normalized_data = np.array([all_data[key] for key in all_data]) * 100
		bottom = np.zeros(1)

		for info_id, height in zip(all_data.keys(), normalized_data):
			color = color_mapping.get(info_id, '#7f7f7f')

			# Přidání bílého ohraničení kolem hodnot
			plt.bar(pos, height, width, color=color, bottom=bottom, edgecolor='white', linewidth=0.5)

			# Přidání popisku do legendy pouze pokud tam ještě není
			if info_id not in legend_labels:
				plt.bar(pos, height, width, label=f'{info_id}', color=color, bottom=bottom, edgecolor='white', linewidth=0.5)
				legend_labels.add(info_id)
			bottom += height

			relative_abundance_storage[info_id].append(height)

	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/all_groups_MOB_specific_SR.txt', 'a') as data:
		data.write('The relative abundances  in "%" are:\n')
		data.write(f'{relative_abundance_storage}\n')

	plt.xlabel("Houses and Collection Points")
	plt.ylabel("Relative Abundance (%)")
	plt.title("Relative abundances of MOB groups in metagenome")
	plt.xticks(x, xtick_labels, rotation=45, ha='right')
	plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	plt.tight_layout()

	plt.savefig(f"/media/ja/data/Mikra/Jane/plasmidome/new/Results/combined_plot_MOB_SR_final.svg")

visualisation_combined_bar(info_names, data_sets, time_labels)