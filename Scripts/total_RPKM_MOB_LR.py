import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
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
	lines = open(mapstat_table, "r").readlines()
	total_length = float(lines[3].split()[2].strip().replace("#", ""))

	readCount = []
	IDs = []
	for line in lines[7:]: 
		IDs.append(line.split()[0])
		readCount.append(float(line.split()[14]))

	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/all_totalRPKM_MOB.txt', 'a') as data:
		data.write(f'mapstat readcounts for {sampel} are:\n')
		data.write(f'{IDs}\n')
		data.write(f'{readCount}\n')
		data.write(f'{total_length}\n')

	return [IDs, readCount], total_length

def cov80(ids, covs, identity, lengths):
	#find hits with coverage greater than 80

	covs_80_find = np.where((covs >= 80) & (identity >= 90))
	ids_80 = ids[covs_80_find]
	covs_80 = covs[covs_80_find]
	identity_90 = identity[covs_80_find]
	lengths_80 = lengths[covs_80_find]

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

	total_rpkm = sum(RPKMs)

	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/all_totalRPKM_MOB.txt', 'a') as data:
		data.write('RKPM values for %s are:\n' % (sampel))
		data.write(f'{rpkm_id}\n')
		data.write(f'{RPKMs}\n')
		data.write(f'{RPKM_readsCount}\n')
		data.write(f'{RPKM_Length}\n')
		data.write(f'{totalNumreads}\n')
		data.write(f'{total_rpkm}\n')

	return total_rpkm
		
##Execute the analysis of RPKM
#Read kma .res files and assign position for ID, coverage and length
H1_13d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H1_13d_o200_mob.res")
H1_20d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H1_20d_o200_mob.res")
H1_27d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H1_27d_o200_mob.res")
H1_35d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H1_35d_o200_mob.res")
H2_10d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H2_10d_o200_mob.res")
H2_17d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H2_17d_o200_mob.res")
H2_24d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H2_24d_o200_mob.res")
H2_31d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H2_31d_o200_mob.res")
H3_9d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H3_9d_o200_mob.res")
H3_16d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H3_16d_o200_mob.res")
H3_23d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H3_23d_o200_mob.res")
H3_30d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H3_30d_o200_mob.res")

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
H1_13_id_80, H1_13_cov_80, H1_13_IDidentity_90, H1_13_IDlength_80 = cov80(H1_13_id, H1_13_cov, H1_13_IDidentity, H1_13_IDlength)
H1_20_id_80, H1_20_cov_80, H1_20_IDidentity_90, H1_20_IDlength_80 = cov80(H1_20_id, H1_20_cov, H1_20_IDidentity, H1_20_IDlength)
H1_27_id_80, H1_27_cov_80, H1_27_IDidentity_90, H1_27_IDlength_80 = cov80(H1_27_id, H1_27_cov, H1_27_IDidentity, H1_27_IDlength)
H1_35_id_80, H1_35_cov_80, H1_35_IDidentity_90, H1_35_IDlength_80 = cov80(H1_35_id, H1_35_cov, H1_35_IDidentity, H1_35_IDlength)

H2_10_id_80, H2_10_cov_80, H2_10_IDidentity_90, H2_10_IDlength_80 = cov80(H2_10_id, H2_10_cov, H2_10_IDidentity, H2_10_IDlength)
H2_17_id_80, H2_17_cov_80, H2_17_IDidentity_90, H2_17_IDlength_80 = cov80(H2_17_id, H2_17_cov, H2_17_IDidentity, H2_17_IDlength)
H2_24_id_80, H2_24_cov_80, H2_24_IDidentity_90, H2_24_IDlength_80 = cov80(H2_24_id, H2_24_cov, H2_24_IDidentity, H2_24_IDlength)
H2_31_id_80, H2_31_cov_80, H2_31_IDidentity_90, H2_31_IDlength_80 = cov80(H2_31_id, H2_31_cov, H2_31_IDidentity, H2_31_IDlength)

H3_9_id_80, H3_9_cov_80, H3_9_IDidentity_90, H3_9_IDlength_80 = cov80(H3_9_id, H3_9_cov, H3_9_IDidentity, H3_9_IDlength)
H3_16_id_80, H3_16_cov_80, H3_16_IDidentity_90, H3_16_IDlength_80 = cov80(H3_16_id, H3_16_cov, H3_16_IDidentity, H3_16_IDlength)
H3_23_id_80, H3_23_cov_80, H3_23_IDidentity_90, H3_23_IDlength_80 = cov80(H3_23_id, H3_23_cov, H3_23_IDidentity, H3_23_IDlength)
H3_30_id_80, H3_30_cov_80, H3_30_IDidentity_90, H3_30_IDlength_80 = cov80(H3_30_id, H3_30_cov, H3_30_IDidentity, H3_30_IDlength)


with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/all_totalRPKM_MOB.txt', 'w') as data:
	data.write('The filtered (cov80, id90) genes are:\n')
	data.write('H1_13d\n')
	data.write(f'{H1_13_id_80}\t{H1_13_cov_80}\t{H1_13_IDidentity}\t{H1_13_IDlength_80}\n')
	data.write('H1_20d\n')
	data.write(f'{H1_20_id_80}\t{H1_20_cov_80}\t{H1_20_IDidentity}\t{H1_20_IDlength_80}\n')
	data.write('H1_27d\n')
	data.write(f'{H1_27_id_80}\t{H1_27_cov_80}\t{H1_27_IDidentity}\t{H1_27_IDlength_80}\n')
	data.write('H1_35d\n')
	data.write(f'{H1_35_id_80}\t{H1_35_cov_80}\t{H1_35_IDidentity}\t{H1_35_IDlength_80}\n')

	data.write('H2_10d\n')
	data.write(f'{H2_10_id_80}\t{H2_10_cov_80}\t{H2_10_IDidentity}\t{H2_10_IDlength_80}\n')
	data.write('H2_17d\n')
	data.write(f'{H2_17_id_80}\t{H2_17_cov_80}\t{H2_17_IDidentity}\t{H2_17_IDlength_80}\n')
	data.write('H2_24d\n')
	data.write(f'{H2_24_id_80}\t{H2_24_cov_80}\t{H2_24_IDidentity}\t{H2_24_IDlength_80}\n')
	data.write('H2_31d\n')
	data.write(f'{H2_31_id_80}\t{H2_31_cov_80}\t{H2_31_IDidentity}\t{H2_31_IDlength_80}\n')

	data.write('H3_9d\n')
	data.write(f'{H3_9_id_80}\t{H3_9_cov_80}\t{H3_9_IDidentity}\t{H3_9_IDlength_80}\n')
	data.write('H3_16d\n')
	data.write(f'{H3_16_id_80}\t{H3_16_cov_80}\t{H3_16_IDidentity}\t{H3_16_IDlength_80}\n')
	data.write('H3_23d\n')
	data.write(f'{H3_23_id_80}\t{H3_23_cov_80}\t{H3_23_IDidentity}\t{H3_23_IDlength_80}\n')
	data.write('H3_30d\n')
	data.write(f'{H3_30_id_80}\t{H3_30_cov_80}\t{H3_30_IDidentity}\t{H3_30_IDlength_80}\n')

#Read kma .mapstat files and assing position of ID + read count, and total length
H1_13_numread_pair, H1_13_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H1_13d_o200_mob.mapstat", 'H1_13')
H1_20_numread_pair, H1_20_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H1_20d_o200_mob.mapstat", 'H1_20')
H1_27_numread_pair, H1_27_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H1_27d_o200_mob.mapstat", 'H1_27')
H1_35_numread_pair, H1_35_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H1_35d_o200_mob.mapstat", 'H1_35')

H2_10_numread_pair, H2_10_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H2_10d_o200_mob.mapstat", 'H2_10')
H2_17_numread_pair, H2_17_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H2_17d_o200_mob.mapstat", 'H2_17')
H2_24_numread_pair, H2_24_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H2_24d_o200_mob.mapstat", 'H2_24')
H2_31_numread_pair, H2_31_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H2_31d_o200_mob.mapstat", 'H2_31')


H3_9_numread_pair, H3_9_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H3_9d_o200_mob.mapstat", 'H3_9')
H3_16_numread_pair, H3_16_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H3_16d_o200_mob.mapstat", 'H3_16')
H3_23_numread_pair, H3_23_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H3_23d_o200_mob.mapstat", 'H3_23')
H3_30_numread_pair, H3_30_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/mob/H3_30d_o200_mob.mapstat", 'H3_30')

#Count the RPKMs
H1_13_RPKM = RPKM_count(H1_13_id_80, H1_13_numread_pair, H1_13_IDlength_80, H1_13_totalLength, 'H1_13')
H1_20_RPKM = RPKM_count(H1_20_id_80, H1_20_numread_pair, H1_20_IDlength_80, H1_20_totalLength, 'H1_20')
H1_27_RPKM = RPKM_count(H1_27_id_80, H1_27_numread_pair, H1_27_IDlength_80, H1_27_totalLength, 'H1_27')
H1_35_RPKM = RPKM_count(H1_35_id_80, H1_35_numread_pair, H1_35_IDlength_80, H1_35_totalLength, 'H1_35')

H2_10_RPKM = RPKM_count(H2_10_id_80, H2_10_numread_pair, H2_10_IDlength_80, H2_10_totalLength, 'H2_10')
H2_17_RPKM = RPKM_count(H2_17_id_80, H2_17_numread_pair, H2_17_IDlength_80, H2_17_totalLength, 'H2_17')
H2_24_RPKM = RPKM_count(H2_24_id_80, H2_24_numread_pair, H2_24_IDlength_80, H2_24_totalLength, 'H2_24')
H2_31_RPKM = RPKM_count(H2_31_id_80, H2_31_numread_pair, H2_31_IDlength_80, H2_31_totalLength, 'H2_31')

H3_9_RPKM = RPKM_count(H3_9_id_80, H3_9_numread_pair, H3_9_IDlength_80, H3_9_totalLength, 'H3_9')
H3_16_RPKM = RPKM_count(H3_16_id_80, H3_16_numread_pair, H3_16_IDlength_80, H3_16_totalLength, 'H3_16')
H3_23_RPKM = RPKM_count(H3_23_id_80, H3_23_numread_pair, H3_23_IDlength_80, H3_23_totalLength, 'H3_23')
H3_30_RPKM = RPKM_count(H3_30_id_80, H3_30_numread_pair, H3_30_IDlength_80, H3_30_totalLength, 'H3_30')

#Create dictionaries for each hall
def visualise_bar(house_data, house_labels):
	house_colors = ['#ea9999', '#9edae5', '#FFE69F']
	fig, axs = plt.subplots(1, 3, figsize=(20,10), sharey=True)

	for house_idx, (house_data_set, house_label) in enumerate(zip(house_data, house_labels)):
		RPKM_vals = [data for data in house_data_set]
		time_points = range(1, 5)
			
		ax = axs[house_idx]
		ax.bar(time_points, RPKM_vals, align='center', color=house_colors[house_idx])
		
		for time_point, total_RPKM in enumerate(RPKM_vals):
			ax.text(time_point + 1, total_RPKM + max(RPKM_vals) * 0.05, f'{total_RPKM:.2f}', va='bottom', fontsize=10, ha='center')

		ax.set_title(f'House {house_label}')
		ax.set_xlabel('Time point')
		ax.set_xticks(time_points)
		ax.set_xticklabels(([f'T{tp}' for tp in time_points]))
		ax.set_ylabel('Total RPKM' if house_idx == 0 else "")

	fig.suptitle(f"Total abundance of plasmid types")
	plt.tight_layout(rect=[0, 0, 1, 0.96])
	plt.savefig(f"/media/ja/data/Mikra/Jane/plasmidome/new/Results/combined_total_RPKM_MOB.svg")

#Data
H1_data = [H1_13_RPKM, H1_20_RPKM, H1_27_RPKM, H1_35_RPKM]
H2_data = [H2_10_RPKM, H2_17_RPKM, H2_24_RPKM, H2_31_RPKM]
H3_data = [H3_9_RPKM, H3_16_RPKM, H3_23_RPKM, H3_30_RPKM]

all_data = [H1_data, H2_data, H3_data]
house_names = ["1", "2", "3"]

#Visualize the data in plots
visualise_bar(all_data, house_names)