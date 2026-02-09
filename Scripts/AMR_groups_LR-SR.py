import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

##Define functions you will use
def read_AMR(AMR_table, sampel):
	#read table with specified AMRs
	lines = open(AMR_table, "r", encoding="ISO-8859-1").readlines()

	ID_AMR = []
	group_AMR = []

	for line in lines[1:]:
		ID_AMR.append(line.split("\t")[0])
		group_AMR.append(line.split("\t")[3].strip().replace("\n", ""))

	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/LR-SR_combined_AMR-only_new.txt', 'w') as data:
		data.write(f'The Gene-Group pairs for LR of {sampel} are:\n')
		data.write(f'{ID_AMR}\t{group_AMR}\n')

	return [ID_AMR, group_AMR]

def read_AMR_SR(AMR_table, sampel):
	#read table with specified AMRs
	lines = open(AMR_table, "r").readlines()

	ID_AMR = []
	group_AMR = []

	for line in lines[1:]:
		ID_AMR.append(line.split("\t")[0])
		group_AMR.append(line.split("\t")[3].strip().replace("\n", ""))

	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/LR-SR_combined_AMR-only_new.txt', 'a') as data:
		data.write(f'The Gene-Group pairs for SR of {sampel} are:\n')
		data.write(f'{ID_AMR}\t{group_AMR}\n')

	return [ID_AMR, group_AMR]

def kma_res(res_table, AMR_list):
	#read the kma results, specify where the values are positioned in res file
	lines = open(res_table, "r").readlines()
	
	ID, cov, identity, length = [], [], [], []
	for line in lines[1:]:
		if line.split()[0] in AMR_list:
			ID.append(line.split("\t")[0])
			cov.append(float(line.split("\t")[5]))
			identity.append(float(line.split("\t")[4]))			
			length.append(float(line.split("\t")[3]))

	return np.array(ID), np.array(cov), np.array(identity), np.array(length)
	# return the values in the form of arrays

def kma_mapstat(mapstat_table, ID_List, sampel):
	#return the ids and correponding readcounts (fragmentCountAln) from mapstat file 
	lines1 = open(mapstat_table, "r").readlines()
	total_length = float(lines1[3].split()[2].strip().replace("#", ""))

	readCount = []
	IDs = []
	for line in lines1[7:]:
		if line.split()[0] in ID_List:
			IDs.append(line.split()[0])
			readCount.append(float(line.split()[14]))

	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/LR-SR_combined_AMR-only_new.txt', 'a') as data:
		data.write('mapstat readcounts for long-read data of %s are:\n' % (sampel))
		data.write(f'{IDs}\n')
		data.write(f'{readCount}\n')
			
	return [IDs, readCount], total_length 

def kma_mapstat_SR(mapstat_table, ID_List, sampel):
	#return the ids and correponding readcounts (fragmentCountAln) from mapstat file 
	lines1 = open(mapstat_table, "r").readlines()
	total_length = float(lines1[3].split()[2].strip().replace("#", ""))

	readCount = []
	IDs = []
	for line in lines1[7:]:
		if line.split()[0] in ID_List:
			IDs.append(line.split()[0])
			readCount.append(float(line.split()[14]))

	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/LR-SR_combined_AMR-only_new.txt', 'a') as data:
		data.write('mapstat readcounts for Short-read data of %s are:\n' % (sampel))
		data.write(f'{IDs}\n')
		data.write(f'{readCount}\n')
			
	return [IDs, readCount], total_length 

def cov80(ids, covs, identity, lengths, sampel):
	#find hits with coverage greater than 80

	covs_80_find = np.where((covs >= 80) & (identity >= 90))
	ids_80 = ids[covs_80_find]
	covs_80 = covs[covs_80_find]
	identity_90 = identity[covs_80_find]
	lengths_80 = lengths[covs_80_find]

	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/LR-SR_combined_AMR-only_new.txt', 'a') as data:
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

	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/LR-SR_combined_AMR-only_new.txt', 'a') as data:
		data.write('RKPM values for %s are:\n' % (sampel))
		data.write(f'{rpkm_id}\n')
		data.write(f'{RPKMs}\n')
		data.write(f'{RPKM_readsCount}\n')
		data.write(f'{RPKM_Length}\n')
		data.write(f'{totalNumreads}\n')
		data.write(f'{total_RPKM}')

	return [rpkm_id, RPKMs], total_RPKM

def AMR80(amr_list, cov80_id, sampel):
	#find only the genes present in cov80
	ID_AMR_80, group_AMR_80 = [], []

	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/LR-SR_combined_AMR-only_new.txt', 'a') as data:
		data.write('The AMR-only lists for %s are:\n' % (sampel))
		data.write(f'{amr_list}\n')
		data.write(f'{cov80_id}\n')

	amr_list_id = amr_list[0]
	amr_list_group = amr_list[1]

	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/LR-SR_combined_AMR-only_new.txt', 'a') as data:
		data.write('The lists which are used to filter out AMR genes only for %s are:\n' % (sampel))
		data.write(f'{amr_list_id}\n')
		data.write(f'{amr_list_group}\n')
		data.write(f'{cov80_id}\n')

	for id_amr, group_amr in zip(amr_list_id, amr_list_group):
		if id_amr in cov80_id:
			ID_AMR_80.append(id_amr)
			group_AMR_80.append(group_amr)

	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/LR-SR_combined_AMR-only_new.txt', 'a') as data:
		data.write('The lists of updated AMR-only table for %s are:\n' % (sampel))
		data.write(f'{ID_AMR_80}\n')
		data.write(f'{group_AMR_80}\n')

	return [ID_AMR_80, group_AMR_80]    

def sort_AMR_groups(RPKM_list, AMR_list, total, samp):
    # Create mappings from gene IDs to RPKM values and AMR groups
    gene_id_to_rpkm = dict(zip(RPKM_list[0], RPKM_list[1]))
    gene_id_to_group = dict(zip(AMR_list[0], AMR_list[1]))
    
    # Initialize dictionary to hold RPKM values for each group
    group_values = {}
    group_sums = {}
    
    # Aggregate RPKM values under their respective AMR groups
    for gene_id, group in gene_id_to_group.items():
        rpkm = gene_id_to_rpkm.get(gene_id)
        if rpkm is not None:
            group_values.setdefault(group, []).append(float(rpkm))
    
    # Calculate counts and average RPKM values for each group
    group_sums = {group: sum(value) for group, value in group_values.items()}
    group_avg = {
        group: sum(values) / total if values else 0.0
        for group, values in group_values.items()
    }
    
    # Write results to the specified file
    with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/LR-SR_combined_AMR-only_new.txt', 'a') as data:
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

    for grup, val in zip(group_list[0], group_list[1]):
    	with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/final_LR-SR_combined_AMR-only.txt', 'a') as table:
    		table.write(f'{samp}\t{grup}\t{val}\t{group_sums}\t{total}\n')
    
    with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/LR-SR_combined_AMR-only_new.txt', 'a') as data:
        data.write(f'The group list for {samp} is:\n')
        data.write(f'{group_list}\n')
    
    return group_list
		
##Execute the analysis of RPKM
#Read AMR tables
#LR part
H1_13_AMR = read_AMR("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/grep_mod_AMR_only/H1_13d_o200_grep_mod_AMR_only.txt", 'H1_13')
H1_20_AMR = read_AMR("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/grep_mod_AMR_only/H1_20d_o200_grep_mod_AMR_only.txt", 'H1_20')
H1_27_AMR = read_AMR("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/grep_mod_AMR_only/H1_27d_o200_grep_mod_AMR_only.txt", 'H1_27')
H1_35_AMR = read_AMR("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/grep_mod_AMR_only/H1_35d_o200_grep_mod_AMR_only.txt", 'H1_35')
H2_10_AMR = read_AMR("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/grep_mod_AMR_only/H2_10d_o200_grep_mod_AMR_only.txt", 'H2_10')
H2_17_AMR = read_AMR("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/grep_mod_AMR_only/H2_17d_o200_grep_mod_AMR_only.txt", 'H2_17')
H2_24_AMR = read_AMR("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/grep_mod_AMR_only/H2_24d_o200_grep_mod_AMR_only.txt", 'H2_24')
H2_31_AMR = read_AMR("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/grep_mod_AMR_only/H2_31d_o200_grep_mod_AMR_only.txt", 'H2_31')
H3_9_AMR = read_AMR("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/grep_mod_AMR_only/H3_9d_o200_grep_mod_AMR_only.txt", 'H3_9')
H3_16_AMR = read_AMR("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/grep_mod_AMR_only/H3_16d_o200_grep_mod_AMR_only.txt", 'H3_16')
H3_23_AMR = read_AMR("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/grep_mod_AMR_only/H3_23d_o200_grep_mod_AMR_only.txt", 'H3_23')
H3_30_AMR = read_AMR("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/grep_mod_AMR_only/H3_30d_o200_grep_mod_AMR_only.txt", 'H3_30')

H1_13_ID_AMR, H1_13_gr_AMR = H1_13_AMR[0], H1_13_AMR[1]
H1_20_ID_AMR, H1_20_gr_AMR = H1_20_AMR[0], H1_20_AMR[1]
H1_27_ID_AMR, H1_27_gr_AMR = H1_27_AMR[0], H1_27_AMR[1]
H1_35_ID_AMR, H1_35_gr_AMR = H1_35_AMR[0], H1_35_AMR[1]
H2_10_ID_AMR, H2_10_gr_AMR = H2_10_AMR[0], H2_10_AMR[1]
H2_17_ID_AMR, H2_17_gr_AMR = H2_17_AMR[0], H2_17_AMR[1]
H2_24_ID_AMR, H2_24_gr_AMR = H2_24_AMR[0], H2_24_AMR[1]
H2_31_ID_AMR, H2_31_gr_AMR = H2_31_AMR[0], H2_31_AMR[1]
H3_9_ID_AMR, H3_9_gr_AMR = H3_9_AMR[0], H3_9_AMR[1]
H3_16_ID_AMR, H3_16_gr_AMR = H3_16_AMR[0], H3_16_AMR[1]
H3_23_ID_AMR, H3_23_gr_AMR = H3_23_AMR[0], H3_23_AMR[1]
H3_30_ID_AMR, H3_30_gr_AMR = H3_30_AMR[0], H3_30_AMR[1]

#SR part
h1_13_AMR = read_AMR_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/02F_grep_mod_AMR_only.txt", 'H1_13')
h1_20_AMR = read_AMR_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/03F_grep_mod_AMR_only.txt", 'H1_20')
h1_27_AMR = read_AMR_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/04F_grep_mod_AMR_only.txt", 'H1_27')
h1_35_AMR = read_AMR_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/05F_grep_mod_AMR_only.txt", 'H1_35')
h2_10_AMR = read_AMR_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/07F_grep_mod_AMR_only.txt", 'H2_10')
h2_17_AMR = read_AMR_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/08F_grep_mod_AMR_only.txt", 'H2_17')
h2_24_AMR = read_AMR_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/09F_grep_mod_AMR_only.txt", 'H2_24')
h2_31_AMR = read_AMR_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/10F_grep_mod_AMR_only.txt", 'H2_31')
h3_9_AMR = read_AMR_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/12F_grep_mod_AMR_only.txt", 'H3_9')
h3_16_AMR = read_AMR_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/13F_grep_mod_AMR_only.txt", 'H3_16')
h3_23_AMR = read_AMR_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/14F_grep_mod_AMR_only.txt", 'H3_23')
h3_30_AMR = read_AMR_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/15F_grep_mod_AMR_only.txt", 'H3_30')

h1_13_ID_AMR, h1_13_gr_AMR = h1_13_AMR[0], h1_13_AMR[1]
h1_20_ID_AMR, h1_20_gr_AMR = h1_20_AMR[0], h1_20_AMR[1]
h1_27_ID_AMR, h1_27_gr_AMR = h1_27_AMR[0], h1_27_AMR[1]
h1_35_ID_AMR, h1_35_gr_AMR = h1_35_AMR[0], h1_35_AMR[1]
h2_10_ID_AMR, h2_10_gr_AMR = h2_10_AMR[0], h2_10_AMR[1]
h2_17_ID_AMR, h2_17_gr_AMR = h2_17_AMR[0], h2_17_AMR[1]
h2_24_ID_AMR, h2_24_gr_AMR = h2_24_AMR[0], h2_24_AMR[1]
h2_31_ID_AMR, h2_31_gr_AMR = h2_31_AMR[0], h2_31_AMR[1]
h3_9_ID_AMR, h3_9_gr_AMR = h3_9_AMR[0], h3_9_AMR[1]
h3_16_ID_AMR, h3_16_gr_AMR = h3_16_AMR[0], h3_16_AMR[1]
h3_23_ID_AMR, h3_23_gr_AMR = h3_23_AMR[0], h3_23_AMR[1]
h3_30_ID_AMR, h3_30_gr_AMR = h3_30_AMR[0], h3_30_AMR[1]

#Read kma .res files and assign position for ID, coverage and length
#LR part
H1_13d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H1_13d_o200_panres.res", H1_13_ID_AMR)
H1_20d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H1_20d_o200_panres.res", H1_20_ID_AMR)
H1_27d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H1_27d_o200_panres.res", H1_27_ID_AMR)
H1_35d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H1_35d_o200_panres.res", H1_35_ID_AMR)
H2_10d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H2_10d_o200_panres.res", H2_10_ID_AMR)
H2_17d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H2_17d_o200_panres.res", H2_17_ID_AMR)
H2_24d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H2_24d_o200_panres.res", H2_24_ID_AMR)
H2_31d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H2_31d_o200_panres.res", H2_31_ID_AMR)
H3_9d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H3_9d_o200_panres.res", H3_9_ID_AMR)
H3_16d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H3_16d_o200_panres.res", H3_16_ID_AMR)
H3_23d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H3_23d_o200_panres.res", H3_23_ID_AMR)
H3_30d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H3_30d_o200_panres.res", H3_30_ID_AMR)

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

#SR part
h1_13d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/02F_panres.res", h1_13_ID_AMR)
h1_20d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/03F_panres.res", h1_20_ID_AMR)
h1_27d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/04F_panres.res", h1_27_ID_AMR)
h1_35d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/05F_panres.res", h1_35_ID_AMR)
h2_10d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/07F_panres.res", h2_10_ID_AMR)
h2_17d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/08F_panres.res", h2_17_ID_AMR)
h2_24d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/09F_panres.res", h2_24_ID_AMR)
h2_31d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/10F_panres.res", h2_31_ID_AMR)
h3_9d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/12F_panres.res", h3_9_ID_AMR)
h3_16d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/13F_panres.res", h3_16_ID_AMR)
h3_23d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/14F_panres.res", h3_23_ID_AMR)
h3_30d = kma_res("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/15F_panres.res", h3_30_ID_AMR)

h1_13_id, h1_13_cov, h1_13_IDidentity, h1_13_IDlength = h1_13d[0], h1_13d[1], h1_13d[2], h1_13d[3]
h1_20_id, h1_20_cov, h1_20_IDidentity, h1_20_IDlength = h1_20d[0], h1_20d[1], h1_20d[2], h1_20d[3]
h1_27_id, h1_27_cov, h1_27_IDidentity, h1_27_IDlength  = h1_27d[0], h1_27d[1], h1_27d[2], h1_27d[3]
h1_35_id, h1_35_cov, h1_35_IDidentity, h1_35_IDlength = h1_35d[0], h1_35d[1], h1_35d[2], h1_35d[3]
h2_10_id, h2_10_cov, h2_10_IDidentity, h2_10_IDlength = h2_10d[0], h2_10d[1], h2_10d[2], h2_10d[3]
h2_17_id, h2_17_cov, h2_17_IDidentity, h2_17_IDlength = h2_17d[0], h2_17d[1], h2_17d[2], h2_17d[3]
h2_24_id, h2_24_cov, h2_24_IDidentity, h2_24_IDlength = h2_24d[0], h2_24d[1], h2_24d[2], h2_24d[3]
h2_31_id, h2_31_cov, h2_31_IDidentity, h2_31_IDlength = h2_31d[0], h2_31d[1], h2_31d[2], h2_31d[3]
h3_9_id, h3_9_cov, h3_9_IDidentity, h3_9_IDlength = h3_9d[0], h3_9d[1], h3_9d[2], h3_9d[3]
h3_16_id, h3_16_cov, h3_16_IDidentity, h3_16_IDlength = h3_16d[0], h3_16d[1], h3_16d[2], h3_16d[3]
h3_23_id, h3_23_cov, h3_23_IDidentity, h3_23_IDlength = h3_23d[0], h3_23d[1], h3_23d[2], h3_23d[3]
h3_30_id, h3_30_cov, h3_30_IDidentity, h3_30_IDlength = h3_30d[0], h3_30d[1], h3_30d[2], h3_30d[3]

#Filter out the data with coverage =< 80
#LR
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

#SR
h1_13_id_80, h1_13_cov_80, h1_13_IDidentity_90, h1_13_IDlength_80 = cov80(h1_13_id, h1_13_cov, h1_13_IDidentity, h1_13_IDlength, 'H1_13')
h1_20_id_80, h1_20_cov_80, h1_20_IDidentity_90, h1_20_IDlength_80 = cov80(h1_20_id, h1_20_cov, h1_20_IDidentity, h1_20_IDlength, 'H1_20')
h1_27_id_80, h1_27_cov_80, h1_27_IDidentity_90, h1_27_IDlength_80 = cov80(h1_27_id, h1_27_cov, h1_27_IDidentity, h1_27_IDlength, 'H1_27')
h1_35_id_80, h1_35_cov_80, h1_35_IDidentity_90, h1_35_IDlength_80 = cov80(h1_35_id, h1_35_cov, h1_35_IDidentity, h1_35_IDlength, 'H1_35')
h2_10_id_80, h2_10_cov_80, h2_10_IDidentity_90, h2_10_IDlength_80 = cov80(h2_10_id, h2_10_cov, h2_10_IDidentity, h2_10_IDlength, 'H2_10')
h2_17_id_80, h2_17_cov_80, h2_17_IDidentity_90, h2_17_IDlength_80 = cov80(h2_17_id, h2_17_cov, h2_17_IDidentity, h2_17_IDlength, 'H2_17')
h2_24_id_80, h2_24_cov_80, h2_24_IDidentity_90, h2_24_IDlength_80 = cov80(h2_24_id, h2_24_cov, h2_24_IDidentity, h2_24_IDlength, 'H2_24')
h2_31_id_80, h2_31_cov_80, h2_31_IDidentity_90, h2_31_IDlength_80 = cov80(h2_31_id, h2_31_cov, h2_31_IDidentity, h2_31_IDlength, 'H2_31')
h3_9_id_80, h3_9_cov_80, h3_9_IDidentity_90, h3_9_IDlength_80 = cov80(h3_9_id, h3_9_cov, h3_9_IDidentity, h3_9_IDlength, 'H3_9')
h3_16_id_80, h3_16_cov_80, h3_16_IDidentity_90, h3_16_IDlength_80 = cov80(h3_16_id, h3_16_cov, h3_16_IDidentity, h3_16_IDlength, 'H3_16')
h3_23_id_80, h3_23_cov_80, h3_23_IDidentity_90, h3_23_IDlength_80 = cov80(h3_23_id, h3_23_cov, h3_23_IDidentity, h3_23_IDlength, 'H3_23')
h3_30_id_80, h3_30_cov_80, h3_30_IDidentity_90, h3_30_IDlength_80 = cov80(h3_30_id, h3_30_cov, h3_30_IDidentity, h3_30_IDlength, 'H3_30')

#Read kma .mapstat files and assing position of ID + read count, and total length
#LR
H1_13_numread_pair, H1_13_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H1_13d_o200_panres.mapstat", H1_13_id_80, 'H1_13')
H1_20_numread_pair, H1_20_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H1_20d_o200_panres.mapstat", H1_20_id_80, 'H1_20')
H1_27_numread_pair, H1_27_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H1_27d_o200_panres.mapstat", H1_27_id_80, 'H1_27')
H1_35_numread_pair, H1_35_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H1_35d_o200_panres.mapstat", H1_35_id_80, 'H1_35')
H2_10_numread_pair, H2_10_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H2_10d_o200_panres.mapstat", H2_10_id_80, 'H2_10')
H2_17_numread_pair, H2_17_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H2_17d_o200_panres.mapstat", H2_17_id_80, 'H2_17')
H2_24_numread_pair, H2_24_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H2_24d_o200_panres.mapstat", H2_24_id_80, 'H2_24')
H2_31_numread_pair, H2_31_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H2_31d_o200_panres.mapstat", H2_31_id_80, 'H2_31')
H3_9_numread_pair, H3_9_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H3_9d_o200_panres.mapstat", H3_9_id_80, 'H3_9')
H3_16_numread_pair, H3_16_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H3_16d_o200_panres.mapstat", H3_16_id_80, 'H3_16')
H3_23_numread_pair, H3_23_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H3_23d_o200_panres.mapstat", H3_23_id_80, 'H3_23')
H3_30_numread_pair, H3_30_totalLength = kma_mapstat("/media/ja/data/Mikra/Jane/plasmidome/new/Marky/Long_read/H3_30d_o200_panres.mapstat", H3_30_id_80, 'H3_30')

#SR
h1_13_numread_pair, h1_13_totalLength = kma_mapstat_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/02F_panres.mapstat", h1_13_id_80, 'h1_13')
h1_20_numread_pair, h1_20_totalLength = kma_mapstat_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/03F_panres.mapstat", h1_20_id_80, 'h1_20')
h1_27_numread_pair, h1_27_totalLength = kma_mapstat_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/04F_panres.mapstat", h1_27_id_80, 'h1_27')
h1_35_numread_pair, h1_35_totalLength = kma_mapstat_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/05F_panres.mapstat", h1_35_id_80, 'h1_35')
h2_10_numread_pair, h2_10_totalLength = kma_mapstat_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/07F_panres.mapstat", h2_10_id_80, 'h2_10')
h2_17_numread_pair, h2_17_totalLength = kma_mapstat_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/08F_panres.mapstat", h2_17_id_80, 'h2_17')
h2_24_numread_pair, h2_24_totalLength = kma_mapstat_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/09F_panres.mapstat", h2_24_id_80, 'h2_24')
h2_31_numread_pair, h2_31_totalLength = kma_mapstat_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/10F_panres.mapstat", h2_31_id_80, 'h2_31')
h3_9_numread_pair, h3_9_totalLength = kma_mapstat_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/12F_panres.mapstat", h3_9_id_80, 'h3_9')
h3_16_numread_pair, h3_16_totalLength = kma_mapstat_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/13F_panres.mapstat", h3_16_id_80, 'h3_16')
h3_23_numread_pair, h3_23_totalLength = kma_mapstat_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/14F_panres.mapstat", h3_23_id_80, 'h3_23')
h3_30_numread_pair, h3_30_totalLength = kma_mapstat_SR("/media/ja/data/Mikra/Jane/plasmidome/new/Short-read/trimmed/15F_panres.mapstat", h3_30_id_80,'h3_30')

#Count the RPKMs
#LR
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

#SR
h1_13_RPKM, h1_13_tot_RPKM = RPKM_count(h1_13_id_80, h1_13_numread_pair, h1_13_IDlength_80, h1_13_totalLength, 'h1_13')
h1_20_RPKM, h1_20_tot_RPKM = RPKM_count(h1_20_id_80, h1_20_numread_pair, h1_20_IDlength_80, h1_20_totalLength, 'h1_20')
h1_27_RPKM, h1_27_tot_RPKM = RPKM_count(h1_27_id_80, h1_27_numread_pair, h1_27_IDlength_80, h1_27_totalLength, 'h1_27')
h1_35_RPKM, h1_35_tot_RPKM = RPKM_count(h1_35_id_80, h1_35_numread_pair, h1_35_IDlength_80, h1_35_totalLength, 'h1_35')
h2_10_RPKM, h2_10_tot_RPKM = RPKM_count(h2_10_id_80, h2_10_numread_pair, h2_10_IDlength_80, h2_10_totalLength, 'h2_10')
h2_17_RPKM, h2_17_tot_RPKM = RPKM_count(h2_17_id_80, h2_17_numread_pair, h2_17_IDlength_80, h2_17_totalLength, 'h2_17')
h2_24_RPKM, h2_24_tot_RPKM = RPKM_count(h2_24_id_80, h2_24_numread_pair, h2_24_IDlength_80, h2_24_totalLength, 'h2_24')
h2_31_RPKM, h2_31_tot_RPKM = RPKM_count(h2_31_id_80, h2_31_numread_pair, h2_31_IDlength_80, h2_31_totalLength, 'h2_31')
h3_9_RPKM, h3_9_tot_RPKM = RPKM_count(h3_9_id_80, h3_9_numread_pair, h3_9_IDlength_80, h3_9_totalLength, 'h3_9')
h3_16_RPKM, h3_16_tot_RPKM = RPKM_count(h3_16_id_80, h3_16_numread_pair, h3_16_IDlength_80, h3_16_totalLength, 'h3_16')
h3_23_RPKM, h3_23_tot_RPKM = RPKM_count(h3_23_id_80, h3_23_numread_pair, h3_23_IDlength_80, h3_23_totalLength, 'h3_23')
h3_30_RPKM, h3_30_tot_RPKM = RPKM_count(h3_30_id_80, h3_30_numread_pair, h3_30_IDlength_80, h3_30_totalLength, 'h3_30')

#Put together AMR table information for id, cov 80 
#LR
H1_13_AMR_80 = AMR80(H1_13_AMR, H1_13_id_80, "H1 13")
H1_20_AMR_80 = AMR80(H1_20_AMR, H1_20_id_80, "H1 20")
H1_27_AMR_80 = AMR80(H1_27_AMR, H1_27_id_80, "H1 27")
H1_35_AMR_80 = AMR80(H1_35_AMR, H1_35_id_80, "H1 35")
H2_10_AMR_80 = AMR80(H2_10_AMR, H2_10_id_80, "H2 10")
H2_17_AMR_80 = AMR80(H2_17_AMR, H2_17_id_80, "H2 17")
H2_24_AMR_80 = AMR80(H2_24_AMR, H2_24_id_80, "H2 24")
H2_31_AMR_80 = AMR80(H2_31_AMR, H2_31_id_80, "H2 31")
H3_9_AMR_80 = AMR80(H3_9_AMR, H3_9_id_80, "H3 9")
H3_16_AMR_80 = AMR80(H3_16_AMR, H3_16_id_80, "H3 16")
H3_23_AMR_80 = AMR80(H3_23_AMR, H3_23_id_80, "H3 23")
H3_30_AMR_80 = AMR80(H3_30_AMR, H3_30_id_80, "H3 30")

#SR
h1_13_AMR_80 = AMR80(h1_13_AMR, h1_13_id_80, "h1 13")
h1_20_AMR_80 = AMR80(h1_20_AMR, h1_20_id_80, "h1 20")
h1_27_AMR_80 = AMR80(h1_27_AMR, h1_27_id_80, "h1 27")
h1_35_AMR_80 = AMR80(h1_35_AMR, h1_35_id_80, "h1 35")
h2_10_AMR_80 = AMR80(h2_10_AMR, h2_10_id_80, "h2 10")
h2_17_AMR_80 = AMR80(h2_17_AMR, h2_17_id_80, "h2 17")
h2_24_AMR_80 = AMR80(h2_24_AMR, h2_24_id_80, "h2 24")
h2_31_AMR_80 = AMR80(h2_31_AMR, h2_31_id_80, "h2 31")
h3_9_AMR_80 = AMR80(h3_9_AMR, h3_9_id_80, "h3 9")
h3_16_AMR_80 = AMR80(h3_16_AMR, h3_16_id_80, "h3 16")
h3_23_AMR_80 = AMR80(h3_23_AMR, h3_23_id_80, "h3 23")
h3_30_AMR_80 = AMR80(h3_30_AMR, h3_30_id_80, "h3 30")

#Put together all RPKMs and sort them alphabetically
with open('/media/ja/data/Mikra/Jane/plasmidome/new/Results/final_LR-SR_combined_AMR-only.txt', 'w') as table:
	table.write(f'Sample ID\tATB group\tGroup AVG (%)\tRPKM sum\ttotal RPKM\n')

#LR
H1_13_groups = sort_AMR_groups(H1_13_RPKM, H1_13_AMR_80, H1_13_tot_RPKM, "H1 13")
H1_20_groups = sort_AMR_groups(H1_20_RPKM, H1_20_AMR_80, H1_20_tot_RPKM, "H1 20")
H1_27_groups = sort_AMR_groups(H1_27_RPKM, H1_27_AMR_80, H1_27_tot_RPKM, "H1 27")
H1_35_groups = sort_AMR_groups(H1_35_RPKM, H1_35_AMR_80, H1_35_tot_RPKM, "H1 35")
H2_10_groups = sort_AMR_groups(H2_10_RPKM, H2_10_AMR_80, H2_10_tot_RPKM, "H2 10")
H2_17_groups = sort_AMR_groups(H2_17_RPKM, H2_17_AMR_80, H2_17_tot_RPKM, "H2 17")
H2_24_groups = sort_AMR_groups(H2_24_RPKM, H2_24_AMR_80, H2_24_tot_RPKM, "H2 24")
H2_31_groups = sort_AMR_groups(H2_31_RPKM, H2_31_AMR_80, H2_31_tot_RPKM, "H2 31")
H3_9_groups = sort_AMR_groups(H3_9_RPKM, H3_9_AMR_80, H3_9_tot_RPKM, "H3 9")
H3_16_groups = sort_AMR_groups(H3_16_RPKM, H3_16_AMR_80, H3_16_tot_RPKM, "H3 16")
H3_23_groups = sort_AMR_groups(H3_23_RPKM, H3_23_AMR_80, H3_23_tot_RPKM, "H3 23")
H3_30_groups = sort_AMR_groups(H3_30_RPKM, H3_30_AMR_80, H3_30_tot_RPKM, "H3 30")

#SR
h1_13_groups = sort_AMR_groups(h1_13_RPKM, h1_13_AMR_80, h1_13_tot_RPKM, "h1 13")
h1_20_groups = sort_AMR_groups(h1_20_RPKM, h1_20_AMR_80, h1_20_tot_RPKM, "h1 20")
h1_27_groups = sort_AMR_groups(h1_27_RPKM, h1_27_AMR_80, h1_27_tot_RPKM, "h1 27")
h1_35_groups = sort_AMR_groups(h1_35_RPKM, h1_35_AMR_80, h1_35_tot_RPKM, "h1 35")
h2_10_groups = sort_AMR_groups(h2_10_RPKM, h2_10_AMR_80, h2_10_tot_RPKM, "h2 10")
h2_17_groups = sort_AMR_groups(h2_17_RPKM, h2_17_AMR_80, h2_17_tot_RPKM, "h2 17")
h2_24_groups = sort_AMR_groups(h2_24_RPKM, h2_24_AMR_80, h2_24_tot_RPKM, "h2 24")
h2_31_groups = sort_AMR_groups(h2_31_RPKM, h2_31_AMR_80, h2_31_tot_RPKM, "h2 31")
h3_9_groups = sort_AMR_groups(h3_9_RPKM, h3_9_AMR_80, h3_9_tot_RPKM, "h3 9")
h3_16_groups = sort_AMR_groups(h3_16_RPKM, h3_16_AMR_80, h3_16_tot_RPKM, "h3 16")
h3_23_groups = sort_AMR_groups(h3_23_RPKM, h3_23_AMR_80, h3_23_tot_RPKM, "h3 23")
h3_30_groups = sort_AMR_groups(h3_30_RPKM, h3_30_AMR_80, h3_30_tot_RPKM, "h3 30")

'''def get_colors(n):
    """Generate a list of n unique colors."""
    cmap = cm.get_cmap('tab20', n)
    return [cmap(i) for i in range(n)]'''

color_mapping = {
    'tetracyclines': '#A9D4A0', 
    'aminocoumarins': '#dfbac9',  
    'aminoglycosides': '#B8DECC', 
    'betalactams': '#ea9999',  
    'polymyxins': '#b4a7d6',  
    'diaminopyrimidines': '#B37C71', 
    'elfamycins': '#C78BA5',  
    'fluoroquinolones': '#FFE69F',  
    'glycopeptides': '#bcbd22', 
    'ionofors': '#d2d7e7',  
    'lincosamides': '#ffcebe',  
    'macrolides': '#d5a6bd',  
    'fusidanes': '#9fc5e8',  
    'MDR': '#dda298',  
    'nitrofurans': '#eec5b7',  
    'nitroimidazoles': '#C0E6F5',   
    'nucleosides': '#6CB5E6',    
    'peptides': '#F2CEEF',  
    'phenicols': '#84D6B3',  
    'phosphonic acid': '#a05195',  
    'pleuromutilins': '#f95d6a',  
    'rifamycins': '#DF6F84',  
    'streptogramins': '#BCBD22',  
    'sulfonamides': '#ffe3c0',
    'orthosomycins': '#6CB5E6',  
    'fosfomycins': '#F95D6A'
    }

def visualisation_combined_bar(info_names, data_sets_LR, data_sets_SR, time_labels, time_points, color_mapping):
	plt.figure(figsize=(16, 10))
	n_halls = len(info_names)
	width = 0.5
	gap_between_houses = 1
	gap_between_points = 1.25

	categories = []
	x_positions = []

	base_x = 0
	for house_idx, house in enumerate(info_names):
		for point_idx, time_label in enumerate(time_labels[house_idx]):
			categories.append(f'{house}_{time_label} LR')
			x_positions.append(base_x - width / 2)

			categories.append(f'{house}_{time_label} SR')
			x_positions.append(base_x + width / 2)

			base_x += gap_between_points

		base_x += gap_between_houses

	x = np.array(x_positions)

	# Track which AMR groups have already been added to the legend
	legend_labels = set()

	all_data = {key: np.zeros(len(categories)) for key in color_mapping.keys()}
	
	# Prepare to plot each time point
	index = 0
	for house_idx in range(n_halls):
		for point_idx in range(len(time_labels[house_idx])):
			LR_hall_data = data_sets_LR[point_idx][house_idx]
			SR_hall_data = data_sets_SR[point_idx][house_idx]

			for info_id, value in zip(LR_hall_data[0], LR_hall_data[1]):
				if info_id in all_data:
					all_data[info_id][index] += value

			index += 1 

			for info_id, value in zip(SR_hall_data[0], SR_hall_data[1]):
				if info_id in all_data:
					all_data[info_id][index] += value

			index += 1

	normalized_data = np.array([all_data[key] for key in all_data]) * 100

	bottom = np.zeros(len(categories))
	for i, (info_id, heights) in enumerate(zip(all_data.keys(), normalized_data)):
	    color = color_mapping.get(info_id, '#000000')

	    # Přidání bílého ohraničení kolem hodnot
	    plt.bar(x, heights, width, color=color, bottom=bottom, edgecolor='white', linewidth=0.5)
	    
	    # Přidání popisku do legendy pouze pokud tam ještě není
	    if info_id not in legend_labels:
	        plt.bar(x, heights, width, label=f'{info_id}', color=color, bottom=bottom, edgecolor='white', linewidth=0.5)
	        legend_labels.add(info_id)
	    bottom += heights

	plt.xlabel("Houses, Collection Points and Data Type")
	plt.ylabel("Relative Abundance (%)")
	plt.title("Relative Abundance of AMR Groups: Plasmidome vs Metagenome")
	plt.xticks(x, categories, rotation=45, ha='right')
	plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	plt.tight_layout()

	plt.savefig(f"/media/ja/data/Mikra/Jane/plasmidome/new/Results/LR-SR_combined_panres_new.svg")  

#Data
info_names = ["H1", "H2", "H3"]
time_labels = [
    ["13d", "20d", "27d", "35d"],  # House 1
    ["10d", "17d", "24d", "31d"],  # House 2
    ["9d", "16d", "23d", "30d"]    # House 3
]

time_points = ["1", "2", "3", "4"]
data_sets_LR = [
    [H1_13_groups, H2_10_groups, H3_9_groups],
    [H1_20_groups, H2_17_groups, H3_16_groups],
    [H1_27_groups, H2_24_groups, H3_23_groups],
    [H1_35_groups, H2_31_groups, H3_30_groups],
]
data_sets_SR = [
    [h1_13_groups, h2_10_groups, h3_9_groups],
    [h1_20_groups, h2_17_groups, h3_16_groups],
    [h1_27_groups, h2_24_groups, h3_23_groups],
    [h1_35_groups, h2_31_groups, h3_30_groups],
]

visualisation_combined_bar(info_names, data_sets_LR, data_sets_SR, time_labels, time_points, color_mapping)