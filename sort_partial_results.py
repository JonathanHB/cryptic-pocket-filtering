import numpy as np
import os

#take the pairs saved to pairing_index/ and sort them in descending order by minimum rmsd to each holo structure
#this produces a list which begins with pairs with the largest conformational change not captured by experimental apo variation
#when multiple holo structures have the same apo structure, only the highest-rmsd pair is returned
#experimental holo variation can be attributable to differences in ligand size, and this code does not attempt to separate this from other factors

manual_removal = ["5x93", "1q6m", "1q6j", "4acx", "4b70", "4b05", "1yp9", "3ufl"]
#structures that are messed up in rare ways such that they fail silently and have to be removed manually
#5x93 has an incomplete remark 465
#"1q6m", "1q6j", "4acx", 4b70", 1yp9, 3ufl, and 4b05 have their residues numbered out of order
#in the case of 4b05 this also causes the gap-checker to fail since it believes
#that the remark 465 residues are all before the start of the structure after reading 498 as the first residue

#2wkp and 1opk are each missing half the structure

directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/new_pockets/iofiles"

holo_list = os.listdir(f"{directory}/pairing_index")

exclude_serial = -1 #exclude structures already compiled into a previous list; negative number to exclude nothing

if exclude_serial >= 0: #negative numbers compile complete lists
    already_screened = np.load(f"{directory}/output_indices/some_lowrmsd_ligand_pairs_v{exclude_serial}.npy")
    #print(already_screened)
    screened_hids = [i[0][0:4] for i in already_screened]
else:
    screened_hids = []
pair_ligands = []

for i in holo_list:#[2:]:
    if i[0:6] == "ligand":
        pair = np.load(f"{directory}/pairing_index/{i}")
        if len(pair[0][5])!=0 and (pair[0][0] not in manual_removal) and (pair[0][0] not in screened_hids): #the latter condition filters out fake structures included due to a now-resolved bug
            pair_ligands.append(np.load(f"{directory}/pairing_index/{i}"))


#print(pair_ligands)
pair_ligands = sorted(pair_ligands, key = lambda x: x[0][4], reverse = True)

#copy the first (and therefore highest RMSD since the array is sorted) apo-holo pair with each apo structure into a new array
apo_ids = []
lowrmsd_pairs = []
for j in pair_ligands:
    if j[0][2] not in apo_ids:
        lowrmsd_pairs.append(j[0])
        apo_ids.append(j[0][2])

    #print(j)
    #print(j[0])

serial = 7 #REMEMBER to update
#print(f"pairing_index/some_lowrmsd_ligand_pairs_v{serial}")
np.save(f"{directory}/output_indices/some_lowrmsd_ligand_pairs_v{serial}", np.array(lowrmsd_pairs, dtype = object))
