import numpy as np
import os
import random

path = "/project/bowmanlab/borowsky.jonathan/FAST-cs/new_pockets/apo_structures_clean"

files = os.listdir("/project/bowmanlab/borowsky.jonathan/FAST-cs/new_pockets/pocketresis_refined")

output_test = []
output_validation = []

n_test = 26
n_validation = 17

vindex = 2 #for keeping track of different copies of the split while testing

drawing_indices = list(range(n_test+n_validation))
test_inds = random.sample(drawing_indices, n_test)
#print(test_inds)

for i in range(len(files)):
    filepath = f"{path}/{files[i][0:5]}_clean.pdb"
    #print(filepath)

    resis = np.load(f"/project/bowmanlab/borowsky.jonathan/FAST-cs/new_pockets/pocketresis_refined/{files[i]}")
    #print(resis)

    #will fail if numbering actually extends to -999
    init = -999
    final = -999

    buffer = ""
    for line in open(f'{path}/{files[i][0:5]}_clean.pdb'): #get pdb apo indices
        if init == -999 and line[0:4] == "ATOM":
            #print(line)
            init = int(line[22:26])
        elif line[0:3] == "TER":
            final = int(buffer)
            break

        buffer = line[22:26]

    length = final-init+1
    label_arr = np.zeros(length)
    resis_shifted = [i-init for i in resis]
    label_arr[resis_shifted] = 1
    #print(label_arr)

    if i in test_inds:
        output_test.append([filepath, label_arr])
    else:
        output_validation.append([filepath, label_arr])

#print(output_test)
#print("_---------------------------------------------------------------------------------------")
#print("_---------------------------------------------------------------------------------------")
#print(output_validation)

#print(len(np.load("/project/bowmanlab/borowsky.jonathan/FAST-cs/new_pockets/labels/new_pocket_labels_test.npy")))

#uncomment to save new labels
#np.save(f"new_pocket_labels_test_v{vindex}", output_test)
#np.save(f"new_pocket_labels_validation_v{vindex}", output_validation)
