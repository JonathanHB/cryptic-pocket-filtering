import numpy as np
import os

#get the index (and holo chain) which filter_blast2pairs.py associates with a given apo or holo chain

directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/new_pockets/iofiles"
blast_out = os.listdir(f"{directory}/blast_output/") #filenames for blastp results of candidate holo (MOAD) structures with at least one candidate apo (PDB blastp hit)structure

#look for a holo id
holo_query = "1glg"
print(blast_out.index(f"{holo_query.upper()}_hits.npy")) #get index from holo_id

#look for the first occurrence of an apo id
apo_query = "4pyk"
apo_query = apo_query.lower()

#index range to search
rstart = 0
rend = len(blast_out)

for holo_candidate in range(rstart,rend): #loop through candidate holo structures
    holo_id = blast_out[holo_candidate][0:4].upper() #get id at index
    apo_ids = [id.lower() for id in np.load(f"{directory}/blast_output/{holo_id}_hits.npy")] #load blast results
    if apo_query in apo_ids:
        print(holo_id.lower())
        print(holo_candidate)
        break
