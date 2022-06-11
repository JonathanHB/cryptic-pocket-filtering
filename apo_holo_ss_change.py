import numpy as np
import mdtraj as md
import itertools
import scipy
import sys

serial_in = 6 #which version of the pocket set to load
serial_out = 7 #which version of the pocket set to save
input_directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/new_pockets/iofiles"

save = False

distance_threshold = 0.5

all_pockets = np.load(f"{input_directory}/output_indices/some_lowrmsd_ligand_pairs_v{serial_in}.npy")

pdb_apo_ids = ["1kx9","1tvq","2hq8","2j1x","2oy4","2zby","2zku","3cab","3nx1","3p53","3rwv","4i92","6e5d","1y1a","3fvj","3qxw","4ic4","4w51","5h9a","6ypk"]

for testind in range(137):#len(all_pockets)):

    if all_pockets[testind][2] not in pdb_apo_ids: #only examine screened structures
        continue

    print(f"{testind}: {all_pockets[testind]}")

    #load structures
    holo_xtal = md.load(f"../monomer_holo/{all_pockets[testind][0]}_chain{all_pockets[testind][1]}.pdb")
    apo_xtal = md.load(f"../monomer_apo/{all_pockets[testind][2]}_chain{all_pockets[testind][3]}.pdb")

    holo_resseqs = np.unique([holo_xtal.top.atom(i).residue.resSeq for i in holo_xtal.top.select("protein")])
    apo_resseqs = np.unique([apo_xtal.top.atom(i).residue.resSeq for i in apo_xtal.top.select("protein")])

    shared_resseqs = np.intersect1d(holo_resseqs, apo_resseqs)

    holo_iis = np.concatenate([holo_xtal.top.select(f"resSeq {i}") for i in shared_resseqs])
    apo_iis = np.concatenate([apo_xtal.top.select(f"resSeq {i}") for i in shared_resseqs])

    #compute secondary structures
    holo_ss = md.compute_dssp(holo_xtal.atom_slice(holo_iis))[0]
    apo_ss = md.compute_dssp(apo_xtal.atom_slice(apo_iis))[0]

    #compute secondary structure changes where both residues were resolved
    delta_ss = [(apo_ss[x] == i or i == 'NA' or apo_ss[x] == 'NA') for x, i in enumerate(holo_ss)]

    print("+".join([str(i) for i in shared_resseqs[np.where(np.array(delta_ss) == False)]]))

    #compute the number of secondary structure changes between resolved residues
    print(len(holo_ss) - sum(delta_ss))
