import numpy as np
import mdtraj as md
import itertools
import scipy
import os
import sys

serial_in = "2b" #which version of the pocket set to load
serial_out = "3b-nr" #which version of the pocket set to save
input_directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/new_pockets_2" #also used for some outputs

distance_threshold = 0.5

all_pockets = np.load(f"{input_directory}/filter_output/all_lowrmsd_ligand_pairs_v{serial_in}.npy")

#notes
#6rsk has messed up numbering that skips 0 and starts with a negative numbered residues

#select residue backbones and sidechains which coordinate the ligand using heavy atoms
def select_resi_parts(ligand_resn, holo_xtal, apo_xtal, distance_threshold, ind_debug):

    #--------------------ligand selection string formatting---------------------

    individual_resns_all = ["'"+i[0].split(" ")[0]+"'" for i in ligand_resn]

    #note that one ligand goes by 986 in rcsb, 098 in moad, and 98 in mdtraj
    #L9K and L9W are actually two different tautomers of the same compound, but rcsb only includes L9W, presumably because the model has no hydrogens so the tautomer is not resolved
    #see http://polymorph.sgc.utoronto.ca/drugged_human_proteome/pages/col_12_IPR001680.html for more information
    #MOAD's TNR is broken into two residues in the modern pdb; the TNR residue is listed as obsolete in the european PDB

    # [MOAD residue name]:[mdtraj residue name]
    mdtraj_moad_ligand_mismatches = {"ADE":"A","URA":"U","BAM":"BEN","DCY":"DC","098":"98","L9K":"L9W","TNR":"A2G' or (resname SER and resSeq 906) or resname 'Null"}
    #note that "or resname 'Null" exists to make the single quotes match up properly when the ligand select string is assembled without having to modify the code below

    for x, i in enumerate(individual_resns_all):
        if i[1:-1] in mdtraj_moad_ligand_mismatches.keys():
            individual_resns_all[x] = mdtraj_moad_ligand_mismatches[i[1:-1]]

    #print(individual_resns_all)
    ligand_select_str = "resname "+" or resname ".join(individual_resns_all)
    #print(ligand_select_str)

    #---------------------------------------------------------------------------

    #note that heavy atom selection does not account for deuterated proteins or ligands
    #select ligand heavy atoms
    heavy_lig = holo_xtal.top.select(f"({ligand_select_str}) and not element H")

    #select protein heavy atoms
    heavy_prot = holo_xtal.top.select("protein and not element H")

    #compute protein-ligand distances
    lig_prot_pair_iis = np.array(list(itertools.product(heavy_lig, heavy_prot)))
    prot_lig_dists = md.compute_distances(holo_xtal,lig_prot_pair_iis, periodic=False).flatten()

    #select protein atoms within distance threshold
    lig_coord_iis = np.where(prot_lig_dists < distance_threshold)[0]
    prot_iis = np.unique([lig_prot_pair_iis[i][-1] for i in lig_coord_iis])

    #obtain all backbone and sidechain atom indices
    sidechain_iis = holo_xtal.top.select("sidechain")
    bb_iis = holo_xtal.top.select("backbone")

    #select backbones and sidechains containing the ligand-coordinating atoms
    sele = []
    for i in prot_iis:

        resi = holo_xtal.top.atom(i).residue.resSeq #get pdb residue number
        if resi < 0:
            i_odd.append(ind_debug)
            print("negative residue number encountered; skipping because mdtraj can't handle it")
            continue

        #separate sidechains and backbones
        if i in sidechain_iis:
            sele.append("sidechain and resSeq %s" % str(resi))
        elif i in bb_iis:
            sele.append("backbone and resSeq %s" % str(resi))
        else:
            print(f'error: atom {i} is in neither sidechain nor backbone')
            break

    sele = np.unique(sele)

    prot_iis_holo = []
    prot_iis_apo = []

    prot_iis_holo_matching = []
    prot_iis_apo_matching = []

    for sel in sele:
        #include only indices of atoms in residues present in apo and holo structures
        if len(holo_xtal.top.select(f"{sel} and not element H")) > 0 and len(apo_xtal.top.select(f"{sel} and not element H")) > 0:
            prot_iis_holo.append(holo_xtal.top.select(f"{sel} and not element H"))
            prot_iis_apo.append(apo_xtal.top.select(f"{sel} and not element H"))

            #filter out incomplete residues with different numbers of atoms resolved for RMSD calculations
            #this is not robust to different atom numberings but this should be standardized
            if len(holo_xtal.top.select(f"{sel} and not element H")) == len(apo_xtal.top.select(f"{sel} and not element H")):
                prot_iis_holo_matching.append(holo_xtal.top.select(f"{sel} and not element H"))
                prot_iis_apo_matching.append(apo_xtal.top.select(f"{sel} and not element H"))

    #catch cases where apo is highly truncated and there are no pocket residues resolved
    #or where there are no well-resolved pocket residues to calculate RMSD from
    try:
        prot_iis_holo = np.concatenate(prot_iis_holo).ravel()
        prot_iis_apo = np.concatenate(prot_iis_apo).ravel()
        prot_iis_holo_matching = np.concatenate(prot_iis_holo_matching).ravel()
        prot_iis_apo_matching = np.concatenate(prot_iis_apo_matching).ravel()
    except ValueError:
        return False

    holo_lining = holo_xtal.atom_slice(prot_iis_holo_matching.astype(int)) #atoms lining the cryptic pocket
    apo_lining = apo_xtal.atom_slice(prot_iis_apo_matching.astype(int)) #atoms lining the cryptic pocket

    cs_rmsd = md.rmsd(apo_lining, holo_lining) #calculate the cryptic site rmsd

    return [sele, prot_iis_holo, prot_iis_apo, cs_rmsd]

#get average apo and holo active site atom distances
def get_lining_distances(holo_xtal, apo_xtal, holo_iis, apo_iis):

    holo_lining_coords = holo_xtal.atom_slice(holo_iis).xyz[0]

    holo_size = np.mean(scipy.spatial.distance.cdist(holo_lining_coords, holo_lining_coords))
    #print(holo_lining_coords.xyz[0])
    apo_lining_coords = apo_xtal.atom_slice(apo_iis).xyz[0]
    #print(apo_lining_coords.xyz[0])
    apo_size = np.mean(scipy.spatial.distance.cdist(apo_lining_coords, apo_lining_coords))

    if apo_size > holo_size: #reverse pockets
        return 0
    elif holo_size > apo_size: #forward pockets
        return 1
    else:
        print("pockets are the same; probable filtering error")
        return 2

#---------------------------main filtering--------------------------------------

#get a list or previously filtered pockets
already_filtered = np.load("/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/new_pockets/iofiles/output_indices/some_lowrmsd_ligand_pairs_v6.npy")
already_filtered_holos = [i[0] for i in already_filtered]


x_forward = 0
x_reverse = 0
i_odd = []

forward_pockets = []
reverse_pockets = []

#assume that apo and holo structures are both in the same iofiles_directory
iofiles1 = "iofiles"
iofiles2 = "iofiles_na"

holo_list1 = os.listdir(f"{input_directory}/{iofiles1}/monomer_holo")
holo_list2 = os.listdir(f"{input_directory}/{iofiles2}/monomer_holo")

for testind in range(len(all_pockets)):

    print(f"{testind}: {all_pockets[testind]}")

    #print(all_pockets[testind][0])
    #print(already_filtered_holos)
    if all_pockets[testind][0] in already_filtered_holos:
        print("skipped")
        continue

    if f"{all_pockets[testind][0]}_chain{all_pockets[testind][1]}.pdb" in holo_list1:
        iofiles_switch = iofiles1
    else:
        iofiles_switch = iofiles2

    holo_xtal = md.load(f"{input_directory}/{iofiles_switch}/monomer_holo/{all_pockets[testind][0]}_chain{all_pockets[testind][1]}.pdb")
    apo_xtal = md.load(f"{input_directory}/{iofiles_switch}/monomer_apo/{all_pockets[testind][2]}_chain{all_pockets[testind][3]}.pdb")

    sele_iis = select_resi_parts(all_pockets[testind][5], holo_xtal, apo_xtal, distance_threshold, testind)
    if sele_iis == False:
        print("no suitable pocket residues were resolved, probably due to a truncated apo structure or very poorly resolved cryptic site")
        i_odd.append(testind)
        continue

    lining_result = get_lining_distances(holo_xtal, apo_xtal, sele_iis[1], sele_iis[2])

    if lining_result == 1:
        forward_pockets.append(np.append(all_pockets[testind], sele_iis[3]))
        x_forward += 1
    elif lining_result == 0:
        reverse_pockets.append(np.append(all_pockets[testind], sele_iis[3]))
        x_reverse += 1
    elif lining_result == 2:
        i_odd.append(testind)
    else:
        print("error; invalid lining comparison result")

print(reverse_pockets)

np.save(f"{input_directory}/filter_output/forward_lowrmsd_ligand_pairs_v{serial_out}.npy", forward_pockets)
np.save(f"{input_directory}/filter_output/reverse_lowrmsd_ligand_pairs_v{serial_out}.npy", reverse_pockets)

print("---------------------------------------------------------------------------------------")
print(f"{len(all_pockets)} pockets")
print(f"{x_forward} forward pockets")
print(f"{x_reverse} reverse pockets")

print(f"indices of odd pockets: {np.unique(i_odd)}")


#indices of odd pockets: [  12   79   80   99  147  193  218  276  401  422  499  540  722  777 1008 1071 1193 1277 1332 1499 1555]

#1672 pockets
#347 forward pockets
#513 reverse pockets
#indices of odd pockets: [  12   79   99  147  193  218  276  401  499  540  722  777 1008 1071 1277 1332 1499]


#---------------------------------------------------------------------------TRIMMINGS------------------------------------------------------------------------

    #sanity check:
    #    a = np.array([-0.9117, 0.5899, -1.5034])
    #    b = np.array([-0.8598, 0.5826, -1.6403])
    #    print(np.sqrt(np.dot(a-b,a-b)))
    #works

    #note that this code does not introduce any duplicate indices and hence no application of np.unique is needed



    #more consise but does not address check that residues are present in both apo and holo
    #prot_iis_holo = np.concatenate([holo_xtal.top.select(f"{sel} and not element H") for sel in sele]).ravel()
    #prot_iis_apo = np.concatenate([apo_xtal.top.select(f"{sel} and not element H") for sel in sele]).ravel()

    #return [sele, prot_iis_holo, prot_iis_apo]
