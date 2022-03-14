import os
import xml.etree.ElementTree as ET
import mdtraj as md
import numpy as np
import csv
import scipy
import itertools
import re
from operator import itemgetter
import sys

#--------------------------------------------------------------------------------------------
#                           Structure Filtering Methods
#--------------------------------------------------------------------------------------------

#using bowmore nodes
#1.	Enter server with access to nodes: 	ssh borowsky.jonathan@ssh.engr.wustl.edu
#2.	Enter node: 				bsub -q bowman -Is /bin/bash
#to source anaconda on a node with no native conda
#1. source /project/bowmore/borowsky.jonathan/anaconda3/bin/activate
#2. conda activate snakes

cs_only = True
if cs_only:
    directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/new_pockets/rmg_files"
else:
    directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/new_pockets/iofiles"


#make output directories, presently unused
def makedirectories():
    directories = ["moad_xml", "rcsb_pdb_holo", "rcsb_pdb_apo", "monomer_apo", "monomer_holo", "pairing_index"]
    for i in directories:
        os.system(f"mkdir /{i}")

def getligs(struct, chain):

    #WARNING: The use of --no-check-certificate compromises security
    #ligand information (contains validity information not found in pdb structure)
    if struct not in existing_xids_moad:
        os.system(f"wget https://www.bindingmoad.org/files/csv/{struct}.csv -O {directory}/moad_xml/{struct}.csv --no-check-certificate")
        existing_xids_moad.append(struct)
        #keeps the list of existing xml files up to date for when the code encounters apo candidates which are in moad and were previously loaded as holo candidates
        #note that this does not track the extant argument and may require manual adjustment?

    bio_ligands = [] #get all the valid ligands from the moad structure

    chain_elems = chain.split("_")

    with open(f'{directory}/moad_xml/{struct}.csv') as csvfile:
        reader = csv.reader(csvfile)
        firstrow = True
        for row in reader:
            if firstrow: #skip first row, which lacks ligand information
                firstrow = False
                continue

            contents = row[3].split(":")
            if (contents[1] in chain_elems or chain == "") and row[4]=="valid":
            #the latter case (chain == "") is for use of this method in get_struct(),
            #where it is used to obtain a list of all chains containing biologically relevant ligands
                bio_ligands.append(contents)

    return bio_ligands #[keep_ligands, bio_ligands, all_ligands, nst_ligands]


#extract the specified chain(s), complete with heteroatoms, remark 465, and conect records, into a new pdb file
#the code could be made significantly more efficient (but very hard to read)
#by integrating this into the loop in getstruct() (and modifying it to deal with all chains in parallel) so the pdb file is only looped-through once
def extract_chain(struct, intype, outtype, chains):

    reading = False
    read = False #should speed up the code

    buffer = []
    atoms = []

    for line in open(f'{directory}/rcsb_pdb_{intype}/{struct}.pdb'): #read each line in file

        #copy remark 465
        if not read:
            pdblist = line.split()

            #this section was adapted from check_gaps() and could probably be streamlined
            if pdblist == ['REMARK', '465', 'M', 'RES', 'C', 'SSSEQI']: #find end of header and initiate reading
                reading = True
                buffer.append(line)
            elif reading and len(pdblist)==5: #read data
                if pdblist[3] in chains:
                    buffer.append(line) #keep the relevant remark 465 lines

            elif reading and (pdblist[0] != 'REMARK' or pdblist[1] != '465'): #find end of section and stop reading
                if len(buffer) == 1: #if remark 465 exists but none of the missing residues were in the chain being saved here, remove the header
                    #print(buffer)
                    buffer.remove('REMARK 465   M RES C SSSEQI                                                     \n')
                reading = False
                read = True

        #copy atoms, heteroatoms, ter lines, and conect records
        if (line[0:4] == "ATOM" or line[0:6] == "HETATM" or line[0:3] == "TER") and line[21:22] in chains:
            atoms.append(line[6:11]) #used to identify the correct CONECT records
            buffer.append(line) #get the protein atoms in the chain
        elif line[0:6] == "CONECT":
            if line[6:11] in atoms:
            #set(atoms).isdisjoint(line.split(" ")[1:]): #if the CONECT record contains any atoms from the chain of interest;
            #unnecessary because all CONECT records should contain atoms only from one monomer or the other, so looking only at the first element should be sufficient
                buffer.append(line)

    buffer.append("END")

    #save file
    with open(f"/project/bowmanlab/borowsky.jonathan/FAST-cs/new_pockets/rmg_files/monomer_{outtype}/{struct}_chain{'_'.join(chains)}.pdb", 'w') as f:
        f.writelines(buffer)

#get structure without downloading files repeatedly, will check resolution
def getstruct(struct, extant, type):

    #download the pdb structure if necessary
    if struct not in extant:
        try: #download and load .pdb file
            os.system(f"wget https://files.rcsb.org/download/{struct}.pdb -O {directory}/rcsb_pdb_{type}/{struct}.pdb")
            #return md.load(f"rcsb_pdb_{type}/{struct}.pdb")

        except IndexError:
            print("error downloading file; it may be unavailable in .pdb format")
            #large structures cannot be downloaded as pdb files and are too big to be desirable for our present analysis, so I'll skip them instead of writing a .cif parser
            os.system(f"rm {directory}/rcsb_pdb_{type}/{struct}.pdb")
            #remove the empty file generated in the failed download attempt so that subsequent runs don't try to load it in mdtraj and subsequently crash
            return False

    #--------------------------------------------check resolution and break structure into monomeric chains, if any---------------------------------------------------------------------------------------------

    res_checked = False #make sure the structure contains REMARK 2; might be unnecessary depending on whether it's ever missing
    #biomolecule = 0 #counts only monomeric ones
    author_found = False #prevent the software from reading any software-determined-only biomolecules
    #if any author-determined ones are present; assumes that the author-determined opnes come first
    look_for_mer = False #read only the first line after each biomolecule line to avoid reading the
    #software-determined biological unit for a given biomolecule if the author-determined one has already been read
    look_for_chains = False
    rem350 = False
    any_chains = False
    #chains = []
    chain_buffer = [] #chains for a single biomolecule are sometimes spread onto multiple lines

    true_type = []
    false_type = [] #only used for auxilliary holo chains with no biological (valid) ligand

    lig_chains = []
    if type == "holo":
        holo_ligs = getligs(struct, "")
        lig_chains = [l[1] for l in holo_ligs]

    for line in open(f'{directory}/rcsb_pdb_{type}/{struct}.pdb'): #get pdb apo indices

        #print(line)

        #check the resolution, exit if it's inadequate
        if line[0:22] == "REMARK   2 RESOLUTION.":
            #print(line.split(" "))
            pdblist = line.split(" ")
            if pdblist[8] == "NOT" and (pdblist[9] == "AVAILABLE." or pdblist[9] == "APPLICABLE."):
                print("resolution unavailable; skipping structure")
                return False
            else:
                try:
                    if float(pdblist[8]) > 2.5: #if the resolution is inadequate
                        print("inadequate resolution of %s A" % float(pdblist[8]))
                        return False
                    else: #if the resolution is adequate
                        print("resolution of %s A" % float(pdblist[8]))
                        res_checked = True
                except ValueError:
                    print(f"nonstandard REMARK 2: {line}")
                    return False


        #determine if the biological unit is a monomer, look for the chains only if it is
        if look_for_mer and (line[0:45] == "REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT:" or (line[0:52] == "REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE:" and (not author_found))):
            if line[0:45] == "REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT:":
                author_found = True

            asm_type = re.split(": | ", line)[6]
            if asm_type != "MONOMERIC":
                print(f"skipping {asm_type} assembly")
            else:
                print(asm_type)
                look_for_chains = True
                #biomolecule += 1
        elif look_for_mer and not (line[0:52] == "REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE:" and author_found):
            print("missing biological unit; skipping structure")
            return False

        #compile a list of all pdb-chains in the biological unit and then extract them to a new file
        #since non-protein ligands may have their own pdb chains, even monomeric proteins may need to be assembled from multiple chains
        if look_for_chains and line[0:41] == "REMARK 350 APPLY THE FOLLOWING TO CHAINS:":
            pdblist = re.split(',| ', line[41:])
            #about one in ten thousand PDB structures has an extra space before the chain ID so I can't exclude the remark text
            #by splitting on a double space
            for i in range(0,len(pdblist)):
                if pdblist[i] != "" and pdblist[i] != "\n":
                    chain_buffer.append(pdblist[i])
        elif look_for_chains and line[0:19] == "REMARK 350   BIOMT1":
            #chains.append(chain_buffer)
            #extract_chain(struct, type, type, chain_buffer)
            if type == "holo" and set(chain_buffer).isdisjoint(lig_chains): #separate holo monomers with no drug-like ligands and save them as apo structures
                false_type.append('_'.join(chain_buffer))
                if f"monomer_holo/{struct}_chain{'_'.join(chain_buffer)}.pdb" not in existing_hmonomers:
                    extract_chain(struct, type, "apo", chain_buffer)
                print(f"apo chains from holo structures: {chain_buffer}")
            else:                                                           #save apo monomers and holo monomers with drug-like ligands
                true_type.append('_'.join(chain_buffer))
                if f"monomer_apo/{struct}_chain{'_'.join(chain_buffer)}.pdb" not in existing_amonomers:
                    extract_chain(struct, type, type, chain_buffer)
                print(f"regular chains: {chain_buffer}")


            look_for_chains = False
            chain_buffer = [] #reset chain list
            any_chains = True

        look_for_mer = False
        #make sure that only the author determined biological unit is read in cases where software and author units disagree
        #by ensuring that only the line after the biomolecule line is read to check the unit
        if line[0:23] == f"REMARK 350 BIOMOLECULE:": # {biomolecule+1} #did not handle variable numbers of spaces well
            look_for_mer = True
            rem350 = True #check that remark 350 exists

        #exit if remark 350 was missing or if no monomers were found
    if not res_checked:
        print("REMARK 2 (resolution) missing from file; skipping structure")
        return False
    elif not rem350:
        print("no remark 350 found; skipping structure")
        return False
    elif not any_chains:
        print("no monomers found")
        return False

    return [true_type, false_type]


#check that the apo structure is really an apo structure and that its active site is free of crystallographic additives
#[holo crystal structure in mdtraj, holo pdb id, holo chain(s), apo crystal structure in mdtraj, apo pdb id, apo chain(s),]
def get_rmsd(holo_xtal, idh, chh, apo_xtal, ida, cha):

    #-------------------------------------------------------------align apo and holo structures------------------------------------------------

    align_inds = get_indices(ida, cha, idh, chh) #get mdtraj indices of shared residues to use for alignment

    if align_inds[0] == [] or align_inds[1] == []: #exit on encountering mismatched numbering
        print("no indices returned")
        return -1

    #get atom indices in apo and holo structures from the residue indices,
    #using alpha carbons only so the number of atoms matches even if there are missing sidechain atoms
    query_a = np.array([])
    for ind in align_inds[0]:
        query_a = np.concatenate((query_a, apo_xtal.top.select(f"resid {ind} and name CA")))

    query_h = np.array([])
    for ind in align_inds[1]:
        query_h = np.concatenate((query_h, holo_xtal.top.select(f"resid {ind} and name CA")))

    #The 'if' case might not be redundant with the check of align_inds above in the ridiculous corner case where
    #the numbering is shifted such that there's a single residue of overlap but that residue is missing an alpha carbon
    if len(query_a) == 0 or len(query_h) == 0:
        print("numbering mismatch in pdb files; may be useable if renumbered")
        return -1
    elif len(query_a) != len(query_h):
        print("numbering mismatch in mdtraj, probably arising from a missing internal alpha carbon")
        return -1
        #note that even with this, there is a risk that both proteins will be missing alpha carbons in different places,
        #causing the code to proceed silently and return an unreasonably high RMSD and a bad alignment.
        #This ought to be fixed properly in the long term by determining the conditions under which mdtraj will load residues and adjusting indices accordingly

    #the above issue (in the elif case) occurs as follows in something like 1 in 250 structures:
    #1. a pdb structure is missing an internal alpha carbon
    #2. get_indices sees this and puts a gap in the pdb index array
    #3. When converting from pdb to mdtraj indices,
    #   get_indices() assumes that mdtraj will not load this residue either when it loads the structure.
    #   It therefore assumes that the extant residues on either side of the gap will be numbered consecutively,
    #   and does not include a gap in the mdtraj index array.
    #4. mdtraj loads the CA-less internal residue anyway. Thanks mdtraj.
    #5. Since the mdtraj index array has no gap, when this block of code in checkligs() runs,
    #   the consecutively numbered mdtraj index array causes it to try to load the CA for the CA-less residue.
    #   This loads nothing silently but leaves one atom number array one element too short
    #6. superpose crashes

    #align apo and holo structures using alpha carbons
    #consider aligning using only ligand-coordinating residues for structures with large conformational changes
    #apo_xtal.superpose(holo_xtal, atom_indices = query_a.astype(int), ref_atom_indices = query_h.astype(int)) #modifies apo_xtal in place

    apo_calphas = apo_xtal.atom_slice(query_a.astype(int))
    holo_calphas = holo_xtal.atom_slice(query_h.astype(int))
    pair_rmsd = md.rmsd(apo_calphas, holo_calphas)[0] #does not appear to read atom_indices arguments

    #apo_xtal.save_pdb("superpose_test.pdb")
    return pair_rmsd


#identify residues found in both apo and holo structures and determine their indices in mdtraj for use in alignment
def get_indices(apo_id, apo_chain, holo_id, holo_chain):
    #pdb indices
    apo_indices = []
    holo_indices = []
    apo_resns = {}
    holo_resns = {}
    #shared_indices = []
    #md_loadable_apo_indices = []
    #md_loadable_holo_indices = []

    #mdtraj indices
    apo_md_indices = []
    holo_md_indices = []

    for line in open(f'{directory}/monomer_apo/{apo_id}_chain{apo_chain}.pdb'): #get pdb apo indices
        if line[0:5] == "ATOM " and line[13:16] == "CA ":
            apo_indices.append(int(line[22:26]))
            apo_resns[int(line[22:26])] = line[17:20]

    for line in open(f'{directory}/monomer_holo/{holo_id}_chain{holo_chain}.pdb'): #get pdb holo indices
        if line[0:5] == "ATOM " and line[13:16] == "CA ":
            holo_indices.append(int(line[22:26]))
            holo_resns[int(line[22:26])] = line[17:20]

    apo_indices = np.unique(apo_indices).tolist()
    holo_indices = np.unique(holo_indices).tolist()

    #convert pdb indices to mdtraj indices
    for i in range(max(apo_indices[0], holo_indices[0]), min(apo_indices[-1], holo_indices[-1])+1):
        if i in apo_indices and i in holo_indices:
            if holo_resns[i] != apo_resns[i]:
                print("numbering mismatch in pdb files; may be useable if renumbered")
                return [[],[]]
            apo_md_indices.append(apo_indices.index(i))
            holo_md_indices.append(holo_indices.index(i))

    return [apo_md_indices, holo_md_indices]

#look for gaps in the protein structure exceeding three residues and find the residues at the edges of gaps to check their distance to the holo candidate ligand
def gap_length(type, pdb_id, chain):
    reading = False #skip header
    #gap_ends = [] #residue ids of residues on either side of each internal gap
    #lastid = 0
    total_gaps = 0

    #get the indices of the first residue to distinguish between missing termini and internal gaps
    first_resi = -10000 #will fail if any pdb file actually has a residue -10000
    last_resi = 0

    for line in open(f'{directory}/monomer_{type}/{pdb_id}_chain{chain}.pdb'): #find first residue
        pdbresi = line[22:26] #the pdb residue number
        pdbhead = line[0:4] #the first line element
        if pdbhead == "ATOM" and first_resi == -10000:
            first_resi = int(pdbresi) #will continue operating without warning if alphabetic residue insertion codes are encountered; see below
        if pdbhead == "TER ":
            last_resi = int(pdbresi) #will continue operating without warning if alphabetic residue insertion codes are encountered; see below
            break
        #lastline = pdblist

    for line in open(f'{directory}/monomer_{type}/{pdb_id}_chain{chain}.pdb'): #track gap length
        pdblist = line.split()
        #print(pdblist)
        if reading and len(pdblist)==5: #read data
            #detect cases of residue indices with multiple different residue types (not just different conformations) and skip them
            try:
                res_id = int(pdblist[4])
            except ValueError:
                print("residue with [alphabetic] insertion code; skipping structure")
                return -1

            if res_id > first_resi and res_id < last_resi: #don't read missing terminal residues
                total_gaps =+ 1

        elif pdblist == ['REMARK', '465', 'M', 'RES', 'C', 'SSSEQI']: #find end of header and initiate reading
            reading = True
        elif reading and (pdblist[0] != 'REMARK' or pdblist[1] != '465'): #stop reading and return the gap-lining residues
            return [last_resi-first_resi, total_gaps] #total length and gap length

    #print("no REMARK 465 present") #this occurs for proteins with no missing residues
    return [last_resi-first_resi,0] #no gaps

#--------------------------------------------------------------------------------------------
#                           Global Variables and Main Loop
#--------------------------------------------------------------------------------------------

include_resis = True #whether to include holo ligand coordinating residues in determining which apo ligands are too close
#this could be achieved by setting coord_threshold to 0 but that would be computationally inefficient compared with never
#calculating the distances in the first place
dist_threshold = 0.5 #[nanometers] minimum acceptable separation in nm between apo and holo ligands for the apo ligands to be considered outside of the holo binding site
coord_threshold = 0.5 #[nanometers] all residues with atoms within this distance of the holo ligand in the holo structure are considered ligand-coordinating
pept_threshold = 1   #any peptide longer than this is considered a protein chain, any peptide as long or shorter is considered a ligand

existing_xids_moad = [i[0:4] for i in os.listdir(f"{directory}/moad_xml/")] #get a list of already-downloaded ligand lists to avoid downloading extra copies
#existing_xids = [i[0:4] for i in os.listdir("./rcsb_xml/")] #get a list of already-downloaded validation reports to avoid downloading extra copies
existing_hids = [i[0:4] for i in os.listdir(f"{directory}/rcsb_pdb_holo/")] #get a list of already-downloaded candidate holo structures to avoid downloading extra copies
existing_aids = [i[0:4] for i in os.listdir(f"{directory}/rcsb_pdb_apo/")] #get a list of already-downloaded candidate apo structures to avoid downloading extra copies

existing_amonomers = [i for i in os.listdir(f"/project/bowmanlab/borowsky.jonathan/FAST-cs/new_pockets/rmg_files/monomer_apo/")]
existing_hmonomers = [i for i in os.listdir(f"/project/bowmanlab/borowsky.jonathan/FAST-cs/new_pockets/rmg_files/monomer_holo/")]
#existing_saved = os.listdir("./paired_pdb_structures") #unimplemented; the save function is unused

blast_out = os.listdir(f"{directory}/blast_output/") #filenames for blastp results of candidate holo (MOAD) structures with at least one candidate apo (PDB blastp hit)structure

#pair_ligands = [] #all (if any) suitable ligands found for each apo-holo pair
#lowrmsd_pairs = [] #the lowest rmsd (if any) suitable ligands found for each apo-holo pair; used to filter out pockets which open spontaneously

pairs_gaps_rmsd = []

#compares all sub_chains of two chains against each other, or compares liganded and ligandless 'holo' chains
#note that this method relies on side effects and global variables declared immediately above
    #for the first call:  holo_chain_list = holo_chains, apo_chain_list = apo_chains[0],                holo_sid=holo_id, apo_sid=apo_id
    #for the second call: holo_chain_list = holo_chains, apo_chain_list = holo_chains_unfiltered[1],    holo_sid=holo_id, apo_sid=holo_id
def compare_chains(holo_chain_list, apo_chain_list, holo_sid, apo_sid, ha_type, h_gaps):
    print(holo_chain_list)
    for a_chain in apo_chain_list: #for all the apo candidate bioassemblies
        print(f"{ha_type}apo_chain {a_chain}")
        #if check_gaps("apo", apo_sid, a_chain): #check for gaps
        a_gaps = gap_length("apo", apo_sid, a_chain)
        apo = md.load(f"{directory}/monomer_apo/{apo_sid}_chain{a_chain}.pdb")

        for h_chain in holo_chain_list: #for all the holo candidate bioassemblies with ligands and no gaps
            print(f"holo_chain {h_chain}")
                #gaps checked above
            holo = md.load(f"{directory}/monomer_holo/{holo_sid}_chain{h_chain}.pdb")
            rmsd = get_rmsd(holo, holo_sid, h_chain, apo, apo_sid, a_chain)

            if get_rmsd != -1:
                #print(h_gaps)
                #print(a_gaps)
                mean_gap_fraction = (h_gaps[h_chain][1]/h_gaps[h_chain][0] + a_gaps[1]/a_gaps[0])/2
                print([mean_gap_fraction, rmsd])
                #print([holo_sid, h_chain, apo_sid, a_chain, h_gaps[h_chain], a_gaps, mean_gap_fraction, rmsd])
                pairs_gaps_rmsd.append([holo_sid, h_chain, apo_sid, a_chain, h_gaps[h_chain], a_gaps, mean_gap_fraction, rmsd])


#-------------------------------------------------------------Main Loop------------------------------------------------
#--------------------------------------------------loop through holo structures----------------------------------------

debug = True

if debug:
    rstart = 0
    rend = 94
else:
    rstart = 0
    rend = len(blast_out)

#20272 and 20273 (holo 5lwr and 5lws) crash on apo structures 5qb7 and 5oyv respectively
#for reasons which are probably too rare to be worth debugging

#savenums = range(0,50000,1000) #intervals at which to save output so far

all_cs_structures = [["2CGA","B","1AFQ","C","0FG"],["4AKE","B","1ANK","B","ANP"],["1RTC","A","1BR6","A","PT1"],["1ADE","A","1CIB","A","IMP"],["1HOO","B","1CIB","A","HDA"],["1CLL","A","1CTR","A","TFP"],["2WGQ","B","1D6Y","B","HY1"],["1ECJ","D","1ECC","B","PCP"],["1DUB","D","1EY3","F","DAK"],["1NUW","A","1EYJ","B","AMP"],["3PUW","E","1FQC","A","GLO"],["1MY1","C","1FTL","A","DNQ"],["1G4E","B","1G67","B","POP/TZP"],["1HAG","E","1GHY","H","121"],["1EX6","A","1GKY","A","5GP"],["1KZ7","D","1GRN","A","AF3"],["1BSQ","A","1GX8","A","RTL"],["1G24","D","1GZF","C","NIR"],["1E2X","A","1H9G","A","MYR"],["1EXM","A","1HA3","B","MAU"],["1IMF","A","1IMB","B","LIP"],["1RDW","X","1J6Z","A","RHO"],["1B6B","A","1KUV","A","CA5"],["1FA9","A","1L5S","B","URC"],["1ALB","A","1LIC","A","HDS"],["1MY0","B","1N0T","D","AT1"],["1ALV","A","1NX3","A","ISA"],["1OK8","A","1OKE","B","BOG"],["1XCG","B","1OW3","B","GDP"],["1Z92","A","1PY2","A","FRH"],["1JWP","A","1PZO","A","CBT"],["1PZT","A","1PZY","D","UDP"],["3HQD","A","1Q0B","B","NAT"],["1BP5","A","1RYO","A","OXL"],["1RRG","A","1S9D","A","AFB"],["2GPO","A","1S9Q","B","CHD"],["2F6V","A","1T49","A","892"],["1TQO","A","1TR5","A","THP"],["1K3F","B","1U1D","F","181"],["1JBU","H","1WUN","H","P5B"],["1XMG","B","1XVC","A","5BR"],["2AKA","A","1YV3","A","BIT"],["2AIR","H","1ZA1","D","CTP"],["3B7D","E","2AL4","F","CX6"],["3CJ0","A","2BRL","A","POO"],["3CJ0","A","3FQK","B","79Z"],["2BU8","A","2BU2","A","TF1"],["3PEO","G","2BYS","J","LOB"],["1K5H","C","2EGH","B","FOM"],["1SWX","A","2EUM","A","LAT"],["2BRK","A","2GIR","B","NN3"],["1UK2","A","2GZ7","A","D3F"],["2CM2","A","2H4K","A","509"],["1NEP","A","2HKA","C","C3S"],["3L7U","C","2HVD","C","ADP"],["1A8I","A","2IEG","B","FRY"],["3CHE","A","2IUZ","B","D1H"],["1H09","A","2IXU","A","MU2"],["2IYT","A","2IYQ","A","ADP/SKM"],["2GFC","A","2JDS","A","L20"],["2ZB1","A","2NPQ","A","BOG"],["1KS9","A","2OFP","A","PAF"],["2OHG","A","2OHV","A","NHL"],["1FVR","A","2OO8","X","RAJ"],["1ZAH","B","2OT1","D","N3P"],["2AX9","A","2PIQ","A","RB1"],["2Q8F","A","2Q8H","A","TF4"],["2WGB","A","2V57","A","PRL"],["1BNC","B","2V5A","A","LZL"],["1RHB","A","2W5K","B","NDP"],["3GXD","B","2WCG","A","MT5"],["2QFO","B","2WI7","A","2KL"],["1QLW","B","2WKW","B","W22"],["2YQC","A","2YQS","A","UD1"],["3FDL","A","2YXJ","A","N3C"],["3BL9","B","3BL7","A","DD1"],["3F74","C","3BQM","C","BQM"],["2H4E","B","3CFN","B","2AN"],["2QLR","C","3DC1","A","AKG"],["2BF3","A","3DHH","E","BML"],["3MN9","A","3EKS","A","CY9"],["1R1W","A","3F82","A","353"],["1SU4","A","3FGO","B","ACP"],["2BLS","B","3GQZ","A","GF7"],["3H5R","A","3H9J","D","APC"],["3NNU","A","3HL7","A","I46"],["1FXX","A","3HL8","A","BBP"],["1NI6","D","3HOK","B","Q80"],["1PKL","B","3HQP","P","ATP/FDP/OXL"],["3DXN","A","3HZT","A","J60"],["1HKA","A","3IP0","A","HHS"],["1W50","A","3IXJ","C","586"],["3KQA","B","3LTH","A","UD1"],["4HB2","C","4HAT","C","LMB"]]

for holo_candidate in range(rstart,rend): #loop through candidate holo structures

    if not cs_only:
        holo_id = blast_out[holo_candidate][0:4].lower() #get id at index
        apo_ids = [id.lower() for id in np.load(f"{directory}/blast_output/{holo_id.upper()}_hits.npy")] #load blast results
    else:
        holo_id = all_cs_structures[holo_candidate][2].lower()
        apo_ids = ["", all_cs_structures[holo_candidate][0].lower()]
        #print(holo_id)
        #print(apo_ids)

    print(f"------------------------------------------------------------------------------------------------------------")
    print(holo_candidate)
    print(f"holo: {holo_id}")
    #break pdb structures from rcsb into their constituent physiological assemblies and check resolution and oligomerization state
    holo_chains_unfiltered = getstruct(holo_id, existing_hids, "holo")
    if isinstance(holo_chains_unfiltered, list): #avoid crashing when a structure doesn't load

        if len(apo_ids[1:]) == 0 and len(holo_chains_unfiltered[1]) == 0:
            print("no candidate apo structures available")
            continue

        #filter structures with unacceptable gaps out of the list of holo candidate bioassemblies with ligands once rather than for each apo structure
        #holo_ligands = {}
        holo_lig_iis = {} #coordinates of holo ligands and ligand coordinating residues to avoid recomputing them
        holo_chains = [] #holo chains which actually have ligands and lack long gaps
        holo_gaps = {}
        #apo_gaps = {}

        #print(holo_chains_unfiltered)
        for h_chain in holo_chains_unfiltered[0]: #index of 0 is for true holo chains
            print(f"checking gaps in holo chain {h_chain}")
            #This may be the last thing printed for holo chains with no apo or ligandless holo chains;
            #I don't want to waste time checking explicitly for this case since exiting silently is fine
            #print(h_chain)
            holo_gaps[h_chain] = gap_length("holo", holo_id, h_chain)
            #    holo_chains.append(h_chain)

        #if holo_chains != []: #continue only if there are suitable holo chains

            #buffer_pair_ligands = [] #used to collect only the lowest RMSD apo-holo pair from each sequence to approximate Sun et. al.'s true pocket screen

            #--------------------------------------------------loop through apo structures-----------------------------------
        #print(holo_chains_unfiltered[0])
        for apo_id in apo_ids[1:]: #loop through candidate apo structures [=blast results] (1st element is candidate holo id)
            print("------------------------------------------------------")
            print(f"apo: {apo_id}")
                #break pdb structures from rcsb into their constituent physiological assemblies and check resolution and oligomerization state
            apo_chains = getstruct(apo_id, existing_aids, "apo")
            if isinstance(apo_chains, list):

                compare_chains(holo_chains_unfiltered[0], apo_chains[0], holo_id, apo_id, "", holo_gaps)

        compare_chains(holo_chains_unfiltered[0], holo_chains_unfiltered[1], holo_id, holo_id, "h", holo_gaps)

            #sort and save buffer_pair_ligands for analysis in case the script crashes before terminating
            #add the lowest rmsd pair to the list thereof
            #if buffer_pair_ligands != []:
            #    buffer_pair_ligands = sorted(buffer_pair_ligands, key = itemgetter(4), reverse = False)
            #    np.save(f"{directory}/pairing_index/ligand_pairs_{holo_candidate}", np.array(buffer_pair_ligands, dtype = object)) #save the ranked list of apo-holo pairs
            #    lowrmsd_pairs.append(buffer_pair_ligands[0]) #add the lowest rmsd pair

#    if holo_candidate in savenums: #save output so far
#        np.save(f"{directory}/pairing_index/all_ligand_pairs_{rstart}_{holo_candidate}", np.array(pair_ligands, dtype = object))
#        np.save(f"{directory}/pairing_index/all_lowrmsd_ligand_pairs_{rstart}_{holo_candidate}", np.array(lowrmsd_pairs, dtype = object))

#TODO: ought to sort these
#save lists of all apo holo pairs and low rmsd pairs only respectively
np.save(f"/project/bowmanlab/borowsky.jonathan/FAST-cs/new_pockets/rmg_files/pairing_index/rmsd_vs_gaps_all", pairs_gaps_rmsd)#np.array(pair_ligands, dtype = object))
#lowrmsd_pairs = sorted(lowrmsd_pairs, key = itemgetter(4), reverse = True)
#np.save(f"{directory}/pairing_index/all_lowrmsd_ligand_pairs_{rstart}_{rend}", np.array(lowrmsd_pairs, dtype = object))

#print("")
print(pairs_gaps_rmsd)
#print(sorted(rmsd_ranks, key = itemgetter(0), reverse = True)) #this may be redundant with pair_ligands


#TODO:
#look at the pair at holo indicex 112
#DONE: first thing: fix blast processing; see 2ae7

#general
#-1-. check gap-end-residue - holo site distances #on hold since it will probably have to be rewritten or modified significantly for multiple chains
#--> 3. fix multimer processing: compare monomer chains against each other where multiple ones are available
    #read resolution from pdb file while you're at it to avoid loading the xml

#Plan:
#[called as necessary for each structure]:
#1. (down)load the pdb file
#1.5. check resolution and proceed only if it's adequate
#2. read remark 350:
    #if remark 350 is absent, report this to the user and skip structure
    #check the seqres records; proceed only if all sequences are the same; report heteromultimer otherwise
    #   I assume that since multiple proteins are only added to crystallization solutions for the purpose of obtaining complexes thereof,
    #   all monomeric units will be from structures with only one protein sequence
    #   (in practice there could be monomers of different sequences if the complex failed to form or had an unexpected stoichiometry)
    #for each biomolecule in remark 350:
        #if remark 350 is present, use the author determined] biological unit if it's available and the software determined biological unit otherwise
        #if the biological unit is 'monomeric', proceed. Otherwise continue to the next biomolecule.
        #for holo structures, note which biomolecules have chains containing valid ligands; use those that do not as supplementary holo structures
        #   This could be also done for moad structures with no blast hits, which are currently excluded from the analysis, to find more apo-holo pairs
        #   3hok is a good example of this
        #   see below for a kludgier but more convenient alternative
        #save remark 465, atoms, and heteroatoms for all chains in the biomolecule to a (single) new pdb file
        #   (there may be multiple chains per monomeric biomolecule because ligands may have their own chains)

#3. feed the pdb files generated above into the current code in place of an apo or holo file;
    #compare every holo chain against every apo chain for each blast hit
    #i.e. (the order of the middle two loops could probably be reversed, but the current order may be better for allowing the structures to be unpacked into chains only as necessary)
    #for holo:
    #   for monomer chain in holo with moad ligand: [if chain is okay (e.g. in terms of gaps)]
    #       for apo:
    #           for monomer chain in apo: [if chain is okay]
    #               determine whether [candidate] apo chain is really apo, save true pairs
    #       for monomer chain in holo without moad ligand: [if chain is okay]
    #            determine whether pseudo - [candidate] holo chain is really apo, save true pairs
    #       Note that it would be programmatically simpler but computationally slightly less efficient to just...
    #       ...include the holo structure as an apo candidate alongside the regular apo candidates


    #
#
#
#

#crypticness
#-2-. make sure cavities are real: holo volume of at least 10 with ligsite (or use fpocket)
#on hold since all 3/100 current apo-holo pairs have no motion and we can't test this properly until we have more
    #check what fraction of structures are eliminated here --> frequency of cryptic structures vs other binding sites
#4. filter for or rank by larger conformational changes
    #may need to implement better apo-holo alignment

#ligand filtering
#6. check for organic ligands?
#5. filter very small ligands
#9. consider allowing less restrictive ligand-proximity criteria since many holo ligands bind next to other ligands,
#which might be found in both apo and holo structures.
#Of course for ligand-free simulations where we want to track the pocket lining residues,
#having a pocket wall formed by a non-protein ligand might be a problem instead

#misc
#8. clean code and add more comments



#-----------------------------------------------------------------------CODE SCRAPS BELOW--------------------------------------------------------------

foo="""if h_chain not in holo_lig_iis:
        ligands_iis = checkligs(holo, holo_sid, h_chain, apo, apo_sid, a_chain, {})
        if ligands_iis != []:
            holo_lig_iis[h_chain] = ligands_iis[1]
            #update the coordinate dictionary if this is a new holo chain and holo ligand and residue coordinates were obtained
            #avoids errors when upstream errors mean that no coordinates are obtained
        else:
            continue #avoid index error below if there's an issue with the structures
    else:
        ligands_iis = checkligs(holo, holo_sid, h_chain, apo, apo_sid, a_chain, holo_lig_iis[h_chain])
        if ligands_iis == []:
            continue #avoid index error below if there's an issue with the structures

    ligands = ligands_iis[0]

    #save protein structures and pairing information
    if ligands != []:

        #save proteins, recording any anomalous ones
        #saveprot(apo, apo_sid, a_chain, "apo")
        #saveprot(holo, holo_sid, h_chain, "holo")

        #np.save(f"pairing_index/ligand_pairs_{holo_candidate}", np.array([holo_sid, h_chain, apo_sid, a_chain, ligands_iis[2][0], ligands], dtype = object))
        pair_ligands.append([holo_sid, h_chain, apo_sid, a_chain, ligands_iis[2][0], ligands])
        buffer_pair_ligands.append([holo_sid, h_chain, apo_sid, a_chain, ligands_iis[2][0], ligands])
        #rmsd_ranks.append([ligands_iis[2][0], holo_sid, h_chain, apo_sid, a_chain] + np.array(ligands).flatten().tolist())

        print(f"pair found; rmsd {ligands_iis[2][0]}")
    else:
        print("no true apo structure found")
"""
