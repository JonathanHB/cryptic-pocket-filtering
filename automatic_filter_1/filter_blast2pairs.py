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
import scipy.spatial.distance

from time import perf_counter

#--------------------------------------------------------------------------------------------
#                           Structure Filtering Methods
#--------------------------------------------------------------------------------------------

#using bowmore nodes
#1.	Enter server with access to nodes: 	ssh borowsky.jonathan@ssh.engr.wustl.edu
#2.	Enter node: 				bsub -q bowman -Is /bin/bash
#to source anaconda on a node with no native conda
#1. source /project/bowmore/borowsky.jonathan/anaconda3/bin/activate
#2. conda activate snakes

blast_directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/new_pockets/iofiles"
directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/new_pockets_2/iofiles_na"

#make output directories, presently unused; ought to be done in bash script since this is run multiple times
def makedirectories():
    os.system(f"cd {directory}")
    directories = ["moad_xml", "rcsb_pdb_holo", "rcsb_pdb_apo", "monomer_apo", "monomer_holo", "pairing_index", "warning_files", "note_files", "qc_structures"]
    for i in directories:
        os.system(f"mkdir {directory}/{i}")

def log_prot_rejection(reject_reason):
    with open(f"{directory}/note_files/protein_rejection_reasons.log", 'a') as f:
        f.writelines(f"protein_rejected: {reject_reason}\n")

#extract the specified chain(s), complete with heteroatoms, remark 465, and conect records, into a new pdb file
#the code could be made significantly more efficient (but very hard to read)
#by integrating this into the loop in getstruct() (and modifying it to deal with all chains in parallel) so the pdb file is only looped-through once
def extract_chain(struct, intype, outtype, chains):

    reading = False
    read = False #should speed up the code

    buffer = []
    atoms = []

    for line in open(f'{blast_directory}/rcsb_pdb_{intype}/{struct}.pdb'): #read each line in file

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
    with open(f"{directory}/monomer_{outtype}/{struct}_chain{'_'.join(chains)}.pdb", 'w') as f:
        f.writelines(buffer)

#get structure without downloading files repeatedly, will check resolution
def getstruct(struct, extant, type):

    #download the pdb structure if necessary
    if struct not in extant:
        try: #download and load .pdb file
            os.system(f"wget https://files.rcsb.org/download/{struct}.pdb -O {blast_directory}/rcsb_pdb_{type}/{struct}.pdb")

        except IndexError:
            print("error downloading file; it may be unavailable in .pdb format")

            log_prot_rejection(f"{struct} not available from rcsb")
            #large structures cannot be downloaded as pdb files and are too big to be desirable for our present analysis
            os.system(f"rm {blast_directory}/rcsb_pdb_{type}/{struct}.pdb")
            #remove the empty file generated in the failed download attempt so that subsequent runs don't try to load it in mdtraj and subsequently crash
            return False

    #--------------------------------------------check resolution and break structure into monomeric chains, if any---------------------------------------------------------------------------------------------

    res_checked = False #make sure the structure contains REMARK 2; might be unnecessary depending on whether it's ever missing
    author_found = False #prevent the software from reading any software-determined-only biomolecules
    #if any author-determined ones are present; assumes that the author-determined ones come first
    any_chains = False #report if no chains are found
    rem350 = False #report if no remark 350 is found

    look_for_mer = False #read only the first line after each biomolecule line to avoid reading the
    #software-determined biological unit for a given biomolecule if the author-determined one has already been read
    look_for_chains = False

    chain_buffer = [] #chains for a single biomolecule are sometimes spread onto multiple lines

    true_type = [] #all apo chains, and holo chains with a biological ligand
    false_type = [] #only used for holo chains with no biological (valid) ligand

    lig_chains = [] # [holo] chains containing biological ligands

    if type == "holo":
        holo_ligs = getligs(struct, "")
        lig_chains = [l[1] for l in holo_ligs]

    for line in open(f'{blast_directory}/rcsb_pdb_{type}/{struct}.pdb'): #get pdb apo indices

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
                        log_prot_rejection(f"{struct} has inadequate resolution of {float(pdblist[8])} A")
                        return False
                    else: #if the resolution is adequate
                        #print("resolution of %s A" % float(pdblist[8])) #not really useful for debugging
                        res_checked = True
                except ValueError:
                    print(f"nonstandard REMARK 2: {line}")
                    log_prot_rejection(f"{struct} has a nonstandard REMARK 2: {line}")
                    return False


        #determine if the biological unit is a monomer, look for the chains only if it is
        if look_for_mer and (line[0:45] == "REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT:" or (line[0:52] == "REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE:" and (not author_found))):
            if line[0:45] == "REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT:":
                author_found = True

            asm_type = re.split(": | ", line)[6]
            if asm_type != "MONOMERIC":
                print(f"skipping {asm_type} assembly")
                log_prot_rejection(f"{struct} is a {asm_type} assembly")
            else:
                print(asm_type)
                look_for_chains = True
                #biomolecule += 1
        elif look_for_mer and not (line[0:52] == "REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE:" and author_found):
            print("missing biological unit; skipping structure")
            log_prot_rejection(f"{struct} is missing a biological unit assignment")

            return False

        #compile a list of all pdb-chains in the biological unit and then extract them to a new file
        #since non-protein ligands may have their own pdb chains, even monomeric proteins may need to be assembled from multiple chains
        if look_for_chains and line[0:41] == "REMARK 350 APPLY THE FOLLOWING TO CHAINS:":
            #beware that on the order of 1/10000 of pdb structures have double spaces before the first chain ID, which this code should tolerate
            pdblist = re.split(',| ', line[41:])
            for i in range(0,len(pdblist)):
                if pdblist[i] != "" and pdblist[i] != "\n":
                    chain_buffer.append(pdblist[i])
        elif look_for_chains and line[0:19] == "REMARK 350   BIOMT1":
            #chains.append(chain_buffer)
            #extract_chain(struct, type, type, chain_buffer)
            if type == "holo" and set(chain_buffer).isdisjoint(lig_chains): #separate holo monomers with no drug-like ligands and save them as apo structures
                false_type.append('_'.join(chain_buffer))
                extract_chain(struct, type, "apo", chain_buffer)
                print(f"apo chains from holo structures: {chain_buffer}")
            else:                                                           #save apo monomers and holo monomers with drug-like ligands
                true_type.append('_'.join(chain_buffer))
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
        log_prot_rejection(f"{struct} missing remark 2")
        return False
    elif not rem350:
        print("no remark 350 found; skipping structure")
        log_prot_rejection(f"{struct} missing remark 350")
        return False
    elif not any_chains:
        print("no monomers found")
        log_prot_rejection(f"{struct} has no monomeric chains")
        return False

    return [true_type, false_type]


#save protein only; currently unused; could be improved to check for existing structures before saving
def saveprot(xtal, id, chain, btype):

    try:
        xtal.atom_slice(xtal.top.select("protein")).save_pdb(f"paired_pdb_structures/{id}_chain{chain}_{btype}.pdb")
        print(f"PDB id {id} chain {chain} saved")

    except ValueError:
        print("mdtraj failed to extract protein due structure anomaly")

#check that the apo structure is really an apo structure and that its active site is free of crystallographic additives
#[holo crystal structure in mdtraj, holo pdb id, holo chain(s), apo crystal structure in mdtraj, apo pdb id, apo chain(s),]
def checkligs(holo_xtal, idh, chh, apo_xtal, ida, cha, holo_ligand_iis):
    #print("timedebug: ")
    #-----------------------------------------------------------------prepare apo ligands-------------------------------------------------------

    solvent = ["HOH","DOD","NA","CL","ZN"] #molecules which do not count as ligands;
    #zinc is in here since it's usually part of the protein in both apo and holo structures
    #DOD is heavy water, pdb appears to lack molecule types for semiheavy water (HOD/DOH)
    apo_ligand_list = []

    #should report if the a ligand ends up as the 0th chain due to numbering issues
    #(thus far only nucleic acids have done this)
    if len(apo_xtal.top.to_fasta(chain=0)) < 10:
        print("anomalously short protein chain: please inspect")
        log_prot_rejection(f"{idh}_{chh} has an anomalously short protein chain; skipping it")
        return [[],{},0]

    for i_ach in range(1,apo_xtal.top.n_chains): #for each apo chain except the first
        smallch = apo_xtal.atom_slice(apo_xtal.top.select(f"chainid {i_ach}")) #make a trajectory out of that chain

        for r_ach in smallch.top.residues: #for each residue in that chain
            resi_name = str(r_ach)[0:-len(str(r_ach.resSeq))] #separate residue name from number

            if resi_name not in solvent:
                apo_ligand_list.append([str(r_ach)[0:-len(str(r_ach.resSeq))],"",str(r_ach.resSeq)])
                #no chain ID is included since the pdb and mdtraj chain numbering/lettering systems don't match


    #-------------------------------------------------------------align apo and holo structures------------------------------------------------

    #print("timedebug: RMSD 1")

    align_inds = get_indices(ida, cha, idh, chh) #get mdtraj indices of shared residues to use for alignment

    if align_inds == []:#[0] == [] or align_inds[1] == []: #exit if get_indices encountered mismatched numbering
        return [[],{},0]

    #print(align_inds)
    #print("timedebug: RMSD 1.1.0")

    #regex-based ones that failed
    #ca_select_query = "residue =~ '1|2|3'"
    #ca_select_query = "residue =~ '" + "|".join(align_inds) + "' and name CA"

    #ca_select_query = "(residue " + " or residue ".join(align_inds) + ") and name CA"
    #print(ca_select_query)
    #print(apo_xtal.top.select(ca_select_query))
    #print("timedebug: RMSD 1.1.1")
    #t1_start = perf_counter()
    ############################## v most of the program's time is spent here v ####################################

    #print([a for a in apo_xtal.top.atoms])

    #atom ([resn][resi]-[atomname])
    #1. split on "-" -> b
        #atomname is b[1]
        #resi is b[0][3:]

    #atom.split("-")[1]=="CA" and atom.split("-")[0][3:] in align_inds for x



    query_a = np.array([x for x, atom in enumerate(apo_xtal.top.atoms) if str(atom).split("-")[1]=="CA" and str(atom).split("-")[0][3:] in align_inds and re.match("\w\w\w",str(atom).split("-")[0][3:]) != None]) #note that if/else goes before 'for' but if alone goes after
    query_h = np.array([x for x, atom in enumerate(holo_xtal.top.atoms) if str(atom).split("-")[1]=="CA" and str(atom).split("-")[0][3:] in align_inds and re.match("\w\w\w",str(atom).split("-")[0][3:]) != None]) #note that if/else goes before 'for' but if alone goes after

    #print(query_a2)
    #query_a2 = apo_xtal.top.select(ca_select_query)
    #print(query_a)
    #query_h = holo_xtal.top.select(ca_select_query)

    ############################## ^ most of the program's time is spent here ^ ####################################
    #t1_stop = perf_counter()

    #print(t1_stop-t1_start)
    #print("timedebug: RMSD 1.1.2")


    #query_a = np.array([apo_xtal.top.select(f"residue {ind} and name CA")[0] for ind in align_inds])
    #query_h = np.array([holo_xtal.top.select(f"residue {ind} and name CA")[0] for ind in align_inds])



    #get atom indices in apo and holo structures from the residue indices,
    #using alpha carbons only so the number of atoms matches even if there are missing sidechain atoms
    #query_a = np.array([])
    #for ind in align_inds:#[0]:
        #print(ind)
        #if ind == 10:
        #    print(query_a)
        #    print(apo_xtal.top.select(f"residue {ind} and name CA"))
    #    query_a = np.concatenate((query_a, apo_xtal.top.select(f"residue {ind} and name CA")))

    #query_h = np.array([])
    #for ind in align_inds:#[1]:
    #    query_h = np.concatenate((query_h, holo_xtal.top.select(f"residue {ind} and name CA")))

    #print("timedebug: RMSD 1.2")

    #print(query_a)
    #print(query_h)

    #The 'if' case might not be redundant with the check of align_inds above in the ridiculous corner case where
    #the numbering is shifted such that there's a single residue of overlap but that residue is missing an alpha carbon
    if len(query_a) == 0 or len(query_h) == 0:
        print("numbering mismatch in pdb files; may be useable if renumbered")
        log_prot_rejection(f"{idh}_{chh} and {ida}_{cha} have a numbering mismatch in the pdb files")
        return [[],{},0]
    elif len(query_a) != len(query_h):

        #print(align_inds)
        #print([atom for x, atom in enumerate(apo_xtal.top.atoms) if str(atom).split("-")[1]=="CA" and str(atom).split("-")[0][3:] in align_inds and re.match("\w\w\w",str(atom).split("-")[0][3:]) != None])
        #print([atom for x, atom in enumerate(holo_xtal.top.atoms) if str(atom).split("-")[1]=="CA" and str(atom).split("-")[0][3:] in align_inds and re.match("\w\w\w",str(atom).split("-")[0][3:]) != None])
        #print(query_a)
        #print(query_h)

        print("numbering mismatch in mdtraj, probably arising from a missing internal alpha carbon")
        log_prot_rejection(f"{idh}_{chh} and {ida}_{cha} have a numbering mismatch in mdtraj, probably arising from a missing internal alpha carbon")
        return [[],{},0]
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

    #apo_calphas = apo_xtal.atom_slice(query_a.astype(int))
    #apo_calphas.save_pdb(f"./{ida}_chain{cha}_ca_prealign.pdb")

    #print("timedebug: RMSD 2")    #align apo and holo structures using alpha carbons

    #consider aligning using only ligand-coordinating residues for structures with large conformational changes
    apo_xtal.superpose(holo_xtal, atom_indices = query_a.astype(int), ref_atom_indices = query_h.astype(int)) #modifies apo_xtal in place

    apo_calphas = apo_xtal.atom_slice(query_a.astype(int))
    holo_calphas = holo_xtal.atom_slice(query_h.astype(int))

    pair_rmsd = md.rmsd(apo_calphas, holo_calphas) #does not appear to read atom_indices arguments

    #for debuging alignment
    #apo_calphas.save_pdb(f"{directory}/qc_structures/{ida}_chain{cha}_ca.pdb")
    #holo_calphas.save_pdb(f"{directory}/qc_structures/{idh}_chain{chh}_ca.pdb")

    #------------------------------------------------------------------get holo ligand(s)-------------------------------------------------------

    holo_ligand_list = refine_ligand_list(idh, chh)
    if holo_ligand_list == []:
        print("no acceptable holo ligands") #this is a relatively rare problem so reaching it later is okay in terms of efficiency
        log_prot_rejection(f"{idh}_{chh} has no acceptable holo ligands")
        return [[],{},0]

    #---------------------check whether the binding site of each candidate holo ligand is unoccupied in the candidate apo structure------------------

    #beware that "resid" and "residue" in mdtraj are NOT the same!

    all_holo_coords = {}

    okay_ligands = [] #holo ligands with unoccupied binding sites in the apo structure

    #for each holo ligand:
    #   for each apo ligand, check if it overlaps with an apo ligand
    #   if no apo ligands are present at the holo site, add the holo ligand to the list of acceptable ligands

    for i in holo_ligand_list: #run through holo candidate ligands

        if holo_ligand_iis == {}: #calculate holo ligand/pocket indices if none were provided

            #note that single quotes are needed for mdtraj to load ligands which have numerals as the first character of their names

            #select ligand atoms
            if len(i[0].split(" ")) == 1: #select single ligands

                #note that one ligand goes by 986 in rcsb, 098 in moad, and 98 in mdtraj
                #L9K and L9W are actually two different tautomers of the same compound, but rcsb only includes L9W, presumably because the model has no hydrogens so the tautomer is not resolved
                #see http://polymorph.sgc.utoronto.ca/drugged_human_proteome/pages/col_12_IPR001680.html for more information
                #MOAD's TNR is broken into two residues in the modern pdb; the TNR residue is listed as obsolete in the european PDB

                # [MOAD residue name]:[mdtraj residue name]
                mdtraj_moad_ligand_mismatches = {"ADE":"A","URA":"U","BAM":"BEN","DCY":"DC","098":"98","L9K":"L9W","TNR":"A2G' or (resname SER and resSeq 906) or resname 'Null"}
                #note that "or resname 'Null" exists to make the single quotes match up properly when the ligand select string is assembled without having to modify the code below

                if i[0] not in mdtraj_moad_ligand_mismatches.keys():
                    rname = i[0]
                else:
                    rname = mdtraj_moad_ligand_mismatches[i[0]]
                ligand_select_str = f"resname '{rname}'"
            else:                      #select multiple ligands
                #assemble query for all ligands listed
                ligand_resn_list = i[0].split(" ")
                ligand_select_str = f"(resname '{ligand_resn_list[0]}'"
                for lname in ligand_resn_list[1:]:
                    ligand_select_str += f" or resname '{lname}'"
                ligand_select_str += ")"

            #deubgging code for unusual residues--------------------------------
            #print(ligand_select_str)
            #print([i for i in holo_xtal.top.residues])
            #print(holo_xtal.top.select(f"{ligand_select_str} and not element H"))
            #-------------------------------------------------------------------

            #deprecated because the underlying bug is now understood
            #if len(holo_xtal.top.select(f"resname '{i[0]}'")) == 0:
            #    print("mdtraj failed to load the ligand")
            #    log_prot_rejection(f"{idh}_{chh}_{i[0]} could not be loaded by mdtraj")
                #np.save(f"{directory}/warning_files/{idh}_{chh}_{i[0]}_not_loaded", i[0])
                #with open(f"{directory}/note_files/protein_rejection_reasons.log", 'a') as f:
                #    f.writelines(f"{idh}_{chh}_{i[0]}_not_loaded")
            #    return [[],{},0]

            holo_coords = holo_xtal.xyz[0][holo_xtal.top.select(f"{ligand_select_str} and not element H")] #add ligand coordinates; [0] is frame of mdtraj trajectory;

            if len(holo_coords) == 0:
                print("mdtraj failed to load the ligand")
                log_prot_rejection(f"{idh}_{chh}_{ligand_select_str} could not be loaded by mdtraj")

                return [[],{},0]

            if include_resis: #if ligand-coordinating residues are included in distance calculations

                #Distance computations to identify ligand-coordinating holo atoms
                ligands = holo_xtal.top.select(f"{ligand_select_str} and not element H")
                protein = holo_xtal.top.select("protein and not element H")
                #heavy = holo_xtal.top.select_atom_indices("heavy")
                #heavy_prot = np.intersect1d(protein,heavy)

                lig_prot_pair_iis = np.array(list(itertools.product(ligands, protein)))
                prot_lig_dists = md.compute_distances(holo_xtal,lig_prot_pair_iis).flatten()

                lig_coord_iis = np.where(prot_lig_dists<coord_threshold)[0] #within 5 angstroms of ligand
                prot_iis = np.unique([lig_prot_pair_iis[i][-1] for i in lig_coord_iis]) #should save these here to save time later

                #create array of ligand and ligand-binding residue coordinates
                #could run aground on ligands which have numerals as the first character of their names
                lining_resis = []

                for k in prot_iis: #loop through atom indices
                    resi = holo_xtal.top.atom(k).residue.resSeq #get residue containing the atom
                    if resi not in lining_resis: #avoid adding the same residues multiple times
                        #print(holo_xtal.top.select(f"residue {str(resi)}"))
                        holo_coords = np.vstack((holo_coords, holo_xtal.xyz[0][holo_xtal.top.select(f"residue {str(resi)}")])) #add coordinates of each ligand-binding residue
                        lining_resis.append(resi)

                #print("+".join([str(i) for i in lining_resis]))
                #print("+".join([str(i) for i in lining_resis]))


            all_holo_coords[i[0]] = holo_coords #add values to global table thereof to avoid recalculating them

        else:
            holo_coords = holo_ligand_iis[i[0]] #use existing holo ligand/pocket indices for the ligand at hand if available
        #print(i)
        #t2_start = perf_counter()
        if ligdists(apo_ligand_list, apo_xtal, holo_coords, i): #check holo-ligand to apo-ligand distances for current holo ligand
             okay_ligands.append(i)
        #t2_stop = perf_counter()
        #print(f"distance calculations: {t2_stop-t2_start}")

    return [okay_ligands, all_holo_coords, pair_rmsd] #if no issues are found; okay_ligands may be empty


def ligdists(apoligandlist, apo_xtal, holo_coords, i):

    for j in apoligandlist: #loop through all nonstructural ligands in the apo structure
        apo_coords = apo_xtal.xyz[0][apo_xtal.top.select(f"resname '{j[0]}' and not element H")]
        interligdist = np.min(scipy.spatial.distance.cdist(apo_coords, holo_coords))
        if interligdist < dist_threshold:
            print(f"{j[0]} too close to {i[0]}; {interligdist} nm")
            return False

        #print(j)
        #for ac in np.array(apo_coords):
        #    for hc in np.array(holo_coords):
        #        dist = np.linalg.norm(ac-hc) #scipy.spatial.distance.cdist(ac, hc)
                #print(dist)
                #if any apo ligand is in the holo site, exit and do not add the holo ligand to the list
        #        if dist < dist_threshold:
        #            print(f"{j[0]} too close to {i[0]}; distance = {dist}")
        #            return False

    #If there were no apo ligands in the holo site, add the holo ligand to the list
    return True


#get structural and biological ligands for getprot()
def getligs(struct, chain):

    #WARNING: The use of --no-check-certificate compromises security
    #ligand information (contains validity information not found in pdb structure)
    if struct not in existing_xids_moad:
        os.system(f"wget https://www.bindingmoad.org/files/csv/{struct}.csv -O {blast_directory}/moad_xml/{struct}.csv --no-check-certificate")
        existing_xids_moad.append(struct)
        #keeps the list of existing xml files up to date for when the code encounters apo candidates which are in moad and were previously loaded as holo candidates
        #note that this does not track the extant argument and may require manual adjustment?

    bio_ligands = [] #get all the valid ligands from the moad structure

    chain_elems = chain.split("_")

    with open(f'{blast_directory}/moad_xml/{struct}.csv') as csvfile:
        reader = csv.reader(csvfile)
        firstrow = True
        for row in reader:
            if firstrow: #skip first row, which lacks ligand information
                firstrow = False
                continue

            contents = row[3].split(":")
            if (contents[1] in chain_elems or chain == "") and row[4]=="valid":
            #the latter case (chain == "") is for use of this method in get_struct(),
            #where it is used to obtain a list of all chains in a pdb structure containing biologically relevant ligands
                bio_ligands.append(contents)

    return bio_ligands #[keep_ligands, bio_ligands, all_ligands, nst_ligands]


#for holo ligands
def refine_ligand_list(idh, chh):

    ligandlist = getligs(idh, chh) #get holo ligand information

    aa_resns = ["ASN", "ASP", "GLN", "GLU", "THR", "SER", "LYS", "ARG", "HIS", "PRO", "GLY", "CYS", "MET", "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP"]

    #check holo structure for duplicate and amino acid ligands which cryptosite excluded, and peptide ligands which we don't want*
    shift = 0
    for i in range(len(ligandlist)):
        lig_a = ligandlist[i-shift] #avoid skipping elements since indices change

        if len(lig_a[0].split(" ")) > pept_threshold:
            print(f"removing peptide or other polymeric ligand {lig_a[0]} of length {len(lig_a[0].split(' '))} exceeding cutoff {pept_threshold}")
            ligandlist.pop(i-shift) #remove() wasn't working for lists within lists
            shift+=1 #avoid skipping elements and running off the end since indices change; just changing i doesn't work since it skips 0
        elif lig_a[0] in aa_resns:
            print(f"removing {lig_a[0]} amino acid ligand")
            ligandlist.remove(lig_a)
            shift+=1 #avoid skipping elements and running off the end since indices change; just changing i doesn't work since it skips 0
        #else: #remove identical duplicate ligands, which is currently unused since it's not clear why duplicates should be a problem
        #    for j in ligandlist:
        #        if lig_a != j:
        #            if lig_a[0] == j[0]:
        #                ligandlist.remove(lig_a)
        #                ligandlist.remove(j)
        #                print("removing duplicate ligands")


    return ligandlist


#identify residues found in both apo and holo structures and determine their indices in mdtraj for use in alignment
def get_indices(apo_id, apo_chain, holo_id, holo_chain):
    #pdb indices
    apo_indices = []
    holo_indices = []
    apo_resns = {}
    holo_resns = {}
    shared_indices = []

    #mdtraj indices
    apo_md_indices = []
    holo_md_indices = []

    for line in open(f'{directory}/monomer_apo/{apo_id}_chain{apo_chain}.pdb'): #get pdb apo indices
        if line[0:5] == "ATOM " and line[13:16] == "CA ":
            apo_indices.append(int(line[22:26]))
            apo_resns[int(line[22:26])] = line[17:20]

        elif line[0:6] == "HETATM": #handle nonstandard protein resiudes so they don't mess up the alignment
            if line[13:16] == "CA ":
                if line[17:20] == "MSE":
                    apo_resns[int(line[22:26])] = "MET"
                elif line[17:20] == "SEC":
                    apo_resns[int(line[22:26])] = "CYS"
                else:
                    print(f"nonstandard residue {line[17:20]}; skipping structure")
                    log_prot_rejection(f"{apo_id} contains nonstandard residue {line[17:20]}")
                    return [[],[]]
                apo_indices.append(int(line[22:26]))

            elif line[17:20] != "MSE" and line[17:20] != "SEC":
                print(f"nonstandard residue {line[17:20]}; skipping structure")
                log_prot_rejection(f"{apo_id} contains nonstandard residue {line[17:20]}")
                return [[],[]]

        elif line[0:4] == "TER ": #avoid reading post-protein heteroatoms
            break


    for line in open(f'{directory}/monomer_holo/{holo_id}_chain{holo_chain}.pdb'): #get pdb holo indices
        if line[0:5] == "ATOM " and line[13:16] == "CA ":
            holo_indices.append(int(line[22:26]))
            holo_resns[int(line[22:26])] = line[17:20]

        elif line[0:6] == "HETATM": #handle nonstandard protein resiudes so they don't mess up the alignment
            if line[13:16] == "CA ":
                if line[17:20] == "MSE":
                    holo_resns[int(line[22:26])] = "MET"
                elif line[17:20] == "SEC":
                    holo_resns[int(line[22:26])] = "CYS"
                else:
                    print(f"nonstandard residue {line[17:20]}; skipping structure")
                    log_prot_rejection(f"{holo_id} contains nonstandard residue {line[17:20]}")
                    return [[],[]]
                holo_indices.append(int(line[22:26]))

            elif line[17:20] != "MSE" and line[17:20] != "SEC":
                print(f"nonstandard residue {line[17:20]}; skipping structure")
                log_prot_rejection(f"{holo_id} contains nonstandard residue {line[17:20]}")
                return [[],[]]

        elif line[0:4] == "TER ": #avoid reading post-protein heteroatoms
            break

    #eliminate duplicate resi numbers from residues with multiple conformations
    apo_indices = np.unique(apo_indices).tolist()
    holo_indices = np.unique(holo_indices).tolist()

    print(apo_indices == holo_indices)
    print(apo_resns == holo_resns)
    #print(holo_indices)
    #print(holo_resns)

    #convert pdb indices to mdtraj indices
    for i in range(max(apo_indices[0], holo_indices[0]), min(apo_indices[-1], holo_indices[-1])+1):
        if i in apo_indices and i in holo_indices:
            if holo_resns[i] != apo_resns[i]:
                print(f"numbering mismatch in pdb files at residue {i}; may be useable if renumbered")

                #debugging code
                #print(apo_resns)
                #print(holo_resns)
                #for i in range(max(apo_indices[0], holo_indices[0]), min(apo_indices[-1], holo_indices[-1])+1):
                #    if i in apo_indices and i in holo_indices:
                #        print(i, holo_resns[i], apo_resns[i])
                #sys.exit(0)

                log_prot_rejection(f"{apo_id} and {holo_id} differ at residue {i} and may have mismatched numbering")
                return []#[[],[]]
            shared_indices.append(str(i))
            #apo_md_indices.append(apo_indices.index(i))
            #holo_md_indices.append(holo_indices.index(i))

    return shared_indices #[apo_md_indices, holo_md_indices]

#look for gaps in the protein structure exceeding three residues and find the residues at the edges of gaps to check their distance to the holo candidate ligand
def check_gaps(type, pdb_id, chain):

    max_gap_len = 999 #essentially any gap is accepted; this value was formerly 3

    reading = False #skip header
    gap_ends = [] #residue ids of residues on either side of each internal gap
    lastid = 0
    gaplen = 0

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

    for line in open(f'{directory}/monomer_{type}/{pdb_id}_chain{chain}.pdb'): #track gap length
        pdblist = line.split()

        if reading and len(pdblist)==5: #read data
            #detect cases of residue indices with multiple different residue types (not just different conformations) and skip them
            try:
                res_id = int(pdblist[4])
            except ValueError:
                print("residue with [alphabetic] insertion code; skipping structure")
                log_prot_rejection(f"{pdb_id}_{chain} contains an [alphabetic] insertion code")
                return False

            if res_id > first_resi and res_id < last_resi: #don't read missing terminal residues
                if res_id != lastid+1: #start a new gap
                    gap_ends.append(lastid+1) #record the residue located after the end of the previous gap
                    lastid = res_id
                    gaplen = 1
                    gap_ends.append(res_id-1) #record the residue located before the start of the current gap
                else: #track length of existing gap
                    lastid += 1
                    gaplen += 1
                    if gaplen == max_gap_len+1:
                        print("gap too long")
                        log_prot_rejection(f"{pdb_id}_{chain} has a gap exceeding {max_gap_len} residues")
                        return False

        elif pdblist == ['REMARK', '465', 'M', 'RES', 'C', 'SSSEQI']: #find end of header and initiate reading
            reading = True
        elif reading and (pdblist[0] != 'REMARK' or pdblist[1] != '465'): #stop reading and return the gap-lining residues
            if gap_ends != []:
                gap_ends = np.unique((gap_ends[1:]).append(lastid+1))
                #the first element was just the initial lastid value and not a real post-gap residue
                #adds the residue after the final gap since it won't have been added above
                #if gaps are separated by only one residue some residues will be listed twice

                #TODO: add code to inspect gap-end to ligand distance

            return True

    #print("no REMARK 465 present") #this occurs for proteins with no missing residues
    return True


#--------------------------------------------------------------------------------------------
#                           Global Variables and Main Loops
#--------------------------------------------------------------------------------------------

include_resis = False #whether to include holo ligand coordinating residues in determining which apo ligands are too close
#this could be achieved by setting coord_threshold to 0 but that would be computationally inefficient compared with never
#calculating the distances in the first place
dist_threshold = 0.5 #[nanometers] minimum acceptable separation in nm between apo and holo ligands for the apo ligands to be considered outside of the holo binding site
coord_threshold = 0.5 #[nanometers] all residues with atoms within this distance of the holo ligand in the holo structure are considered ligand-coordinating
pept_threshold = 3   #any peptide longer than this is considered a protein chain, any peptide as long or shorter is considered a ligand

existing_xids_moad = [i[0:4] for i in os.listdir(f"{blast_directory}/moad_xml/")] #get a list of already-downloaded ligand lists to avoid downloading extra copies
existing_hids = [i[0:4] for i in os.listdir(f"{blast_directory}/rcsb_pdb_holo/")] #get a list of already-downloaded candidate holo structures to avoid downloading extra copies
existing_aids = [i[0:4] for i in os.listdir(f"{blast_directory}/rcsb_pdb_apo/")] #get a list of already-downloaded candidate apo structures to avoid downloading extra copies
#existing_saved = os.listdir("./paired_pdb_structures") #unimplemented; the save function is unused

blast_out = os.listdir(f"{blast_directory}/blast_output/") #filenames for blastp results of candidate holo (MOAD) structures with at least one candidate apo (PDB blastp hit) structure not itself from MOAD

pair_ligands = [] #all (if any) suitable ligands found for each apo-holo pair
lowrmsd_pairs = [] #the lowest rmsd (if any) suitable ligands found for each apo-holo pair; used to filter out pockets which open spontaneously

#compares all sub_chains of two pdb structures against each other, or compares liganded and ligandless 'holo' chains of one holo structure
#note that this method relies on side effects and global variables declared immediately above
    #for the first call:  holo_chain_list = holo_chains, apo_chain_list = apo_chains[0],                holo_sid=holo_id, apo_sid=apo_id
    #for the second call: holo_chain_list = holo_chains, apo_chain_list = holo_chains_unfiltered[1],    holo_sid=holo_id, apo_sid=holo_id
def compare_chains(holo_chain_list, apo_chain_list, holo_sid, apo_sid, ha_type):
    for a_chain in apo_chain_list: #for all the apo candidate bioassemblies
        print(f"{ha_type}apo_chain {a_chain}")
        if check_gaps("apo", apo_sid, a_chain):
            apo = md.load(f"{directory}/monomer_apo/{apo_sid}_chain{a_chain}.pdb")

            for h_chain in holo_chain_list: #for all the holo candidate bioassemblies with ligands and no gaps
                print(f"holo_chain {h_chain}")
                #gaps checked above
                #print("timdebug: startload")
                holo = md.load(f"{directory}/monomer_holo/{holo_sid}_chain{h_chain}.pdb")
                #print("timedebug: endload")

                #if the holo chain has not previously been compared to an apo chain, identify and save the ligand-coordinating residues, if any
                if h_chain not in holo_lig_iis:
                    #print("timdebug: startcheckligs") # <<---- the problem is in here
                    ligands_iis = checkligs(holo, holo_sid, h_chain, apo, apo_sid, a_chain, {})
                    #print("timdebug: endcheckligs")
                    if ligands_iis != [[],{},0]:
                        holo_lig_iis[h_chain] = ligands_iis[1]
                        #update the coordinate dictionary if this is a new holo chain and holo ligand and residue coordinates were obtained
                        #ligand_iis may be empty due to numbering mismatches in addition to the absence of valid holo ligands,
                        #so holo_lig_iis should be left empty when this is the case to allow for the possibility that future chains
                        #will match and produce nonempty holo_indices
                        #this would be better handled if checkligs had different output values for different pair rejection conditions
                    else:
                        continue #proceed to the next holo chain if no valid ligands can be identified; redundant with if statement below

                else: #if the holo chain has previously been compared to an apo chain, load the previously-identified ligand-coordinating residues
                    ligands_iis = checkligs(holo, holo_sid, h_chain, apo, apo_sid, a_chain, holo_lig_iis[h_chain])
                    #if ligands_iis == [[],{},0]:
                    #    continue #proceed to the next holo chain if no valid ligands can be identified


                #save protein structures and pairing information if valid ligands are found
                if ligands_iis[0] != []: #if the structures are okay but the apo canidate actually has a ligand the other indices of ligand_iis will not be empty

                    #ligands = ligands_iis[0]

                    pair_ligands.append([holo_sid, h_chain, apo_sid, a_chain, ligands_iis[2][0], ligands_iis[0]])
                    buffer_pair_ligands.append([holo_sid, h_chain, apo_sid, a_chain, ligands_iis[2][0], ligands_iis[0]])

                    print(f"pair found; rmsd {ligands_iis[2][0]}")
                else:
                    print("no true apo structure found")
                    log_prot_rejection(f"{apo_sid}_{a_chain} is not a true apo structure of {holo_sid}_{h_chain} due to ligand overlap")

#-------------------------------------------------------------Main Loop------------------------------------------------

debug = 0

if debug == 0: #production runs; not for debugging
    rstart = int(sys.argv[1])
    rend = int(sys.argv[2])
elif debug == 1: #debug at a single holo index
    rindex = 9274
    rstart = rindex
    rend = rindex+1
else: #arbitrary range
    rstart = 157
    rend = 500

#nonstandard residue test case: #barrel with lots of MSE: 12434 #calmodulin: #15820

#1q6m has screwy numbering and appears to have too high of an rmsd
#5x93 has an incomplete remark 465

#the following were removed from the blast_output folder
#bad_holo_ids = [5lwr, 5lws, 4P3Q, 4P3R, 4PSS, 4PST, 4PTH, 4PTJ] #the last 6 are from the same paper
#7591 crashes on holo 4p3q/apo 4pdj with a TypeError, probably because the former is a multi-conformation x-ray structure with over 100 states which nearly killed pymol
#not worth figuring out how to process it. 4p3r has the same problem
#20272 and 20273 (holo 5lwr and 5lws) crash on apo structures 5qb7 and 5oyv respectively with an AttributeError for reasons which are probably too rare to be worth debugging
#the offending file looks fine and loads in pymol

if debug == 0 and rstart == 0:
    makedirectories() #generate output subdirectories in iofiles

#--------------------------------------------------loop through holo structures----------------------------------------

for holo_candidate in range(rstart,rend): #loop through candidate holo structures
    holo_id = blast_out[holo_candidate][0:4].lower() #get id at index
    apo_ids = [id.lower() for id in np.load(f"{blast_directory}/blast_output/{holo_id.upper()}_hits.npy")] #load blast results

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
        #holo_gaps = {}
        #apo_gaps = {}

        for h_chain in holo_chains_unfiltered[0]: #index of 0 is for true holo chains
            print(f"checking gaps in holo chain {h_chain}")
            #This may be the last thing printed for holo chains with no apo or ligandless holo chains;
            #I don't want to waste time checking explicitly for this case since exiting silently is fine
            if check_gaps("holo", holo_id, h_chain):
                holo_chains.append(h_chain)

        if holo_chains != []: #continue only if there are suitable holo chains

            buffer_pair_ligands = [] #used to collect only the lowest RMSD apo-holo pair from each sequence to approximate Sun et. al.'s true pocket screen

            #--------------------------------------------------loop through apo structures-----------------------------------

            for apo_id in apo_ids[1:]: #loop through candidate apo structures [=blast results] (1st element is candidate holo id)

#------------------debugging--------------------------
                #if apo_id.lower() != "4o8z":
                #    continue
#-----------------------------------------------------

                print("------------------------------------------------------")
                print(f"apo: {apo_id}")
                #break pdb structures from rcsb into their constituent physiological assemblies and check resolution and oligomerization state
                apo_chains = getstruct(apo_id, existing_aids, "apo")
                if isinstance(apo_chains, list):

                    compare_chains(holo_chains, apo_chains[0], holo_id, apo_id, "")

            compare_chains(holo_chains, holo_chains_unfiltered[1], holo_id, holo_id, "h")

            #sort and save buffer_pair_ligands for analysis in case the script crashes before terminating
            #add the lowest rmsd pair to the list thereof
            if buffer_pair_ligands != []:
                buffer_pair_ligands = sorted(buffer_pair_ligands, key = itemgetter(4), reverse = False)
                np.save(f"{directory}/pairing_index/ligand_pairs_{holo_candidate}", np.array(buffer_pair_ligands, dtype = object)) #save the ranked list of apo-holo pairs
                lowrmsd_pairs.append(buffer_pair_ligands[0]) #add the lowest rmsd pair


#kind of useless because of how I've been assembling the set out of the output files for each holo structure
#save lists of all apo holo pairs and low rmsd pairs only respectively
np.save(f"{directory}/pairing_index/all_ligand_pairs_{rstart}_{rend}", np.array(pair_ligands, dtype = object))
lowrmsd_pairs = sorted(lowrmsd_pairs, key = itemgetter(4), reverse = True)
np.save(f"{directory}/pairing_index/all_lowrmsd_ligand_pairs_{rstart}_{rend}", np.array(lowrmsd_pairs, dtype = object))

print("")
print(pair_ligands)


#TODO:
#rewrite the list below

#sve ligand-coordinating residues for each hit
#auto generate output folders in bash script


#look at the pair at holo index 112
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
