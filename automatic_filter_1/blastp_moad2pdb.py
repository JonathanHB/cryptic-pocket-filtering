import re
import numpy as np
import os
import Bio.Blast
from Bio.Blast import NCBIXML

#take the lists of PDB structures downloaded from MOAD (https://www.bindingmoad.org/Home/download)
#for each structure:
    #download the fasta sequence
    #blastp it against PDB structures of at least 50 nucleotides
    #save the resulting structures with single hsps 100% sequence identity and sufficient coverage
    #   note that having a coverage requirement below 1 leaves open the possibility of mutations outside the hsp
    #   but having a coverage requirement of 1 excludes many identical structures of slightly different sequence lengths
    #   (i.e. a disordered terminal tail might be truncated in different places in different structures)

#path to directory
directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/new_pockets/iofiles"

coverage_limit = 0.90 #minimum fractional overlap of allowable hits; filter_blast2pairs.py subsequently makes sure no internal mismatches get through this way
part = "all" #part a or b or both ("all") of the MOAD database, which is broken in two on account of its large size

#get the PDB IDs of all the MOAD structures
structures = np.hstack([np.unique([i[0:4].upper() for i in os.listdir(f"{directory}/every_part_a/BindingMOAD_2020/")]), np.unique([i[0:4].upper() for i in os.listdir(f"{directory}/every_part_b/BindingMOAD_2020/")])])
print(len(structures))
#structures for which fasta files have already been downloaded
existing_ids = [i[0:4] for i in os.listdir(f"{directory}/rcsb_fasta_{part}/")]

#loop through moad structures
for i in structures:

    hits = []
    print(i)

    #download fasta files if necessary
    if i not in existing_ids:
        os.system(f"wget https://www.rcsb.org/fasta/entry/{i}/display -O {directory}/rcsb_fasta_{part}/{i}.fasta")

    #blast the MOAD sequence against PDB
    #note that the ftp version of blastp must be used because the apt-get version is out of date and can't load the pdbaa database.
    #-outfmt 5 is required to produce a biopython-readable .xml file
    os.system(f"/project/bowmanlab/borowsky.jonathan/installations/ncbi-blast-2.12.0+/bin/blastp -db /project/bowmanlab/borowsky.jonathan/installations/pdbaa -query {directory}/rcsb_fasta_{part}/{i}.fasta -outfmt 5 -out {directory}/blast_xml_out/{i}-results.xml")

    #load the blastp results
    result_handle = open(f"{directory}/blast_xml_out/{i}-results.xml")

    #skip structures containing multiple different protein sequences
    #As multiple proteins are usually added to crystallization solutions for the purpose of studying complexes thereof,
    #there are likely relatively few multi protein PDB structures with monomeric bioassemblies
    #so I'm skipping them for the time being for ease of processing
    #The following can be used to read xml results for multimers, but further downstream adjustments would be required
    #  "blast_records = NCBIXML.parse(result_handle)
    #   blast_record = next(blast_records)"
    try:
        blast_record = NCBIXML.read(result_handle)
    except ValueError as err:
        print(err)
        continue

    #for each structure aligned to the moad structure
    for alignment in blast_record.alignments:
        print(alignment.title)
        hsp = alignment.hsps[0] #there should be only one hsp for any alignment good enough to be useful; this is checked below
        #print(hsp.identities)
        #check that there is one hsp with 100% identity and good coverage
        if len(alignment.hsps) == 1 and hsp.identities == hsp.align_length and hsp.align_length/alignment.length > coverage_limit:

            entries = re.split("pdb\||>pdb\|", alignment.title) #separate out all pdb files with a given sequence
            pdb_ids = [title[0:4] for title in entries[1:]] #the first entry is "" from the left of the first pdb| header
            hits = hits+pdb_ids

    #make the query id, if any, the first element in the array (I think it is before np.unique() rearranges things,
    #but I'm not sure and this way if the blast results are in an unusual order it still works)

    hits_u = np.unique(hits)

    hits_apo_a  = [] #all non-moad hits
    hits_apo_m = np.delete(hits_u, np.where(i.upper()==hits_u)) #all hits except the holo structure
    for acand in hits_u: #for all unique hits
        if acand not in structures:
            hits_apo_a.append(acand) #add all non-moad structures

    #add the holo structure id to the start of the list
    hits_a = np.hstack([i.upper(), hits_apo_a])
    hits_m = np.hstack([i.upper(), hits_apo_m])

    print(hits_a)

    #save the output with and without structures with no apo candidates
    if len(hits_a) > 1:
        np.save(f"{directory}/blast_output/{i}_hits.npy", hits_a)
    else:
        np.save(f"{directory}/blast_output_ho/{i}_hits_ho.npy", hits_a)

    if len(hits_m) > 1:
        np.save(f"{directory}/blast_output_m/{i}_hits_m.npy", hits_m)
    else:
        np.save(f"{directory}/blast_output_m_ho/{i}_hits_m_ho.npy", hits_m)

#---------------------------------------------------End of Code---------------------------------------------------

#this script ran successfully from 8/15/21 to 8/16/21, processing all of MOAD
