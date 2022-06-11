import numpy as np
import os
import sys
import csv
import glob

import Bio
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist

#Jonathan Borowsky
#1/24/22
#This progarm compares the sequence identities of every pair of proteins used in
#the GVP cryptic pocket detection project.

################################################################################
#get FASTA amino acid sequences for all proteins
################################################################################

#get lists of proteins
#-------------------------------------------------------------------------------

#proteins taken directly fron PDB
new_bricks_cryptosite = ['5o2k','6hb0','4p0i','2lao','4r72','1urp','2gg4','3ttl','1gud','5za4','1j8f','1s2o','3ugk','5g1m','4wh1','6h8v','2cey','3gyy','5uxa','5nzm','1tm2','1jej','1kmo','2cgk','3ppn','1ezm','4v38','5nia','3kje','6rvm','2fjy','3p53','1y1a','2w9t','1brq','3fvj','6ypk','2hq8','1kx9','1tvq','2oy4','2zku','3nx1','3qxw','3rwv','4i92','4ic4','4w51','5h9a','6e5d','1ofv','1qys','5bvl','2fd7','4hjk','1igd','4tql','4ake','1bsq','1alb','1ex6','1nep','1ni6','2bls','2qfo','3f74','1exm','1ade','1my0','1rhb','1rtc','2cm2']

test_val_set = ['5o2k','6hb0','4p0i','2lao','4r72','1urp','2gg4','3ttl','1gud','5za4','1j8f','1s2o','3ugk','5g1m','4wh1','6h8v','2cey','3gyy','5uxa','5nzm','1tm2','1jej','1kmo','2cgk','3ppn','1ezm','4v38','5nia','3kje','6rvm','2fjy','3p53','1y1a','2w9t','1brq','3fvj','6ypk','2hq8','1kx9','1tvq','2oy4','2zku','3nx1','3qxw','3rwv','4i92','4ic4','4w51','5h9a','6e5d']
bricks = ['1ofv','1qys','5bvl','2fd7','4hjk','1igd','4tql']
#loop modelling and selenomethionine -> methionine conversion was done for
#bricks, but both selenomethionine and methionine have one letter codes M and
#the gaps are a product of limited structural resolution rather than the
#sequence used

#proteins used in five fold cross validation
filenames_ffcv = glob.glob("/project/bowmanlab/ameller/gvp/msft-share/*.pdb")
prot_names = [fname.split("/")[-1].split(".")[0].lower() for fname in filenames_ffcv]
#this set of proteins overlaps with the proteins above
#Some of these proteins are SARS proteins with structures from swiss model
#homology modelling rather than PDB

#list of all proteins for the main loop:
proteins_all = np.unique(new_bricks_cryptosite+prot_names)

#get FASTA sequences
#-------------------------------------------------------------------------------

fasta_dir = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/fasta_files"

existing_ids = [i.split(".")[0] for i in os.listdir(fasta_dir)]

#get FASTA sequences from pdb files
for fname, pname in zip(filenames_ffcv, prot_names):
    if pname not in new_bricks_cryptosite and pname not in existing_ids:
        with open(fname, 'r') as pdb_file:
            #running SeqIO.parse repeatedly causes it to crash with an empty file error
            SeqIO.write(SeqIO.parse(pdb_file, 'pdb-atom'), f"{fasta_dir}/{pname}.fasta", "fasta")
            existing_ids.append(i)

#download cryptosite/brick/new set fasta files
for i in new_bricks_cryptosite:
    if i not in existing_ids:
        os.system(f"wget https://www.rcsb.org/fasta/entry/{i}/display -O {fasta_dir}/{i}.fasta")
        existing_ids.append(i)

################################################################################
#check sequence identity of each pair of proteins
################################################################################

id_cutoff = 0.4 #maximum sequence identity

#eliminate low-identity hits which commonly occur by chance when aligning a
#short sequence to a much longer one
#set to 1 by default, set to 0 to disable this adjustment
len_scale = 0 #1

zmax = len(proteins_all)*(len(proteins_all)-1)/2
z = 0

print(f"aligning {int(zmax)} protein pairs")

for x, prot1 in enumerate(proteins_all):

    print(f"{x}: {prot1}, {'{:.4f}'.format(z/zmax)} complete")
    fasta_seq1 = SeqIO.read(open(f"{fasta_dir}/{prot1}.fasta"),'fasta').seq

    for y, prot2 in enumerate(proteins_all[x+1:]):
            fasta_seq2 = SeqIO.read(open(f"{fasta_dir}/{prot2}.fasta"),'fasta').seq

            #matrix = matlist.blosum62
            alignment = pairwise2.align.globalds(fasta_seq1, fasta_seq2, matlist.blosum62, -11, -1) #parameters matching blastp
            #setting d uses the matrix, s uses the gap creation and extension penalties
            #note that format_alignment() can be used to display alignments nicely as long as they aren't wider than the terminal

            minlen = min(len(fasta_seq1), len(fasta_seq2))
            #use the shorter of the two sequences as the denominator to determine the percent identity
            #as a fraction of the maximum achieveable percent identity given the sequence lengths
            #this approach is more conservative

            maxlen = max(len(fasta_seq1), len(fasta_seq2))

            identical_resis = 0

            for k in range(len(alignment[0].seqA)):
                if alignment[0].seqA[k] == alignment[0].seqB[k]:
                    identical_resis += 1

            if len(alignment[0].seqA)/len(alignment[0].seqB) != 1:
                print("----------------------------------------------------------------\n---------------------------------------eep------------------------------")

            seqid = identical_resis/minlen #fraction of identical residues in aligned sequences

            #only works if the sequence space is only sparsely populated, which it appears to be for the set of filtered pairs
            #the second condition eliminates low-identity hits which commonly occur by chance when aligning a short sequence to a much longer one
            if seqid > id_cutoff: #and seqid > (1-minlen/maxlen)*len_scale:
                if prot1 in test_val_set:
                    p1_class=" (validation/test)"
                elif prot1 in bricks:
                    p1_class=" (brick)"
                else:
                    p1_class=" (5fcv/train)"

                if prot2 in test_val_set:
                    p2_class=" (validation/test)"
                elif prot2 in bricks:
                    p2_class=" (brick)"
                else:
                    p2_class=" (5fcv/train)"

                print(f"{prot1}{p1_class} and {prot2}{p2_class} have {seqid*100}% identity")

            z+=1 ##counter
