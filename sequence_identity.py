import numpy as np
import os
import Bio
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
import csv

all_cs_structures = [["2CGA","B","1AFQ","C","0FG"],["4AKE","B","1ANK","B","ANP"],["1RTC","A","1BR6","A","PT1"],["1ADE","A","1CIB","A","IMP"],["1HOO","B","1CIB","A","HDA"],["1CLL","A","1CTR","A","TFP"],["2WGQ","B","1D6Y","B","HY1"],["1ECJ","D","1ECC","B","PCP"],["1DUB","D","1EY3","F","DAK"],["1NUW","A","1EYJ","B","AMP"],["3PUW","E","1FQC","A","GLO"],["1MY1","C","1FTL","A","DNQ"],["1G4E","B","1G67","B","POP/TZP"],["1HAG","E","1GHY","H","121"],["1EX6","A","1GKY","A","5GP"],["1KZ7","D","1GRN","A","AF3"],["1BSQ","A","1GX8","A","RTL"],["1G24","D","1GZF","C","NIR"],["1E2X","A","1H9G","A","MYR"],["1EXM","A","1HA3","B","MAU"],["1IMF","A","1IMB","B","LIP"],["1RDW","X","1J6Z","A","RHO"],["1B6B","A","1KUV","A","CA5"],["1FA9","A","1L5S","B","URC"],["1ALB","A","1LIC","A","HDS"],["1MY0","B","1N0T","D","AT1"],["1ALV","A","1NX3","A","ISA"],["1OK8","A","1OKE","B","BOG"],["1XCG","B","1OW3","B","GDP"],["1Z92","A","1PY2","A","FRH"],["1JWP","A","1PZO","A","CBT"],["1PZT","A","1PZY","D","UDP"],["3HQD","A","1Q0B","B","NAT"],["1BP5","A","1RYO","A","OXL"],["1RRG","A","1S9D","A","AFB"],["2GPO","A","1S9Q","B","CHD"],["2F6V","A","1T49","A","892"],["1TQO","A","1TR5","A","THP"],["1K3F","B","1U1D","F","181"],["1JBU","H","1WUN","H","P5B"],["1XMG","B","1XVC","A","5BR"],["2AKA","A","1YV3","A","BIT"],["2AIR","H","1ZA1","D","CTP"],["3B7D","E","2AL4","F","CX6"],["3CJ0","A","2BRL","A","POO"],["3CJ0","A","3FQK","B","79Z"],["2BU8","A","2BU2","A","TF1"],["3PEO","G","2BYS","J","LOB"],["1K5H","C","2EGH","B","FOM"],["1SWX","A","2EUM","A","LAT"],["2BRK","A","2GIR","B","NN3"],["1UK2","A","2GZ7","A","D3F"],["2CM2","A","2H4K","A","509"],["1NEP","A","2HKA","C","C3S"],["3L7U","C","2HVD","C","ADP"],["1A8I","A","2IEG","B","FRY"],["3CHE","A","2IUZ","B","D1H"],["1H09","A","2IXU","A","MU2"],["2IYT","A","2IYQ","A","ADP/SKM"],["2GFC","A","2JDS","A","L20"],["2ZB1","A","2NPQ","A","BOG"],["1KS9","A","2OFP","A","PAF"],["2OHG","A","2OHV","A","NHL"],["1FVR","A","2OO8","X","RAJ"],["1ZAH","B","2OT1","D","N3P"],["2AX9","A","2PIQ","A","RB1"],["2Q8F","A","2Q8H","A","TF4"],["2WGB","A","2V57","A","PRL"],["1BNC","B","2V5A","A","LZL"],["1RHB","A","2W5K","B","NDP"],["3GXD","B","2WCG","A","MT5"],["2QFO","B","2WI7","A","2KL"],["1QLW","B","2WKW","B","W22"],["2YQC","A","2YQS","A","UD1"],["3FDL","A","2YXJ","A","N3C"],["3BL9","B","3BL7","A","DD1"],["3F74","C","3BQM","C","BQM"],["2H4E","B","3CFN","B","2AN"],["2QLR","C","3DC1","A","AKG"],["2BF3","A","3DHH","E","BML"],["3MN9","A","3EKS","A","CY9"],["1R1W","A","3F82","A","353"],["1SU4","A","3FGO","B","ACP"],["2BLS","B","3GQZ","A","GF7"],["3H5R","A","3H9J","D","APC"],["3NNU","A","3HL7","A","I46"],["1FXX","A","3HL8","A","BBP"],["1NI6","D","3HOK","B","Q80"],["1PKL","B","3HQP","P","ATP/FDP/OXL"],["3DXN","A","3HZT","A","J60"],["1HKA","A","3IP0","A","HHS"],["1W50","A","3IXJ","C","586"],["3KQA","B","3LTH","A","UD1"],["4HB2","C","4HAT","C","LMB"]]


#all_cs_apos = [i[0] for i in all_cs_structures]
all_cs_holos_withdimers = [i[2] for i in all_cs_structures]
all_cs_holos = []
existing_ids = [i[0:4] for i in os.listdir(f"rcsb_fasta_cs/")]

#download cryptosite fasta files and prepare a list of the monomeric ones
for i in all_cs_holos_withdimers:
    if i not in existing_ids:
        os.system(f"wget https://www.rcsb.org/fasta/entry/{i}/display -O rcsb_fasta_cs/{i}.fasta")
        existing_ids.append(i.upper())
    try:
        fasta_seq2 = SeqIO.read(open(f"rcsb_fasta_cs/{i}.fasta"),'fasta').seq
        all_cs_holos.append(i)
    except ValueError:
        print(f"{i} has multiple chains and will be skipped") #if there are two different chains it's not in the new set to begin with

serial_in = 2
serial_out = 3
pairing_index = np.load(f"output_indices/screened-pairs-{serial_in}.npy")
#pairing_index = np.load(f"output_indices/screened-pairs-{serial}.npy")

output_pairs = [] #pairs which are not too close to other pairs
remove_inds = [] #indices of structures to remove
closest_cs = [] #cryptosite structure with the most similar sequence to each ns structure; [identity; length ratio; cryptosite info]


all_ns_holos = [i[2] for i in pairing_index]

id_cutoff = 0.4 #maximum sequence identity
len_scale = 1 #eliminate low-identity hits which commonly occur by chance when aligning a short sequence to a much longer one

#too_similar = []

for i in range(len(all_ns_holos)):
    print(i)
    id_1 = all_ns_holos[i]
    #print(f"holo1: {id_1}")
    fasta_seq1 = SeqIO.read(open(f"rcsb_fasta_all/{id_1.upper()}.fasta"),'fasta').seq
    #print(fasta_seq1.seq)
    for j in range(i+1, len(all_ns_holos)):
        id_2 = all_ns_holos[j]
        #print(f"holo2: {id_2}")
        fasta_seq2 = SeqIO.read(open(f"rcsb_fasta_all/{id_2.upper()}.fasta"),'fasta').seq

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

        seqid = identical_resis/minlen #fraction of identical residues in aligned sequences

        #only works if the sequence space is only sparsely populated, which it appears to be for the set of filtered pairs
        #the second condition eliminates low-identity hits which commonly occur by chance when aligning a short sequence to a much longer one
        if seqid > id_cutoff and seqid > (1-minlen/maxlen)*len_scale:
            #too_similar.append([pairing_index[i], pairing_index[j], seqid])
            print(pairing_index[i])
            print(pairing_index[j])
            print(f"holo {id_1} and {id_2} have {seqid*100}% identity")

            if seqid == 1: #if the structures are identical, remove the higher-rmsd pair
                if pairing_index[i][4] > pairing_index[j][4]:
                    remove_inds.append(i)
                else:
                    remove_inds.append(j)
            else: #remove the lower-rmsd pair for nonidentical structures
                if pairing_index[i][4] < pairing_index[j][4]:
                    remove_inds.append(i)
                else:
                    remove_inds.append(j)

    max_cs_id = [0,[]]

    for j in range(len(all_cs_holos)):
        id_2 = all_cs_holos[j]

        fasta_seq2 = SeqIO.read(open(f"rcsb_fasta_cs/{id_2.upper()}.fasta"),'fasta').seq

        #matrix = matlist.blosum62
        alignment = pairwise2.align.globalds(fasta_seq1, fasta_seq2, matlist.blosum62, -11, -1) #parameters matching blastp
        #setting d uses the matrix, s uses the gap creation and extension penalties
        #note that format_alignment() can be used to display alignments nicely as long as they aren't wider than the terminal

        minlen = min(len(fasta_seq1), len(fasta_seq2))
        #use the shorter of the two sequences as the denominator to determine the percent identity
        #as a fraction of the maximum achieveable percent identity given the sequence lengths
        #this approach is more conservative

        identical_resis = 0

        for k in range(len(alignment[0].seqA)):
            if alignment[0].seqA[k] == alignment[0].seqB[k]:
                identical_resis += 1

        seqid = identical_resis/minlen

        if seqid > max_cs_id[0]:
            max_cs_id = [seqid, id_2]

    closest_cs.append(max_cs_id)


#produce a human-readable csv file
f = open(f"output_indices/screened-seq-pairs-{serial_out}.csv", 'w', newline='')
writer = csv.writer(f, quoting = csv.QUOTE_NONE)

#create a new array without the too-similar structures
for i in range(len(pairing_index)):
    if i not in remove_inds:

        outline = list(pairing_index[i])+list(closest_cs[i])
        output_pairs.append(outline)
        writer.writerow(outline)

print(len(pairing_index)-len(output_pairs))

f.close()

np.save(f"output_indices/screened-seq-pairs-{serial_out}.npy", output_pairs)
#print(too_similar)

#outline
"""
1. load pairing index
1.5. load fasta sequences of cryptosite apo structures
2. for each apo structure:
    3. for each other apo structure [do a triangular matrix of them to avoid duplicates since seqid is commutative]: #ns-ns comparison
        calculate sequence identity
        if seqid > cutoff:
            remove the lower-rmsd pair (or the higher one if seqid = 1)
            append the structures to a list of too-close pairs along with their seqid

    4. for each cryptosite apo structure: #cs-ns comparison
        calculate sequence identity
        if seqid > cutoff:
            #remove the lower-rmsd pair?
            append the structures to a list of too-close pairs along with their seqid


"""
