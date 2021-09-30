import numpy as np
import os
import csv

#used to compile a list of pockets from the structures and data saved by the pymol manual inspection macro
 
#check for structures also found in cryptosite
all_cs_structures = [["2CGA","B","1AFQ","C","0FG"],["4AKE","B","1ANK","B","ANP"],["1RTC","A","1BR6","A","PT1"],["1ADE","A","1CIB","A","IMP"],["1HOO","B","1CIB","A","HDA"],["1CLL","A","1CTR","A","TFP"],["2WGQ","B","1D6Y","B","HY1"],["1ECJ","D","1ECC","B","PCP"],["1DUB","D","1EY3","F","DAK"],["1NUW","A","1EYJ","B","AMP"],["3PUW","E","1FQC","A","GLO"],["1MY1","C","1FTL","A","DNQ"],["1G4E","B","1G67","B","POP/TZP"],["1HAG","E","1GHY","H","121"],["1EX6","A","1GKY","A","5GP"],["1KZ7","D","1GRN","A","AF3"],["1BSQ","A","1GX8","A","RTL"],["1G24","D","1GZF","C","NIR"],["1E2X","A","1H9G","A","MYR"],["1EXM","A","1HA3","B","MAU"],["1IMF","A","1IMB","B","LIP"],["1RDW","X","1J6Z","A","RHO"],["1B6B","A","1KUV","A","CA5"],["1FA9","A","1L5S","B","URC"],["1ALB","A","1LIC","A","HDS"],["1MY0","B","1N0T","D","AT1"],["1ALV","A","1NX3","A","ISA"],["1OK8","A","1OKE","B","BOG"],["1XCG","B","1OW3","B","GDP"],["1Z92","A","1PY2","A","FRH"],["1JWP","A","1PZO","A","CBT"],["1PZT","A","1PZY","D","UDP"],["3HQD","A","1Q0B","B","NAT"],["1BP5","A","1RYO","A","OXL"],["1RRG","A","1S9D","A","AFB"],["2GPO","A","1S9Q","B","CHD"],["2F6V","A","1T49","A","892"],["1TQO","A","1TR5","A","THP"],["1K3F","B","1U1D","F","181"],["1JBU","H","1WUN","H","P5B"],["1XMG","B","1XVC","A","5BR"],["2AKA","A","1YV3","A","BIT"],["2AIR","H","1ZA1","D","CTP"],["3B7D","E","2AL4","F","CX6"],["3CJ0","A","2BRLA/3FQKB","0","POO/79Z"],["2BU8","A","2BU2","A","TF1"],["3PEO","G","2BYS","J","LOB"],["1K5H","C","2EGH","B","FOM"],["1SWX","A","2EUM","A","LAT"],["2BRK","A","2GIR","B","NN3"],["1UK2","A","2GZ7","A","D3F"],["2CM2","A","2H4K","A","509"],["1NEP","A","2HKA","C","C3S"],["3L7U","C","2HVD","C","ADP"],["1A8I","A","2IEG","B","FRY"],["3CHE","A","2IUZ","B","D1H"],["1H09","A","2IXU","A","MU2"],["2IYT","A","2IYQ","A","ADP/SKM"],["2GFC","A","2JDS","A","L20"],["2ZB1","A","2NPQ","A","BOG"],["1KS9","A","2OFP","A","PAF"],["2OHG","A","2OHV","A","NHL"],["1FVR","A","2OO8","X","RAJ"],["1ZAH","B","2OT1","D","N3P"],["2AX9","A","2PIQ","A","RB1"],["2Q8F","A","2Q8H","A","TF4"],["2WGB","A","2V57","A","PRL"],["1BNC","B","2V5A","A","LZL"],["1RHB","A","2W5K","B","NDP"],["3GXD","B","2WCG","A","MT5"],["2QFO","B","2WI7","A","2KL"],["1QLW","B","2WKW","B","W22"],["2YQC","A","2YQS","A","UD1"],["3FDL","A","2YXJ","A","N3C"],["3BL9","B","3BL7","A","DD1"],["3F74","C","3BQM","C","BQM"],["2H4E","B","3CFN","B","2AN"],["2QLR","C","3DC1","A","AKG"],["2BF3","A","3DHH","E","BML"],["3MN9","A","3EKS","A","CY9"],["1R1W","A","3F82","A","353"],["1SU4","A","3FGO","B","ACP"],["2BLS","B","3GQZ","A","GF7"],["3H5R","A","3H9J","D","APC"],["3NNU","A","3HL7","A","I46"],["1FXX","A","3HL8","A","BBP"],["1NI6","D","3HOK","B","Q80"],["1PKL","B","3HQP","P","ATP/FDP/OXL"],["3DXN","A","3HZT","A","J60"],["1HKA","A","3IP0","A","HHS"],["1W50","A","3IXJ","C","586"],["3KQA","B","3LTH","A","UD1"],["4HB2","C","4HAT","C","LMB"]]

all_cs_apos = [i[0] for i in all_cs_structures]
all_cs_holos = [i[2] for i in all_cs_structures]

#run this script from
#any directory contained in the same directory as the output directories from cf-new-checkstruct-savepockets.py
serial = 3

#pdb_output_dir = "../apo_structures_clean" #unused
npy_output_dir = "../pair_indices"
csv_output_name = f"../screened-pairs-{serial}" #also used to save a numpy file

files = os.listdir(npy_output_dir)

pairs = [] #[[apo, apo id, holo, holo id, rmsd, type, in cryptosite, comments],...]


# open the file in write mode
f = open(csv_output_name, 'w', newline='')
# create the csv writer
writer = csv.writer(f, quoting = csv.QUOTE_NONE)

for i in files:
    if i[6:11] == "index": #avoid extraneous files
        file = np.load(f"{npy_output_dir}/{i}", allow_pickle = True)
        print(file)
        if file[2].upper() in all_cs_apos or file[2].upper() in all_cs_apos or file[0].upper() in all_cs_holos or file[0].upper() in all_cs_holos:
            in_cs = "yes"
        else:
            in_cs = "no"
        label = i[11:].split(".")[0] #extract pocket type and comment
        #extract comment, if any
        #positions 0 and 2 are underscores but there are more in the comments themselves in lieu of spaces
        if len(label) > 3:
            comment = label[3:]
        else:
            comment = ""

        print(file[5])
        writer.writerow(list(file[2:4])+list(file[0:2])+list(file[5][0])+[file[4]]+[label[1]]+[in_cs]+[comment])
        # write a row to the csv file

        pairs.append(list(file[2:4])+list(file[0:2])+list(file[5])+[file[4]]+[label[1]]+[in_cs]+[comment])

# close the file
f.close()

np.save(csv_output_name, pairs)
