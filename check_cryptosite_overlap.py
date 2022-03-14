import numpy as np
import os

#the complete structures from the cryptosite paper; [apo pdb id, apo chain, holo pdb id, holo chain, ligand name]
cs_structures=[["2CGA","B","1AFQ","C","0FG"],["4AKE","B","1ANK","B","ANP"],["1RTC","A","1BR6","A","PT1"],["1ADE","A","1CIB","A","IMP"],["1HOO","B","1CIB","A","HDA"],["2WGQ","B","1D6Y","B","HY1"],["1ECJ","D","1ECC","B","PCP"],["1DUB","D","1EY3","F","DAK"],["1EX6","A","1GKY","A","5GP"],["1KZ7","D","1GRN","A","AF3"],["1BSQ","A","1GX8","A","RTL"],["1EXM","A","1HA3","B","MAU"],["1ALB","A","1LIC","A","HDS"],["1MY0","B","1N0T","D","AT1"],["1ALV","A","1NX3","A","ISA"],["1XCG","B","1OW3","B","GDP"],["3HQD","A","1Q0B","B","NAT"],["1BP5","A","1RYO","A","OXL"],["2GPO","A","1S9Q","B","CHD"],["2F6V","A","1T49","A","892"],["2AKA","A","1YV3","A","BIT"],["2AIR","H","1ZA1","D","CTP"],["1K5H","C","2EGH","B","FOM"],["1UK2","A","2GZ7","A","D3F"],["1NEP","A","2HKA","C","C3S"],["3L7U","C","2HVD","C","ADP"],["3CHE","A","2IUZ","B","D1H"],["1H09","A","2IXU","A","MU2"],["1ZAH","B","2OT1","D","N3P"],["3GXD","B","2WCG","A","MT5"],["2QFO","B","2WI7","A","2KL"],["3F74","C","3BQM","C","BQM"],["2BF3","A","3DHH","E","BML"],["1SU4","A","3FGO","B","ACP"],["2BLS","B","3GQZ","A","GF7"],["1NI6","D","3HOK","B","Q80"],["1HKA","A","3IP0","A","HHS"]]

#all blast results
all_out = os.listdir("./blast_output/")

matches = np.zeros(len(cs_structures))

for i in range(len(all_out)):
    for j in range(len(cs_structures)):
        #print(all_out[i][0:4])
        #print(cs_structures[j][2])
        if all_out[i][0:4].lower()==cs_structures[j][2].lower():
            matches[j]=1

print(matches) #note that half the matches are expected in the not-yet-processed part b of MOAD
print(np.count_nonzero(matches == 1))
    #print(np.load(f"blast_output/{all_out[i]}"))

#print(np.load("blast_output/11bg_hits.npy"))
