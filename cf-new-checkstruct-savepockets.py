import numpy as np
import sys
from pymol import cmd
from pymol import util
import os
from operator import itemgetter
import re

serial = 5

def get_names(x, y): #x is the index of the holo structure of interest, y is a label for saving structures

    numpy = False
    if numpy:
        pair = sorted(np.load(f"../some_lowrmsd_ligand_pairs_v{serial}.npy", allow_pickle = True), key = itemgetter(4), reverse = True)[int(x)]
    else:
        #ignore the ligands; they're wrong
        pdb_names = [["5o2j","A","5o2k","E","2PO"],["5ota","A","4p0i","B","AQQ"],["4r74","D","4r72","A","F6P (2 conf)"],["2gga","A","2gg4","A","GPJ:S3P"],["5yyb","A","5za4","A","NGE"],["4ocp","A","4wh1","A","GN1:ADP-Mg"],["6h76","A","2cey","A","SLB"],["2vpn","B","3gyy","C","4CS"],["5igy","A","5uxa","A","GDP-Ca, ERY"],["1kmp","A","1kmo","A","2xCIT-2Fe"],["2cgl","A","2cgk","A","ATP:LFR"],["5ni6","A","5nia","A","DJ3 (2 conf)"],["5xdt","A","6rvm","C","ZI7-Ca"],["2p70","A","2fjy","A","2xPRZ"],["5u7i","A","3itm","D","7XS:Mg"],["6i11","A","3p53","A","H0H*"],["2b03","A","3fvj","A","TUD"],["3ozt","A","5lqj","C","OZZ-Mg"]]
        holos = [i[0] for i in pdb_names]
        apos = [i[2] for i in pdb_names]

        pairs_all = sorted(np.load(f"../some_lowrmsd_ligand_pairs_v{serial}.npy", allow_pickle = True), key = itemgetter(4), reverse = True)

        pairs = []
        for i in pairs_all:
            if i[0] in holos and i[2] in apos:
                pairs.append(i)
        pair = pairs[int(x)]

    print(pair)

    #load structures from pymol
    #holo
    if len(pair[1].split("_")) > 1: #load multi-chain structures; has had bugs in the past
        print("multiple holo chains; process them manually")
        sys.exit(0)
    else: #load single-chain structures
        holo = f"{pair[0]}{pair[1]}"

    #apo
    if len(pair[3].split("_")) > 1:
        print("multiple apo chains; process them manually")
        sys.exit(0)
    else:
        apo = f"{pair[2]}{pair[3]}"

    holo_lig_resis = [[i[0],i[2]] for i in pair[5]] #actually [resn, resi]

    return [apo, holo, holo_lig_resis, pair, y]

def load(apo_holo_lig):

    apo = apo_holo_lig[0]
    holo = apo_holo_lig[1]

    cmd.delete("all")

    cmd.fetch(holo, f"{holo}", 0, 1, -1, -2, -1, 'pdb', 0)
    cmd.fetch(apo, f"{apo}", 0, 1, -1, -2, -1, 'pdb', 0)

    #from https://pymol.org/dokuwiki/doku.php?id=api:cmd:fetch
    #"""cmd.fetch(string code, string name, int state, init finish,
    #      int discrete, int multiplex, int zoom, string type,
    #      int async, string path, string file, int quiet)"""

    #fetch with async = 0 does not appear to avoid bugs

    cmd.align(holo, apo)

    return apo_holo_lig

def display(apo_holo_lig):

    apo = apo_holo_lig[0]
    holo = apo_holo_lig[1]
    holo_lig_resis = apo_holo_lig[2]

    all_cs_structures = [["2CGA","B","1AFQ","C","0FG"],["4AKE","B","1ANK","B","ANP"],["1RTC","A","1BR6","A","PT1"],["1ADE","A","1CIB","A","IMP"],["1HOO","B","1CIB","A","HDA"],["1CLL","A","1CTR","A","TFP"],["2WGQ","B","1D6Y","B","HY1"],["1ECJ","D","1ECC","B","PCP"],["1DUB","D","1EY3","F","DAK"],["1NUW","A","1EYJ","B","AMP"],["3PUW","E","1FQC","A","GLO"],["1MY1","C","1FTL","A","DNQ"],["1G4E","B","1G67","B","POP/TZP"],["1HAG","E","1GHY","H","121"],["1EX6","A","1GKY","A","5GP"],["1KZ7","D","1GRN","A","AF3"],["1BSQ","A","1GX8","A","RTL"],["1G24","D","1GZF","C","NIR"],["1E2X","A","1H9G","A","MYR"],["1EXM","A","1HA3","B","MAU"],["1IMF","A","1IMB","B","LIP"],["1RDW","X","1J6Z","A","RHO"],["1B6B","A","1KUV","A","CA5"],["1FA9","A","1L5S","B","URC"],["1ALB","A","1LIC","A","HDS"],["1MY0","B","1N0T","D","AT1"],["1ALV","A","1NX3","A","ISA"],["1OK8","A","1OKE","B","BOG"],["1XCG","B","1OW3","B","GDP"],["1Z92","A","1PY2","A","FRH"],["1JWP","A","1PZO","A","CBT"],["1PZT","A","1PZY","D","UDP"],["3HQD","A","1Q0B","B","NAT"],["1BP5","A","1RYO","A","OXL"],["1RRG","A","1S9D","A","AFB"],["2GPO","A","1S9Q","B","CHD"],["2F6V","A","1T49","A","892"],["1TQO","A","1TR5","A","THP"],["1K3F","B","1U1D","F","181"],["1JBU","H","1WUN","H","P5B"],["1XMG","B","1XVC","A","5BR"],["2AKA","A","1YV3","A","BIT"],["2AIR","H","1ZA1","D","CTP"],["3B7D","E","2AL4","F","CX6"],["3CJ0","A","2BRL","A","POO"],["3CJ0","A","3FQK","B","79Z"],["2BU8","A","2BU2","A","TF1"],["3PEO","G","2BYS","J","LOB"],["1K5H","C","2EGH","B","FOM"],["1SWX","A","2EUM","A","LAT"],["2BRK","A","2GIR","B","NN3"],["1UK2","A","2GZ7","A","D3F"],["2CM2","A","2H4K","A","509"],["1NEP","A","2HKA","C","C3S"],["3L7U","C","2HVD","C","ADP"],["1A8I","A","2IEG","B","FRY"],["3CHE","A","2IUZ","B","D1H"],["1H09","A","2IXU","A","MU2"],["2IYT","A","2IYQ","A","ADP/SKM"],["2GFC","A","2JDS","A","L20"],["2ZB1","A","2NPQ","A","BOG"],["1KS9","A","2OFP","A","PAF"],["2OHG","A","2OHV","A","NHL"],["1FVR","A","2OO8","X","RAJ"],["1ZAH","B","2OT1","D","N3P"],["2AX9","A","2PIQ","A","RB1"],["2Q8F","A","2Q8H","A","TF4"],["2WGB","A","2V57","A","PRL"],["1BNC","B","2V5A","A","LZL"],["1RHB","A","2W5K","B","NDP"],["3GXD","B","2WCG","A","MT5"],["2QFO","B","2WI7","A","2KL"],["1QLW","B","2WKW","B","W22"],["2YQC","A","2YQS","A","UD1"],["3FDL","A","2YXJ","A","N3C"],["3BL9","B","3BL7","A","DD1"],["3F74","C","3BQM","C","BQM"],["2H4E","B","3CFN","B","2AN"],["2QLR","C","3DC1","A","AKG"],["2BF3","A","3DHH","E","BML"],["3MN9","A","3EKS","A","CY9"],["1R1W","A","3F82","A","353"],["1SU4","A","3FGO","B","ACP"],["2BLS","B","3GQZ","A","GF7"],["3H5R","A","3H9J","D","APC"],["3NNU","A","3HL7","A","I46"],["1FXX","A","3HL8","A","BBP"],["1NI6","D","3HOK","B","Q80"],["1PKL","B","3HQP","P","ATP/FDP/OXL"],["3DXN","A","3HZT","A","J60"],["1HKA","A","3IP0","A","HHS"],["1W50","A","3IXJ","C","586"],["3KQA","B","3LTH","A","UD1"],["4HB2","C","4HAT","C","LMB"]]

    all_cs_apos = [i[0] for i in all_cs_structures]
    all_cs_holos = [i[2] for i in all_cs_structures]

    #the two apo-holo cross-checks are in case something got scrambled
    if apo[0:4].upper() in all_cs_apos or holo[0:4].upper() in all_cs_apos or apo[0:4].upper() in all_cs_holos or holo[0:4].upper() in all_cs_holos:
        print("One or both of these structures is also in cryptosite.")

    cmd.hide("spheres", "name Na+Cl")
    cmd.hide("lines")
    cmd.hide("everything", "resn HOH")

    try:
        util.cbac(holo) #intermittently glitchy
        util.cbag(apo)
    except:
        print("util glitch; try running the code for a different structure and then rerunning for this one")

    cmd.center(holo)
    cmd.zoom()

    cmd.hide("sticks")
    cmd.show("lines", "not polymer.protein")

    #select ligand by residue number
    for i in holo_lig_resis:
        cmd.show("sticks", f"{holo} and resi {i[1]}") #intermittently glitchy?

def save(apo_holo_lig):

    apo = apo_holo_lig[0]
    holo = apo_holo_lig[1]
    holo_lig_resis = apo_holo_lig[2]
    pair_label = apo_holo_lig[4]

    pdb_output_dir = "../apo_structures_clean"
    npy_output_dir = "../pair_indices"

    coord_threshold = 5 #maximum distance in Angstroms between ligand and coordinating residues

    #save ligandless apo structure.
    #Note that structual ions and other non-protein structural groups such as glycans may be missing.
    #It might be good to add a second line to save structures with likely-structural ions (Zn, Ca, etc.)
    cmd.save(f"{pdb_output_dir}/{apo}_clean.pdb", f" {apo} and polymer.protein")

    coord_ions = ["Ca","Mg","Zn","Fe","PO3"]#,"PO4","SO4"] #ions which are treated as parts of nearby ligands
    ion_coord_radius = 3
    coord_cmd = "resn "+" or resn ".join(coord_ions)

    all_pocket_resis = []
    print(holo_lig_resis)
    #the code does not support saving different ligands with different opening motions
    for i in holo_lig_resis:
        #find ligand-coordinating residues
        #checking both resi and resn avoids any issues with duplicate residues and with residues from the apo chain with the same number as the ligand
        #cmd.select(f"byres({holo} and polymer.protein within {coord_threshold} of resi {i[1]})")
        #cmd.select(f"byres((sele) and polymer.protein within {coord_threshold} of resn {i[0]})")
        #print(f"byres((resn {i[0]} and resi {i[1]}) or ({coord_cmd}) within {coord_threshold} of (resn {i[0]} and resi {i[1]}))")

        #select the ligand and ions coordinating it (i.e. Mg and ATP)
        cmd.select(f"byres(({holo} and resn {i[0]} and resi {i[1]}) or ({holo} and ({coord_cmd})) within {ion_coord_radius} of ({holo} and resn {i[0]} and resi {i[1]}))")

        #select the protein residues coordinating both the [organic] ligand and associated ions
        cmd.select(f"byres({holo} and polymer.protein within {coord_threshold} of (sele))")

        #extract residue numbers from the selection
        stored.testobj = [] #note that "stored" is a special pymol helper variable
        cmd.iterate("(sele)", "stored.testobj.append(resi)")
        pocket_resis = np.sort([int(i) for i in np.unique(stored.testobj)])
        all_pocket_resis += list(pocket_resis)

        #save the residue numbers
        print("+".join([str(i) for i in pocket_resis]))
        np.save(f"{npy_output_dir}/{apo}_pocketresis_{i[0]}-{i[1]}_{pair_label}.npy", pocket_resis)

    #save the pairing index element for the pair of interest
    np.save(f"{npy_output_dir}/{apo}_index_{pair_label}.npy", apo_holo_lig[3])

    #save all ligand-coordinating residue numbers if there are multiple ligands
    if len(holo_lig_resis) > 1:
        np.save(f"{npy_output_dir}/{apo}_pocketresis_all_{pair_label}.npy", np.unique(all_pocket_resis))

@cmd.extend
def checkst(a):

    if a[-1].isdigit(): #if only the index is given, load structure for inspection
        display(load(get_names(a, ""))) #the label "" is unused
    else:
        in_args = re.split('(\d+)', a)
        #\d is the regex for digits,
        #+ uses the largest available consecutive set of digits as a delimiter,
        #() includes the delimiter in the output rather than removing it
        #note that having two delimiters produces three strings, the first of which is empty

        try: #if the structure has already been loaded, use it
            save(get_names(in_args[1], in_args[2]))
        except: #if the structure isn't loaded, load it
            print("target structure not loaded; reloading to save (or initial save attempt crashed)")
            save(load(get_names(in_args[1], in_args[2])))



#to use this script, enter "run cf-new-checkstruct-savepockets.py" once into the command line
#one can then run "checkst i" where i is a nonnegative real integer to examine the structure at that index
#running "checkst is", where i is a nonnegative real integer and s is just the letter s,
#saves the ligand-coordinating residues and cleaned apo structures.
#It automatically reloads the structure if necessary
#to examine the protein alone run "hide lines; hide spheres; hide sticks"
