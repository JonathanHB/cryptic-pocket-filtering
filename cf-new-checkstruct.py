import numpy as np
import sys
from pymol import cmd
from pymol import util
import os
from operator import itemgetter

def pairnum(pairs, num):
    if num == 0:
        if len(pairs) == 1:
            print("1 pair")
        else:
            print(f"{len(pairs)} pairs")

@cmd.extend
def checkst(a, *argv):

    start = 48
    end = 323
                #0 = all pairs
    type = 2    #1 = low rmsd pairs
                #2 = all pairs from a particular round

    n = int(a)
    if len(argv) != 0:
        m = argv[0]
    elif type == 2:
        print("missing argument")
        sys.exit(0)

    #note that this code may crash if run on pdb structures where the monomer contains multiple nominal chains

    index = 0

    if type == 0:
        pairs = np.load(f"../pairing_index/all_ligand_pairs_{start}_{end}.npy", allow_pickle = True)
        index = n
        pairnum(pairs, n)
    elif type == 1:
        pairs = np.load(f"../pairing_index/all_lowrmsd_ligand_pairs_{start}_{end}.npy", allow_pickle = True)
        index = n
        pairnum(pairs, n)
    elif type == 2:
        pairs = np.load(f"../pairing_index/ligand_pairs_{n}.npy", allow_pickle = True)
        index = m
        pairnum(pairs, m)
    else:
        print("invalid mode")
        sys.exit(0)

    spairs = sorted(pairs, key = itemgetter(4), reverse = True)
    pair = spairs[index]

    #try: #the regular error message is fine for this
    #    pair = spairs[n]
    #except IndexError:
    #    print(f"there are only {len(pairs)} apo-holo pairs")
    #    sys.exit(0)

    #load structures from pymol
    if len(pair[1].split("_")) > 1:
        holo = f"{pair[0]}{'+'.join(pair[1].split('_'))}"
    else:
        holo = f"{pair[0]}{pair[1]}"

    if len(pair[3].split("_")) > 1:
        apo = f"{pair[2]}{'+'.join(pair[3].split('_'))}"
    else:
        apo = f"{pair[2]}{pair[3]}"

    #if len(pair[3].split("_")) > 1:
    #    pair3 = "+".join(pair[3].split("_"))

    #apo = f"{pair[2]}{pair[3]}"

    holo_lig_resis = [i[2] for i in pair[5]]

    cmd.delete("all")

    #from https://pymol.org/dokuwiki/doku.php?id=api:cmd:fetch
    #"""cmd.fetch(string code, string name, int state, init finish,
    #      int discrete, int multiplex, int zoom, string type,
    #      int async, string path, string file, int quiet)"""

    #fetch with async = 0 to avoid bugs
    cmd.fetch(holo, f"{holo}", 0, 1, -1, -2, -1, 'pdb', 0) #fails
    cmd.fetch(apo, f"{apo}", 0, 1, -1, -2, -1, 'pdb', 0)

    cmd.align(apo, holo)

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

    #print(holo_lig_resis)
    for i in holo_lig_resis:
        cmd.show("sticks", f"{holo} and resi {i}") #intermittently glitchy


#to use this script, enter "run cf-new-checkstruct.py" once into the command line
#one can then run "checkst(i)" or "checkst(i,j)" where i is a nonnegative real integer to examine the structure at that index
#to examine the protein alone run "hide lines; hide spheres; hide sticks"
