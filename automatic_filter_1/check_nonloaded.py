import os

#compile nonrepetitive lists of the residues and ligands from the warning files

names = os.listdir("./warning_files")

chains = ["A", "B", "C", "D", "E", "F"]

nonloaded = []
in_chain_0 = []

for i in names:
    parts = i.split("_")
    #for j in parts[1:]:
    if parts[-1] == "loaded.npy" and parts[-3] not in nonloaded:
        nonloaded.append(parts[-3])
    elif parts[-1] == "0.npy" and parts[-4] not in in_chain_0:
        in_chain_0.append(parts[-4])

print("ligands not loaded by mdtraj")
print(nonloaded)
print("found in chain 0")
print(in_chain_0)
