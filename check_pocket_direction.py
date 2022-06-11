import numpy as np

#report which of the new set pockets open forward

serial_in = 6 #which version of the pocket set to load
input_directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/new_pockets/iofiles"

forward_pockets = np.load(f"{input_directory}/output_indices/forward_lowrmsd_ligand_pairs_v{serial_in}.npy", allow_pickle=True)

test_val_set = ['5o2k','6hb0','4p0i','2lao','4r72','1urp','2gg4','3ttl','1gud','5za4','1j8f','1s2o','3ugk','5g1m','4wh1','6h8v','2cey','3gyy','5uxa','5nzm','1tm2','1jej','1kmo','2cgk','3ppn','1ezm','4v38','5nia','3kje','6rvm','2fjy','3p53','1y1a','2w9t','1brq','3fvj','6ypk','2hq8','1kx9','1tvq','2oy4','2zku','3nx1','3qxw','3rwv','4i92','4ic4','4w51','5h9a','6e5d']

print("forward pockets")
for i in test_val_set:
    if i in forward_pockets[:,2]:
        print(i)
