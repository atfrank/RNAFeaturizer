
with open('tests/output_native.txt') as f:
    native_dists = eval(f.readline())
with open('tests/output_decoy.txt') as f:
    decoy_dists = eval(f.readline())
# Now a python dictionary
print("Native cavitiy distances:\n", native_dists)
# get a list of neighboring atom types
print("\nList of neighbor atom types:\n", native_dists[1].keys())
# get a list of distances between native cavity and all ADE.N6 atoms
print("\nDistances between native cavity and all ADE.N6 atoms: ", native_dists[1][':ADE.N6'])
# get a list of distances between decoy cavity No.2 and all GUA.C1' atoms
print("\nDistances between decoy cavity No.2 and all GUA.C1' atoms:", decoy_dists[2][":GUA.C1'"])
