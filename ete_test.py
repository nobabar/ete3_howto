from ete3 import Tree, PhyloTree, TreeStyle
import pandas as pd

# importer les données additionelles des prots
prot_names = []
prot_dates = []
prot_origin = []

with open("./data/spike_header_cl.txt", "r") as filein:
    for ligne in filein:
        l_spl = ligne.split("|")
        prot_names.append(l_spl[0].split()[0])
        prot_origin.append(l_spl[1].split()[1][1:-1])
        prot_dates.append(l_spl[2].strip())

# organiser as dataframe
spike_df = pd.DataFrame(data = {"name" : prot_names,
                                "origin" : prot_origin,
                                "date" : prot_dates})
spike_df["date"] = pd.to_datetime(spike_df["date"])

# importer l'arbre
t = PhyloTree(newick = "./data/spike_data_708.nwk", alignment = "./data/clustalo_alg.fasta", alg_format = "fasta")

# fooling around
print(f"Est-ce qu'on est à la racine : {t.is_root()}")
print(f"La distance à son premier enfant : {t.children[0].dist}")

leaf_ex, dist_to_leaf = t.get_closest_leaf()
print(f"La distance à la plus courte feuille {leaf_ex} est : {dist_to_leaf}")
ex_up = leaf_ex.up.up.up


# Dessiner les arbres : 
# circulaire?
circ_style = TreeStyle()
circ_style.mode = "c"
circ_style.scale = 20
# marche pas peut-être essayer avec plus petit
ex_up.render("test_render.svg", tree_style = circ_style)

# Imprimer tous les noms de noeuds terminaux
# for node in t.traverse("preorder"):
#     if node.is_leaf():
#         print(node.name)

