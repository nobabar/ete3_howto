from ete3 import Tree, PhyloTree, TreeStyle

t = PhyloTree(newick = "./data/spike_data_708.nwk", alignment = "./data/clustalo_alg.fasta", alg_format = "fasta")
print(f"Est-ce qu'on est à la racine : {t.is_root()}")
print(f"La distance à son premier enfant : {t.children[0].dist}")

# Dessiner les arbres : 
# circulaire?
circ_style = TreeStyle()
circ_style.mode = "c"
circ_style.scale = 20
t.render("test_render.svg", tree_style = circ_style)

# # Imprimer tous les noms de noeuds terminaux
# for node in t.traverse("preorder"):
#     if node.is_leaf():
#         print(node.name)

