import pycircos
import matplotlib.pyplot as plt
import re
import pandas as pd
    """_summary_
    converting a metagenomics dataframe to the circos
    visualization for the enhanced visualization. 
    """
class treePLot:
    """_summary_
    reading and importing the files for the association
    of the tree and the phylogenomic abundances of the 
    species
    """
    def __init__(self, tree_type, num, tree):
        self.tree = input("path_to_tree")
        self.num = num
        self.tree_type = tree_type
        self.tarc = Tarc(tree=self.tree, format="self.tree_type", interspace=1)
        self.tcircle = Tcircle(figzie=(self.num, self.num))
        self.classify = []
        self.color = []
    
    def setPhyla(self,dataframe, columnA, 
                    color,columnB, 
                    edgeA, edgeB, 
                    edgecolor,
                    size, width):
        self.edgecolor = edgecolor
        self.size = size
        self.width = width
        """_summary_

        Args:
            dataframe (_type_): _provide_the_dataframe_
            columnA (_type_): _columnA for selection of the bacterial species_
            color (_type_): _color file as plain text_
            columnB (_type_): _columnB for the selection of the bacterial clade_
            edgeA (_type_): _bacterial ID_
            edgeB (_type_): _bacterial phyla as in the species tree_
        """
        self.dataframe = pd.read_csv("dataframe.csv", sep="\t", header=True)
        self.classify = list(filter(None,self.dataframe[self.columnA].dropna().tolist()))
        with open(self.color, "rb") as f:
            text = f.read()
            self.color.append(re.split(r'\n'),text)
            f.close()
        cross_dict_map = dict(zip(self.classify,self.color))
        visual_dict = {}
        clade_dict = {}
        self.dataframe[[edgeA, edgeB]].to_csv("edgesPopulation.csv")
        with open("edgesPopulation.csv") as edges:
            edges.readline()
            for line in edges:
                line = line.strip().split(",")
                if line[1] in cross_dict_map:
                    clade_dict[line[0]] = {"color":cross_dict_map[line[1]], 
                                  "size":self.size, "linewidth":self.width, 
                                          "edgecolor":str(self.edgecolor)}
                else: 
                    clade_dict[line[0]] = {"color":cross_dict_map["Other"], 
                                  "size":self.size, "linewidth":self.width, 
                                          "edgecolor":str(self.edgecolor)}
self.tcircle.add_tarc(self.tarc)
self.tcircle.set_tarcs()
self.tcircle.plot_tree(self.tarc.arc_id, rlim=(0,550), 
             clade_dict=clade_dict, 
                  linewidth=0.4, linecolor="#606060")
self.tcircle.figure.savefig("filename.pdf")

def highLight(self,species,color, fontsize,alpha):
    """_summary_

    Args:
        species (_type_): _provide the species file_
        color (_type_): _provide the colors file_
        fontsize (_type_): _provide the font size_
        alpha (_type_): _provide the alpha_
    """
    self.fontsize = fontsize
    self.alpha = alpha
with open(species as "rb+") as species:
    highlight_species = species.read()
    species_list = list(re.split("\n"),highlight_species)
    species.close()
with open(color as "rb+") as color:
    color = color.read()
    color_list = list(re.split("\n"),color)
    color.close()
highlight = dict(alpha)
for i in range(len(species)):
    alpha[species_list[i]] = {"color":color[i], 
           "alpha":self.apha, "fontsize":self.fontsize, "label":species_list[i]}
           self.tcircle.plot_highlight(self.tarc.arc_id, highlight)
           self.tcircle.figure.savefig("filename_highlight.pdf")