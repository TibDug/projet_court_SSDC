import sys
import os
import argparse

__authors__ = ("Lara HERRMANN & Thibault DUGAUQUIER")
__contact__ = ("lara.herrma@gmail.com & thibault.dug@gmail.com")
__date__ = "09 / 09 / 2020"


hydrophobe_list = ["LYS","ARG","HIS","GLU","A","V","I","L","M","F","Y","W"]
hydrophile_list = ["P","S","N","Q","D","E","R","K"]

def get_argument(arguments):
    if len(arguments) != 4:  # Verifie le nombre d'arguments en ligne de commande
        sys.exit("Error: 3 parameters needed.")
    if os.path.exists(arguments[1]) != 1:  # Verifie la presence de l'input
        sys.exit("Error: {} does not exist.".format(arguments[1]))
    if arguments[1].split(".")[1] != "pdb":  # Verifie l'extension de l'input
        sys.exit("Error: {} is not a .pdb file.".format(arguments[1]))
    return arguments[1], arguments[2], int(arguments[3])

def calcul_centre_masse(coord_carbones_alphas_dict):
	xmean, ymean, zmean = 0, 0, 0
	for CA in coord_carbones_alphas_dict :
		xmean += coord_carbones_alphas_dict[CA][1]
		ymean += coord_carbones_alphas_dict[CA][2]
		zmean += coord_carbones_alphas_dict[CA][3]
	mass_center = [xmean/len(coord_carbones_alphas_dict), ymean/len(coord_carbones_alphas_dict), zmean/len(coord_carbones_alphas_dict)]
	return(mass_center)
	
def carbones_alphas_infos(PDB_file):
    coord_carbones_alphas_dict = {}
    type_dict = {}
    with open(PDB_file, "r") as fichier_proteine:
        id_CA = 0
        for ligne in fichier_proteine:
            liste_calpha = []
            if ligne[12:16].strip() == 'CA':
                coord_carbones_alphas_dict[id_CA] = (ligne[31:38].strip(), ligne[39:46].strip(), ligne[47:54].strip())
                if ligne[77:78].strip() in hydrophobe_list:
                    type_dict[id_CA] = 0
                elif ligne[77:78].strip() in hydrophile_list:
                    type_dict[id_CA] = 1
                else:
                    type_dict[id_CA] = -1
    return coord_carbones_alphas_dict, type_dict

if __name__ == "__main__":
    # Recuperation et traitement des donnees en entree.
    # Recuperation des arguments en entree.
    input_file, output_file, points_number = get_argument(sys.argv)
