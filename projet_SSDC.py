#ligne de commande : python3 projet_SSDC.py 1t5n.pdb output.txt 3

import sys
import os
import argparse
import math
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

__authors__ = ("Lara HERRMANN & Thibault DUGAUQUIER")
__contact__ = ("lara.herrma@gmail.com & thibault.dug@gmail.com")
__date__ = "09 / 09 / 2020"

hydrophobe_list = ["PHE", "GlY", "ILE", "LEU", "MET", "VAL", "TRP", "TYR"]
hydrophile_list = ["ALA", "CYS", "ASP", "GLU", "HIS", "LYS", "ASN", "PRO", "GLN", "ARG", "SER","THR"]

def get_argument(arguments):
    if len(arguments) != 4:  # Verifie le nombre d'arguments en ligne de commande
        sys.exit("Error: 3 parameters needed.")
    if os.path.exists(arguments[1]) != 1:  # Verifie la presence de l'input
        sys.exit("Error: {} does not exist.".format(arguments[1]))
    if arguments[1].split(".")[1] != "pdb":  # Verifie l'extension de l'input
        sys.exit("Error: {} is not a .pdb file.".format(arguments[1]))
    return arguments[1], arguments[2], int(arguments[3])

def carbones_alphas_infos(PDB_file):
    coord_carbones_alphas_dict = {}
    type_dict = {}
    with open(PDB_file, "r") as fichier_proteine:
        id_CA = 0
        for ligne in fichier_proteine:
            liste_calpha = []
            if ligne[0:6].strip() == 'ATOM' and ligne[12:16].strip() == 'CA':
                coord_carbones_alphas_dict[id_CA] = (float(ligne[30:38].strip()), float(ligne[38:46].strip()), float(ligne[46:54].strip()))
                if ligne[17:20].strip() in hydrophobe_list:
                    type_dict[id_CA] = 0
                elif ligne[17:20].strip() in hydrophile_list:
                    type_dict[id_CA] = 1
                else:
                    type_dict[id_CA] = -1
                id_CA += 1
    return coord_carbones_alphas_dict, type_dict

def calcul_centre_masse(coord_carbones_alphas_dict):
    xmean, ymean, zmean = 0, 0, 0
    for CA in coord_carbones_alphas_dict :
        xmean += coord_carbones_alphas_dict[CA][0]
        ymean += coord_carbones_alphas_dict[CA][1]
        zmean += coord_carbones_alphas_dict[CA][2]
    mass_center = (xmean/len(coord_carbones_alphas_dict), ymean/len(coord_carbones_alphas_dict), zmean/len(coord_carbones_alphas_dict))
    return(mass_center)

def accessible_surface_area(PDB_file):
    ASA_dict = {}
    parser = PDBParser()
    structure_id = PDB_file.split(".")[0]
    structure = parser.get_structure(structure_id, PDB_file)
    model = structure[0]
    dssp = DSSP(model, PDB_file, dssp='mkdssp')
    id_CA = 0
    for CA in list(dssp.keys()):
        if dssp[CA][1] != 'X':
            ASA_dict[id_CA] = dssp[CA][3]
            id_CA += 1
    return ASA_dict


def fibonacci_sphere(samples, mass_center):

    points = []
    phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians

    for i in range(samples):
        y = (1 - (i / float(samples - 1)) * 2) # y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment
        
        y = y + mass_center[1]
        x = (math.cos(theta) * radius) + mass_center[0]
        z = (math.sin(theta) * radius) + mass_center[2]

        points.append((x, y, z))

    return points
      

if __name__ == "__main__":
    # Recuperation et traitement des donnees en entree.
    # Recuperation des arguments en entree.
    PDB_file, output_file, points_number = get_argument(sys.argv)
    coord_carbones_alphas_dict, type_dict = carbones_alphas_infos(PDB_file)
    print(type_dict.values())
    mass_center = calcul_centre_masse(coord_carbones_alphas_dict)
    print(mass_center)
    accessible_surface_area(PDB_file)
    points = fibonacci_sphere(points_number, mass_center)
    print(mass_center)
    print(points)
    
    x = [mass_center[0]]
    y = [mass_center[1]]
    z = [mass_center[2]]
    for point in points:
        x.append(point[0])
        y.append(point[1])
        z.append(point[2])
    # Creating figure 
    fig = plt.figure(figsize = (10, 7)) 
    ax = plt.axes(projection ="3d")     
    # Creating plot 
    ax.scatter3D(x, y, z, color = "green"); 
    plt.title("simple 3D scatter plot")
    plt.show()
