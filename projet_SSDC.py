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
    if len(arguments) != 5:  # Verifie le nombre d'arguments en ligne de commande
        sys.exit("Error: 4 parameters needed.")
    if os.path.exists(arguments[1]) != 1:  # Verifie la presence de l'input
        sys.exit("Error: {} does not exist.".format(arguments[1]))
    if arguments[1].split(".")[1] != "pdb":  # Verifie l'extension de l'input
        sys.exit("Error: {} is not a .pdb file.".format(arguments[1]))
    return arguments[1], arguments[2], int(arguments[3]), float(arguments[4])

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
                    type_dict[id_CA] = -1
                elif ligne[17:20].strip() in hydrophile_list:
                    type_dict[id_CA] = 1
                else:
                    type_dict[id_CA] = 0
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
    
    
def determination_vecteur_normal(mass_center, point):
    vect_normal = (point[0]-mass_center[0], point[1]-mass_center[1], point[2]-mass_center[2])
    return vect_normal
      

def determination_d(vect_normal, point):
    # determination de l'équation cartésienne du plan : 
    # ax + by + cz + d = 0 
    # détermination de d à partir du vecteur 
    d = -(vect_normal[0] * point[0] + vect_normal[1] * point[1] + vect_normal[2] * point[2])
    return d

def distance_residus_plan(vect_normal, coordonees_plan_ref, coordonnees_CA):
    d_plan_reference = determination_d(vect_normal, coordonees_plan_ref)
    distance = abs(vect_normal[0] * coordonnees_CA[0] + vect_normal[1] * coordonnees_CA[1] + vect_normal[2] * coordonnees_CA[2] + d_plan_reference) / math.sqrt(vect_normal[0] ** 2 + vect_normal[1] ** 2+ vect_normal[2] ** 2)
    # print("La distance estde :", distance)
    return distance <= 0.5

# boucle qui englobe tout ça et qui fait "voyager" vers chacun des vecteurs 
# boucle pour parcourir tous les résidus et donne la distance de chaque résidus avec leplanq u'on est entrain de regarder 
# détermination de l'hydrophobicité pour chaque plan et retourner le plan qui est le + hydrophobe 
# boucle qui fait avancer le plan 

############# Est ce qu'on prend TOUS les CA ou seulement ceux dont dssp est élevé ? 
############# En prenant uniquement les CA > seuil 
def calcul_hydrophobicite_plan(ASA_dict, type_dict, seuil_ASA, vecteur_normal, coordonees_plan_ref, coord_carbones_alphas_dict): 
    #Parcours tous les résidus, vérifie qu'ils appartiennent au plan
    #Calcul l'hydrophobicité du plan
    hydrophobicite = 0
    liste_CA = []
    for CA in ASA_dict:
        if ASA_dict[CA] >= seuil_ASA:
            # print("l'ASA est de :", ASA_dict[CA])
            if distance_residus_plan(vecteur_normal, coordonees_plan_ref, coord_carbones_alphas_dict[CA]):
                print("Coordonnees du residus sur le plan : ",coord_carbones_alphas_dict[CA])
                hydrophobicite += type_dict[CA]
                liste_CA.append(CA)
                ############# Est ce qu'on vérifirai pas l'hydrophobicite ici au lieu 
                ############# de tout stocker dans un dico ??? 
    return hydrophobicite, liste_CA

def parcours_plans(ASA_dict, type_dict, seuil_ASA, vecteur_normal, mass_center, coord_carbones_alphas_dict, distance_max_center):
    coordonees_plan_ref = mass_center
    hydrophobicite_max_vecteur = 0
    while(True): 
        print("changement de plan")
        hydrophobicite_plan, liste_CA = calcul_hydrophobicite_plan(ASA_dict, type_dict, seuil_ASA, vecteur_normal, coordonees_plan_ref, coord_carbones_alphas_dict)
        if hydrophobicite_plan > hydrophobicite_max_vecteur:
            hydrophobicite_max_vecteur = hydrophobicite_plan
            liste_CA_max_vecteur = liste_CA
        coordonees_plan_ref = (coordonees_plan_ref[0] + 0.5 * vecteur_normal[0], coordonees_plan_ref[1] + 0.5 * vecteur_normal[1], coordonees_plan_ref[2] + 0.5 * vecteur_normal[2])
        if calcul_distance(coordonees_plan_ref, mass_center) > distance_max_center:
            break
    coordonees_plan_ref = mass_center
    while(True): 
        print("changement de plan")
        hydrophobicite_plan, liste_CA = calcul_hydrophobicite_plan(ASA_dict, type_dict, seuil_ASA, vecteur_normal, coordonees_plan_ref, coord_carbones_alphas_dict)
        if hydrophobicite_plan > hydrophobicite_max_vecteur:
            hydrophobicite_max_vecteur = hydrophobicite_plan
            liste_CA_max_vecteur = liste_CA
        coordonees_plan_ref = (coordonees_plan_ref[0] - 0.5 * vecteur_normal[0], coordonees_plan_ref[1] - 0.5 * vecteur_normal[1], coordonees_plan_ref[2] - 0.5 * vecteur_normal[2])
        if calcul_distance(coordonees_plan_ref, mass_center) > distance_max_center:
            break
    return hydrophobicite_max_vecteur, liste_CA_max_vecteur

def parcours_vecteurs(points, mass_center, coord_carbones_alphas_dict, ASA_dict, type_dict, seuil_ASA, distance_max_center):
    hydrophobicite_max_proteine = 0
    coordonnees_membrane = []
    for point in points:
        print("Coordonnees du point du vecteur: ", point)
        print("Coordonnees du centre de masse: ", mass_center)
        vecteur_normal = determination_vecteur_normal(mass_center, point)
        hydrophobicite_max_vecteur, liste_CA_max_vecteur = parcours_plans(ASA_dict, type_dict, seuil_ASA, vecteur_normal, mass_center, coord_carbones_alphas_dict, distance_max_center)
        if hydrophobicite_max_vecteur > hydrophobicite_max_proteine :
            hydrophobicite_max_proteine = hydrophobicite_max_vecteur
            liste_CA_max_proteine = liste_CA_max_vecteur
    return hydrophobicite_max_proteine, liste_CA_max_proteine

def calcul_distance(coords1, coords2):
    return math.sqrt((coords1[0]-coords2[0])**2 + (coords1[1]-coords2[1])**2 + (coords1[2]-coords2[2])**2)

def distance_max_center(coord_carbones_alphas_dict, mass_center):
    distance_max = 0
    for coords in coord_carbones_alphas_dict.values():
        distance = calcul_distance(coords, mass_center)
        if distance > distance_max:
            distance_max = distance
    return distance_max

if __name__ == "__main__":
    # Recuperation et traitement des donnees en entree.
    # Recuperation des arguments en entree.
    PDB_file, output_file, points_number, seuil_ASA = get_argument(sys.argv)
    coord_carbones_alphas_dict, type_dict = carbones_alphas_infos(PDB_file)
    mass_center = calcul_centre_masse(coord_carbones_alphas_dict)
    distance_max_center = distance_max_center(coord_carbones_alphas_dict, mass_center)
    ASA_dict = accessible_surface_area(PDB_file)
    points = fibonacci_sphere(points_number, mass_center)
    hydrophobicite_max, liste_CA_max_proteine = parcours_vecteurs(points, mass_center, coord_carbones_alphas_dict, ASA_dict, type_dict, seuil_ASA, distance_max_center)
    print(hydrophobicite_max, liste_CA_max_proteine)
    # x = [mass_center[0]]
    # y = [mass_center[1]]
    # z = [mass_center[2]]
    # for point in points:
    #     x.append(point[0])
    #     y.append(point[1])
    #     z.append(point[2])
    # # Creating figure 
    # fig = plt.figure(figsize = (10, 7)) 
    # ax = plt.axes(projection ="3d")     
    # # Creating plot 
    # ax.scatter3D(x, y, z, color = "green"); 
    # plt.title("simple 3D scatter plot")
    # plt.show()
