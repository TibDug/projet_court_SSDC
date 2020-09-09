# HEAD
__authors__ = ("Lara HERRMANN & Thibault DUGAUQUIER")
__contact__ = ("lara.herrma@gmail.com & thibault.dug@gmail.com")
__date__ = "09 / 09 / 2020"

# CODE
hydrophobe_list = ["H","T","C","U","A","V","I","L","M","F","Y","W"]
hydrophile_list = ["P","S","N","Q","D","E","R","K"]

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
