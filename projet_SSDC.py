import sys
import os
import argparse

__authors__ = ("Lara HERRMANN & Thibault DUGAUQUIER")
__contact__ = ("lara.herrma@gmail.com & thibault.dug@gmail.com")
__date__ = "09 / 09 / 2020"

def get_argument(arguments):
    if len(arguments) != 4:  # Verifie le nombre d'arguments en ligne de commande
        sys.exit("Error: 3 parameters needed.")
    if os.path.exists(arguments[1]) != 1:  # Verifie la presence de l'input
        sys.exit("Error: {} does not exist.".format(arguments[1]))
    if arguments[1].split(".")[1] != "pdb":  # Verifie l'extension de l'input
        sys.exit("Error: {} is not a .pdb file.".format(arguments[1]))
    return arguments[1], arguments[2], int(arguments[3])

if __name__ == "__main__":
    # Recuperation et traitement des donnees en entree.
    # Recuperation des arguments en entree.
    input_file, output_file, points_number = get_argument(sys.argv)