#!/usr/bin/python3

import BioSimSpace as BSS
from BioSimSpace import Units as _Units
import os
import sys
from argparse import ArgumentParser

import pipeline

from pipeline.prep import *
from pipeline.utils import *
from pipeline.analysis import *

import pandas as pd
import matplotlib.pyplot as plt


# for each perturbation, get whether or not it is shrinking or growing


def check_arguments(args):
    # pass the checks to the other check functions

    if args.ligands_folder:
        ligands_folder = args.ligands_folder
    else:
        ligands_folder = str(
            input(
                "what is the ligands folder where the ligand files from the network are?: "
            )
        ).strip()

    if args.execution_folder:
        execution_folder = args.execution_folder
    else:
        execution_folder = str(
            input(
                "what is the execution model folder where the results file should go?: "
            )
        ).strip()

    if args.network_file:
        network_file = args.network_file
    else:
        try:
            network_file = validate.file_path(
                f"{execution_folder}/network_combined.dat"
            )
        except:
            network_file = str(input("what is the path to the network file?: ").strip())

    return ligands_folder, execution_folder, network_file


def main():
    # accept all options as arguments
    parser = ArgumentParser(description="featurise the perturbations")
    parser.add_argument(
        "-lf",
        "--ligands_folder",
        type=str,
        default=None,
        help="ligands folder where the ligand sdf files are located",
    )
    parser.add_argument(
        "-ef",
        "--execution_folder",
        type=str,
        default=None,
        help="execution model folder for the system, with the network file if not provided",
    )
    parser.add_argument(
        "-n", "--network_file", type=str, default=None, help="path to network file"
    )
    args = parser.parse_args()

    # check arguments
    print("checking the provided command line arguments...")
    ligands_folder, execution_folder, network_file = check_arguments(args)

    # check all the perturbations in the network
    perts, ligs = get_info_network(net_file=network_file)

    grow_shrink_dict = {}

    for pert in perts:
        lig0_name = pert.split("~")[0]
        lig1_name = pert.split("~")[1]

        try:  # sdf first
            lig0 = BSS.IO.readMolecules(f"{ligands_folder}/{lig0_name}.sdf")[0]
        except:  # mol2 if sdf is not available
            lig0 = BSS.IO.readMolecules(f"{ligands_folder}/{lig0_name}.mol2")[0]
        try:  # sdf first
            lig1 = BSS.IO.readMolecules(f"{ligands_folder}/{lig1_name}.sdf")[0]
        except:  # mol2 if sdf is not available
            lig1 = BSS.IO.readMolecules(f"{ligands_folder}/{lig1_name}.mol2")[0]

        if lig0.nAtoms() < lig1.nAtoms():
            grow_shrink = "grow"
        elif lig0.nAtoms() > lig1.nAtoms():
            grow_shrink = "shrink"
        else:
            grow_shrink = "same"

        print(
            f"lig0 has {lig0.nAtoms()} and lig1 has {lig1.nAtoms()}, so this is a {grow_shrink} perturbation."
        )

        grow_shrink_dict[pert] = grow_shrink

    # save dict to dat file
    df = pd.DataFrame.from_dict(
        grow_shrink_dict, orient="index", columns=["grow/shrink"]
    )
    df.index.name = "pert"
    df.to_csv(f"{execution_folder}/grow_shrink_featurise.dat", sep=",")


if __name__ == "__main__":
    main()
