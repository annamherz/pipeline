#!/usr/bin/python3

from pipeline.prep import *
from pipeline.utils import *
from pipeline.analysis import *
import BioSimSpace as BSS
import sys
import os

BSS.setVerbose = True


def reparam(
    main_dir,
    lig_name,
):
    sys_file_path = f"{main_dir}/prep/{lig_name}_sys_equil_solv"
    lig_file_path = f"{main_dir}/prep/{lig_name}_lig_equil_solv"
    if os.path.exists(f"{sys_file_path}.prm7"):
        if os.path.exists(f"{lig_file_path}.prm7"):
            print(f"Prep files already generated for {lig_name}. Reparameterising.")
            print("Renaming old files:")
            for file in [sys_file_path, lig_file_path]:
                for ext in ["prm7", "rst7"]:
                    os.rename(f"{file}.{ext}", f"{file}_original.{ext}")

        else:
            print("lig and sys files not found!!")
            sys.exit()

    for leg in [sys_file_path, lig_file_path]:
        system_sol = BSS.IO.readMolecules(
            [f"{leg}_original.rst7", f"{leg}_original.prm7"]
        )

        ligand = merge.extract_ligand(system_sol)
        # paramaterise the ligand
        print(f"Parameterising {lig_name}...")
        lig_p = ligprep.lig_paramaterise(
            ligand, validate.lig_ff("sage")
        )  # openff_unconstrained-2.0.0
        print(validate.lig_ff("sage"))

        system_sol.removeMolecules(ligand)
        sys_reparam = lig_p + system_sol

        # saving pre runs
        print(f"Saving solvated for {leg} and {lig_name}")
        BSS.IO.saveMolecules(
            f"{leg}",
            sys_reparam,
            ["PRM7", "RST7"],  # , "PDB"
        )

        print(f"Done for {lig_name}.")


def main():
    reparam(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
