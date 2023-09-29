# import libraries
import BioSimSpace as BSS

import yaml
import csv
import numpy as np

from ..utils import *
from ._dictionaries import *

from typing import Union, Optional

# conversion of constants
from scipy.constants import R, calorie

kJ2kcal = 1 / calorie
R_kJmol = R / 1000
R_kcalmol = R_kJmol * kJ2kcal


class convert:
    """class of static methods for converting data"""

    def __init__(self):
        pass

    @staticmethod
    def yml_into_freenrgworkflows(exp_file: str, exp_file_dat: str):
        """convert yml format into one suitable for freenergworkflows

        Args:
            exp_file (str): yml file. Each key is a ligand, with a 'measurement' that has 'unit', 'value', 'error'. Unit in uM or nM.
            exp_file_dat (str): new file to write experimental data to (fwf format)
        """
        # get the experimental data into a useable format (from yml to csv)
        # for freenergworkflows, want to save as lig, Ki
        # experimental values (e.g. ic50/ki) for all ligands in our set.

        exp_file = validate.file_path(exp_file)
        exp_file_dat = validate.string(exp_file_dat)

        with open(exp_file, "r") as file:
            data = yaml.safe_load(file)  # loads as dictionary

        with open(exp_file_dat, "w") as file:
            writer = csv.writer(file, delimiter=",")
            writer.writerow(["ligand", "value", "error"])

            # the data needs to be IC50, uM
            # am assuming that ki and IC50 are the same

            for key in data.keys():  # write for each ligand that was in yaml file
                if data[key]["measurement"]["unit"] == "uM":
                    writer.writerow(
                        [
                            key,
                            data[key]["measurement"]["value"],
                            data[key]["measurement"]["error"],
                        ]
                    )
                elif data[key]["measurement"]["unit"] == "nM":
                    writer.writerow(
                        [
                            key,
                            "{:.4f}".format(data[key]["measurement"]["value"] / 1000),
                            data[key]["measurement"]["error"] / 1000,
                        ]
                    )

    @staticmethod
    def convert_M_kcal(
        value: float, err: float, temperature: Optional[Union[int, float]] = 300
    ) -> tuple:
        """convert value into kcal/mol

        Args:
            value (float): value
            temperature (int, optional): temperature of the simulation. Defaults to 300.

        Returns:
            string: value in kcal/mol
        """

        # gas constant in kcal per Kelvin per mol, exp val converted into M
        kcal_val = (
            R_kcalmol * temperature * np.log(value)
        )  # value is /1 as the concentration of substrate

        # propagate the error
        kcal_err = abs(R_kcalmol * temperature * (err / (value * np.log(10))))

        return kcal_val, kcal_err

    @staticmethod
    def yml_into_exper_dict(
        exp_file, temperature: Optional[Union[int, float]] = 300
    ) -> dict:
        """convert yml file into experimental dictionary of values.

        Args:
            exp_file (str): yml file. Each key is a ligand, with a 'measurement' that has 'unit', 'value', 'error'. Unit in uM or nM.
            temperature (int, optional): Temperature to use during the conversion. Defaults to 300.

        Returns:
            dict: kcal/mol value for each ligand
        """

        exp_file = validate.file_path(exp_file)
        temperature = validate.is_float(temperature)

        with open(exp_file, "r") as file:
            data = yaml.safe_load(file)  # loads as dictionary

        exper_raw_dict = {}
        for key in data.keys():  # write for each ligand that was in yaml file
            # check what type of data
            if data[key]["measurement"]["type"].lower().strip() == "ki":
                # assuming concentration of substrate in the assay is 1, IC50 is approx Ki*2
                factor = 2
            elif data[key]["measurement"]["type"].lower().strip() == "ic50":
                factor = 1
            else:
                factor = 1
                logging.error(
                    "the type of data was not recognised. Must be IC50 or Ki. Assuming IC50 ..."
                )

            if data[key]["measurement"]["unit"] == "uM":
                magnitude = 10**-6
            elif data[key]["measurement"]["unit"] == "nM":
                magnitude = 10**-9
            else:
                logging.critical("only nM and uM recognised as units!")
                magnitude = 1

            exper_raw_dict[key] = (
                data[key]["measurement"]["value"] * factor * magnitude,
                data[key]["measurement"]["error"] * factor * magnitude,
            )

        exper_val_dict = convert.exper_raw_dict_into_val_dict(
            exper_raw_dict, temperature
        )

        return exper_val_dict

    @staticmethod
    def csv_into_exper_dict(
        exp_file: str, temperature: Optional[Union[int, float]] = 300
    ) -> dict:
        """convert csv file into experimental dictionary of values.

        Args:
            exp_file (str): csv file. Has columns 'ligand', 'unit', 'value', 'error'. Unit in uM or nM.
            temperature (int, optional): Temperature to use during the conversion. Defaults to 300.

        Returns:
            dict: kcal/mol value for each ligand
        """

        exp_file = validate.file_path(exp_file)
        temperature = validate.is_float(temperature)

        res_df = pd.read_csv(exp_file)

        exper_raw_dict = {}

        for index, row in res_df.iterrows():
            if row["type"].lower().strip() == "ki":
                # assuming concentration of substrate in the assay is 1, IC50 is approx Ki*2
                factor = 2
            elif row["type"].lower().strip() == "ic50":
                factor = 1
            else:
                factor = 1
                logging.error(
                    "the type of data was not recognised. Must be IC50 or Ki. Assuming IC50 ..."
                )

            if row["unit"].strip() == "uM":
                magnitude = 10**-6
            elif row["unit"].strip() == "nM":
                magnitude = 10**-9
            else:
                logging.critical("only nM and uM recognised as units!")
                magnitude = 1

            value = float(row["value"].strip()) * factor * magnitude

            try:
                err = (float(row["error"].strip()) * factor * magnitude,)
            except:
                err = None

                exper_raw_dict[row["ligand"]] = (value, err)

        exper_val_dict = convert.exper_raw_dict_into_val_dict(
            exper_raw_dict, temperature
        )

        return exper_val_dict

    @staticmethod
    def exper_raw_dict_into_val_dict(
        exper_raw_dict: dict, temperature: Optional[Union[int, float]] = 300
    ) -> dict:
        """convert raw exp data dict in uM to kcal/mol

        Args:
            exper_raw_dict (dict): the experimental data in uM with {ligand:(value, error)}
            temperature (int, optional): temperature for conversion. Defaults to 300.

        Returns:
            _type_: _description_
        """

        exper_raw_dict = validate.dictionary(exper_raw_dict)
        temperature = validate.is_float(temperature)

        exper_val_dict = {}

        for key in exper_raw_dict.keys():
            lig = str(key)

            exp_val = float(
                exper_raw_dict[key][0]
            )  # this is in IC50 from reading in the files
            err_val = float(exper_raw_dict[key][1])
            # if there is no error, assign based on https://doi.org/10.1016/j.drudis.2009.01.012
            # ie that assay error is usually 0.3 log units standard deviation pIC50 or a factor of 2 in IC50
            if not err_val:  # ie there is no provided error
                # median standard deviation is 0.3
                err_val = 0.3 + np.log(10) * exp_val

            # convert into kcal mol
            kcal_val, kcal_err = convert.convert_M_kcal(
                exp_val, err_val, temperature=temperature
            )

            # add to dict
            exper_val_dict.update({lig: (kcal_val, kcal_err)})

        return exper_val_dict

    @staticmethod
    def cinnabar_file(
        results_files: list,
        exper_val: Union[dict, str],
        output_file: str,
        perturbations: Optional[list] = None,
        method: Optional[str] = None,
    ):
        """convert results files into format needed for cinnabar. If multiple results files, uses the average of a perturbation.

        Args:
            results_files (list): list of results files
            exper_val (dict or str): dict of experimental values or yml file of experimental values.
            output_file (str): output file path
            perturbations (list, optional): list of perturbations to include. Defaults to None.
            method (str, optional): name of the method to consider.
        """
        # files is a list of files
        results_files = validate.is_list(results_files, make_list=True)
        # output file

        if exper_val:
            # first check if the experimental values are a dict or a file
            try:
                exper_val_dict = validate.dictionary(exper_val)
            except:
                validate.file_path(exper_val)
                logging.info("input is a file, will convert this into a dict...")
                logging.info("please check that the conversion of values is okay.")
                exper_val_dict = convert.yml_into_exper_dict(exper_val)

            add_exper_values = True

        else:
            add_exper_values = False

        # write to a csv file
        with open(f"{output_file}.csv", "w") as cinnabar_data_file:
            writer = csv.writer(cinnabar_data_file, delimiter=",")

            if add_exper_values:
                # first, write the experimental data
                writer.writerow(["# Experimental block"])
                writer.writerow(["# Ligand", "expt_DDG", "expt_dDDG"])

                # write in kcal/mol
                for lig in exper_val_dict.keys():
                    writer.writerow(
                        [lig, f"{exper_val_dict[lig][0]}", f"{exper_val_dict[lig][1]}"]
                    )

            # second write the perturbation data
            writer.writerow([" "])
            writer.writerow(["# Calculated block"])
            writer.writerow(
                [
                    "# Ligand1",
                    "Ligand2",
                    "calc_DDG",
                    "calc_dDDG(MBAR)",
                    "calc_dDDG(additional)",
                ]
            )

            # need to write the average of the data, otherwise cinnabar just uses the last entry
            comp_diff_dict = make_dict.comp_results(
                results_files,
                perturbations=perturbations,
                method=method,
            )

            # write to file
            for key in comp_diff_dict:
                lig_0 = key.split("~")[0]
                lig_1 = key.split("~")[1]
                comp_ddG = comp_diff_dict[key][0]
                comp_err = comp_diff_dict[key][1]

                if not comp_ddG:
                    pass
                else:
                    if perturbations:
                        pert = f"{lig_0}~{lig_1}"
                        if pert in perturbations:
                            writer.writerow([lig_0, lig_1, comp_ddG, comp_err, "0.0"])
                        else:
                            pass

                    else:
                        writer.writerow([lig_0, lig_1, comp_ddG, comp_err, "0.0"])
