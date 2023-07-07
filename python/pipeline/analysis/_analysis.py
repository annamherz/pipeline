# import libraries
import warnings
import BioSimSpace as BSS
from BioSimSpace import Units as _Units
import pickle
import numpy as _np
import math as _math
import pandas as _pd
from scipy.stats import sem
import pickle
import matplotlib.pyplot as plt

from ..utils import *
from ..prep import analysis_protocol


class analyse:
    """class to analyse a work dir"""

    def __init__(self, work_dir, pert=None, engine=None, analysis_prot=None):
        """class for analysing results in a directory for a single perturbation

        Args:
            work_dir (str): directory of the perturbation results
            pert (str, optional): Name of the perturbation. lig0~lig1 . Defaults to None.
            engine (str, optional): engine to analyse for. Defaults to None.
            analysis_protocol (pipeline.prep.analysis_protocol, optional): the analysis protocol. Defaults to None.

        Raises:
            ValueError: need to have perturbation written as 'ligand_0~ligand_1'
        """
        # instantiate the class with the work directory

        self._work_dir = validate.folder_path(work_dir)
        self._pickle_dir = validate.folder_path(f"{self._work_dir}/pickle", create=True)
        self._graph_dir = validate.folder_path(self._work_dir + "/graphs", create=True)

        if pert:
            self.perturbation = validate.string(pert)
            try:
                self.ligand_0 = self.perturbation.split("~")[0]
                self.ligand_1 = self.perturbation.split("~")[1]
            except:
                raise ValueError(
                    "please seperate the ligands in the pert using '~' so ligand 0 and ligand 1 can be identified."
                )
        else:
            try:
                self.perturbation = self._work_dir.split("/")[-1]
                self.ligand_0 = self.perturbation.split("~")[0]
                self.ligand_1 = self.perturbation.split("~")[1]
            except Exception as e:
                print(e)
                print(
                    "was unable to get the perturbation name and ligands from the file path.\n please add these when initialising using pert='lig_0~lig_1'."
                )

        if engine:
            self.engine = validate.engine(engine)
        else:
            try:
                try:
                    # first try to find in second position
                    self.engine = validate.engine(
                        self._work_dir.split("/")[-2].replace("_extracted", "")
                    )
                except:
                    # then try to find anywhere in the folder path
                    for eng in BSS.FreeEnergy.engines():
                        if eng.upper() in self._work_dir.upper():
                            self.engine = validate.engine(eng)
                            print(f"found {eng.upper()} as engine in work_dir")
            except Exception as e:
                print(e)
                print(
                    "was unable to get the engine from the file path.\n please add when initialising using engine='ENGINE'."
                )

        # initialise dictionaries
        self._bound_pmf_dict = {}  # for the intial results
        self._free_pmf_dict = {}
        self._bound_matrix_dict = {}  # for the energy matrix
        self._free_matrix_dict = {}
        self._bound_val_dict = {}  # for the actual value for all the windows
        self._free_val_dict = {}
        self._bound_err_dict = {}  # for the final error defined
        self._free_err_dict = {}
        self.repeats_tuple_list = []

        # for convergence plotting
        self.spert_results_dict = {}
        self.spert_bound_dict = {}
        self.spert_free_dict = {}
        self.epert_results_dict = {}
        self.epert_bound_dict = {}
        self.epert_free_dict = {}

        # intialise other things
        self._get_repeat_folders()
        self.options_dict = None
        self._set_default_options()  # will set default options and file extension
        self.is_analysed = False

        # initialise values
        self.free_val = None
        self.free_err = None
        self.bound_val = None
        self.bound_err = None
        self.freenrg_val = None
        self.freenrg_err = None

        # list for the successful calculations
        self.bound_calculated = []
        self.free_calculated = []

        if analysis_prot:
            analysis_prot = analysis_protocol(analysis_prot, auto_validate=True)
            self.set_options(analysis_prot)

    @staticmethod
    def _default_analysis_options_dict():
        """the default analysis options dictionary

        Returns:
            dict: default dictionary
        """

        options_dict = {
            "estimator": "MBAR",
            "method": "alchemlyb",
            "check overlap": True,
            "try pickle": True,
            "save pickle": True,
            "auto equilibration": False,
            "statistical inefficiency": False,
            "truncate percentage": 0,
            "truncate keep": "end",
            "mbar method": None,  # robust or default
            "name": None,
        }

        return options_dict

    def _set_default_options(self):
        """set the default options in the class"""

        default_options = analyse._default_analysis_options_dict()
        self.set_options(default_options)

    def _file_ext(self):
        """set the file extension for the analusis files based on the options dict

        Returns:
            str: the file extension
        """

        file_ext = analyse.file_ext(self.options_dict)
        self.file_extension = file_ext

        return file_ext

    @staticmethod
    def file_ext(options_dict):
        """write a file extension based on a protocol style options dictionary

        Args:
            options_dict (dict): analysis protocol dictionary

        Returns:
            str: the file extension
        """

        # validate any inputs in the dictionary
        options_dict = analyse._update_options_dict(options_dict)

        file_ext = str(
            f"{options_dict['estimator']}_{options_dict['method']}_{options_dict['mbar method']}_"
            + f"eq{str(options_dict['auto equilibration']).lower()}_"
            + f"stats{str(options_dict['statistical inefficiency']).lower()}_"
            + f"truncate{str(options_dict['truncate percentage'])}{options_dict['truncate keep']}"
        )

        return file_ext

    def _pickle_ext(self):
        """the extension for the pickle files based on the file extension, engine, and perturbation

        Returns:
            str: pickle extension
        """

        pickle_ext = analyse.pickle_ext(
            self.options_dict, self.perturbation, self.engine
        )
        self.pickle_extension = pickle_ext

        return pickle_ext

    @staticmethod
    def pickle_ext(options_dict, perturbation, engine):
        # validate any inputs in the dictionary
        options_dict = analyse._update_options_dict(options_dict)

        file_ext = analyse.file_ext(options_dict)

        pickle_ext = f"{perturbation}_{engine}_" + f"{file_ext}"

        return pickle_ext

    @staticmethod
    def _update_options_dict(options_dict, current_options=None):
        """update the options dict, if no current one the default will be used.

        Args:
            options_dict (dict or pipeline.prep.analysis_protocol): dictionary of options or the analysis protocol
            current_options (dict, optional): current options. Defaults to None.

        Returns:
            dict: the new options dict that will be used
        """
        # returns a dict with the considered options, and fills in any not provided with default or previous options.

        # if analysis protocol is supplied, make sure to get the dictionary form
        if isinstance(options_dict, analysis_protocol):
            options_dict = options_dict.dictionary()

        options_dict = validate.dictionary(options_dict)

        # use the default dict as the default options
        if not current_options:
            current_options = analyse._default_analysis_options_dict()
        else:
            current_options = validate.dictionary(current_options)

        new_options_dict = {}

        # when initialising, a default dict is set
        # replace any values in this by the new options dict

        # replace any default values by those passed
        for key, value in current_options.items():
            if key in options_dict:
                new_options_dict[key] = options_dict[key]
            else:
                new_options_dict[key] = value

        # validate to make sure all inputs are acceptable format
        val_options_dict = analyse._validate_options_dict(new_options_dict)

        return val_options_dict

    @staticmethod
    def _validate_options_dict(options_dict):
        """validate the passed dictionary

        Args:
            options_dict (dict): analysis options dict

        Returns:
            dict: validated options dictionary.
        """
        # validate all the new inputs
        # replace as needed in the options dict

        if "estimator" in options_dict:
            estimator = validate.estimator(options_dict["estimator"])
            options_dict["estimator"] = estimator

        if "mbar method" in options_dict:
            mbar_method = validate.mbar_method(options_dict["mbar method"])
            options_dict["mbar method"] = mbar_method
            if options_dict["estimator"] == "TI":
                options_dict["mbar method"] = None

        if "check overlap" in options_dict:
            check_overlap = validate.boolean(options_dict["check overlap"])
            if check_overlap == "True" and estimator != "MBAR":
                check_overlap = False
            options_dict["check overlap"] = check_overlap

        if "method" in options_dict:
            method = validate.analysis_method(options_dict["method"])
            options_dict["method"] = method

        if "save pickle" in options_dict:
            save_pickle = validate.boolean(options_dict["save pickle"])
            options_dict["save pickle"] = save_pickle

        if "try pickle" in options_dict:
            try_pickle = validate.boolean(options_dict["try pickle"])
            options_dict["try pickle"] = try_pickle

        if "auto equilibration" in options_dict:
            auto_equilibration = validate.boolean(options_dict["auto equilibration"])
            options_dict["auto equilibration"] = auto_equilibration

        if "statistical inefficiency" in options_dict:
            statistical_inefficiency = validate.boolean(
                options_dict["statistical inefficiency"]
            )
            options_dict["statistical inefficiency"] = statistical_inefficiency

        if "truncate percentage" in options_dict:
            truncate_percentage = validate.integer(options_dict["truncate percentage"])
            options_dict["truncate percentage"] = truncate_percentage

        if "truncate keep" in options_dict:
            truncate_keep = validate.truncate_keep(options_dict["truncate keep"])
            options_dict["truncate keep"] = truncate_keep

        if "name" in options_dict:
            if options_dict["name"]:
                name = validate.string(options_dict["name"])
                options_dict["name"] = name
            else:
                pass

        return options_dict

    def set_options(self, options_dict):
        """set the analysis options for this object

        Args:
            options_dict (dict or pipeline.prep.analysis_protocol): options for analysis
        """

        # first use staticmethod to get a new options dict
        # if already have an options dict for this object, want to use that
        options_dict = analyse._update_options_dict(
            options_dict, current_options=self.options_dict
        )

        # update self.options_dict for use w the file extension
        self.options_dict = options_dict

        # then set all of these to self options
        self.estimator = options_dict["estimator"]
        self._mbar_method = options_dict["mbar method"]
        self._check_overlap = options_dict["check overlap"]
        self.method = options_dict["method"]
        self._save_pickle = options_dict["save pickle"]
        self._try_pickle = options_dict["try pickle"]
        self._auto_equilibration = options_dict["auto equilibration"]
        self._statistical_inefficiency = options_dict["statistical inefficiency"]
        self._truncate_percentage = options_dict["truncate percentage"]
        self._truncate_keep = options_dict["truncate keep"]
        self.name = options_dict["name"]

        # set the file extensions
        self._file_ext()
        self._pickle_ext()

    def _get_repeat_folders(self):
        """how many of each the free and bound repeat folders there are."""

        self._b_folders, self._f_folders = get_repeat_folders(self._work_dir)

        self._no_of_b_repeats = len(self._b_folders)
        self._no_of_f_repeats = len(self._f_folders)
        self._b_repeats = list(range(self._no_of_b_repeats))
        self._f_repeats = list(range(self._no_of_f_repeats))

        if self._no_of_b_repeats != self._no_of_f_repeats:
            print(
                f"There are a different number of repeats for bound ({self._no_of_b_repeats}) and free ({self._no_of_f_repeats}) for {self._work_dir}."
                + f"these are {self._b_folders} and {self._f_folders}."
            )
        else:
            print(
                f"There are {self._no_of_b_repeats} repeats for each the bound and the free for {self._work_dir}."
                + f"these are {self._b_folders} and {self._f_folders}."
            )

    def _check_pickle(self):
        """check if there are all the pickle files present in the pickle folder.

        Returns:
            boolean: whether to try pickle when analysing
        """

        try_pickle = True
        pickle_ext = self.pickle_extension

        # try loading in if previously calculated
        if try_pickle:
            try:
                print(
                    f"trying to locate pickles in default pickle folder, {self._pickle_dir} for {pickle_ext}..."
                )
                with open(
                    f"{self._pickle_dir}/bound_pmf_{pickle_ext}.pickle", "rb"
                ) as file:
                    self._bound_pmf_dict = pickle.load(file)
                with open(
                    f"{self._pickle_dir}/free_pmf_{pickle_ext}.pickle", "rb"
                ) as file:
                    self._free_pmf_dict = pickle.load(file)
                with open(
                    f"{self._pickle_dir}/bound_matrix_{pickle_ext}.pickle", "rb"
                ) as file:
                    self._bound_matrix_dict = pickle.load(file)
                with open(
                    f"{self._pickle_dir}/free_matrix_{pickle_ext}.pickle", "rb"
                ) as file:
                    self._free_matrix_dict = pickle.load(file)
                with open(
                    f"{self._pickle_dir}/bound_val_{pickle_ext}.pickle", "rb"
                ) as file:
                    self._bound_val_dict = pickle.load(file)
                with open(
                    f"{self._pickle_dir}/free_val_{pickle_ext}.pickle", "rb"
                ) as file:
                    self._free_val_dict = pickle.load(file)
                with open(
                    f"{self._pickle_dir}/bound_err_{pickle_ext}.pickle", "rb"
                ) as file:
                    self._bound_err_dict = pickle.load(file)
                with open(
                    f"{self._pickle_dir}/free_err_{pickle_ext}.pickle", "rb"
                ) as file:
                    self._free_err_dict = pickle.load(file)
                print("pickles found!")
            except:
                print("loading pickle failed. Calculating normally.")
                try_pickle = False

        return try_pickle

    def analyse_all_repeats(self):
        """Analyse all existing free-energy data from a simulation working directory.

        Returns:
            tuple: tuple of result, error, repeats tuple list
            (freenrg_rel[0], freenrg_rel[1], repeats_tuple_list)
        """

        if self._try_pickle:
            do_pickle = self._check_pickle()
        else:
            do_pickle = False

        if do_pickle:
            freenrg_rel, repeats_tuple_list = self._analyse_all_repeats_pickle()
        else:
            freenrg_rel, repeats_tuple_list = self._analyse_all_repeats_normal()

        # set the average and error
        self.freenrg = freenrg_rel[0]
        self.error = freenrg_rel[1]
        self.repeats_tuple_list = repeats_tuple_list
        self.is_analysed = True

        if self._check_overlap:
            self.check_overlap()

        if self._save_pickle:
            if do_pickle:
                print("already using pickles, will not be saving again.")
            else:
                print("saving the pmf dictionaries for bound and free as pickles.")
            self.save_pickle()

        return (freenrg_rel[0], freenrg_rel[1], repeats_tuple_list)

    def _analyse_all_repeats_normal(self):
        """Analyse all existing free-energy data from a simulation working directory.

        Returns:
            tuple: (freenrg_rel, repeats_tuple_list)
        """

        # remake list for the successful calculations
        self.bound_calculated = []
        self.free_calculated = []

        process_dict = {
            "auto equilibration": self._auto_equilibration,
            "statistical inefficiency": self._statistical_inefficiency,
            "truncate percentage": self._truncate_percentage,
            "truncate keep": self._truncate_keep,
            "mbar method": self._mbar_method,
        }

        # Analyse the results for each leg of the transformation.
        for b in self._b_repeats:
            try:
                name = str(b) + "_bound"
                pmf_bound, overlap_matrix_bound = BSS.FreeEnergy.Relative.analyse(
                    f"{self._work_dir}/{self._b_folders[b]}",
                    estimator=self.estimator,
                    method=self.method,
                    **process_dict,
                )
                self._bound_pmf_dict.update({name: pmf_bound})
                self._bound_matrix_dict.update({name: overlap_matrix_bound})
                self.bound_calculated.append(name)
            except Exception as e:
                print(e)
                print(
                    f"Unable to analyse values for {name}, which is repeat {self._b_folders[b]} in {self._work_dir}."
                )

        for f in self._f_repeats:
            try:
                name = str(f) + "_free"
                pmf_free, overlap_matrix_free = BSS.FreeEnergy.Relative.analyse(
                    f"{self._work_dir}/{self._f_folders[f]}",
                    estimator=self.estimator,
                    method=self.method,
                    **process_dict,
                )
                self._free_pmf_dict.update({name: pmf_free})
                self._free_matrix_dict.update({name: overlap_matrix_free})
                self.free_calculated.append(name)
            except Exception as e:
                print(e)
                print(
                    f"Unable to analyse values for {name}, which is repeat {self._f_folders[f]} in {self._work_dir}."
                )

        # create a dictionary of the calculated values, for each bound and free
        for r in self._b_repeats:
            try:
                bound_name = str(r) + "_bound"
                bound_val = (
                    self._bound_pmf_dict[bound_name][-1][1]
                    - self._bound_pmf_dict[bound_name][0][1]
                )
                bound_err = BSS.Types.Energy(
                    _math.sqrt(
                        (
                            self._bound_pmf_dict[bound_name][-1][2].value()
                            * self._bound_pmf_dict[bound_name][-1][2].value()
                        )
                        + (
                            self._bound_pmf_dict[bound_name][0][2].value()
                            * self._bound_pmf_dict[bound_name][0][2].value()
                        )
                    ),
                    self._bound_pmf_dict[bound_name][-1][2].unit(),
                )
                self._bound_val_dict.update({bound_name: bound_val})
                self._bound_err_dict.update({bound_name: bound_err})
            except Exception as e:
                print(e)
                print(
                    f"""Unable to compute values for {bound_name} in {self._work_dir}.\
                    Check earlier error message if these values could be analysed."""
                )
        for r in self._f_repeats:
            try:
                free_name = str(r) + "_free"
                free_val = (
                    self._free_pmf_dict[free_name][-1][1]
                    - self._free_pmf_dict[free_name][0][1]
                )
                free_err = BSS.Types.Energy(
                    _math.sqrt(
                        (
                            self._free_pmf_dict[free_name][-1][2].value()
                            * self._free_pmf_dict[free_name][-1][2].value()
                        )
                        + (
                            self._free_pmf_dict[free_name][0][2].value()
                            * self._free_pmf_dict[free_name][0][2].value()
                        )
                    ),
                    self._free_pmf_dict[free_name][-1][2].unit(),
                )
                self._free_val_dict.update({free_name: free_val})
                self._free_err_dict.update({free_name: free_err})
            except Exception as e:
                print(e)
                print(
                    f"""Unable to compute values for {free_name} in {self._work_dir}.\
                    Check earlier error message if these values could be analysed."""
                )

        freenrg_rel, repeats_tuple_list = self._calculate_freenrg()

        return (freenrg_rel, repeats_tuple_list)

    def _analyse_all_repeats_pickle(self):
        """Analyse all existing free-energy data from a simulation working directory in the pickle directory.

        Returns:
            tuple: (freenrg_rel, repeats_tuple_list)
        """

        # remake list for the successful calculations
        self.bound_calculated = []
        self.free_calculated = []

        # Analyse the results for each leg of the transformation.
        for b in self._b_repeats:
            try:
                name = str(b) + "_bound"
                if name in self._bound_pmf_dict.keys():
                    self.bound_calculated.append(name)
            except Exception as e:
                print(e)
                print(
                    f"Unable to analyse values for {name} from pickle, which should be repeat {self._b_folders[b]} in {self._work_dir}."
                )

        for f in self._f_repeats:
            try:
                name = str(f) + "_free"
                if name in self._free_pmf_dict.keys():
                    self.free_calculated.append(name)
            except Exception as e:
                print(e)
                print(
                    f"Unable to analyse values for {name} from pickle, which should be repeat {self._f_folders[f]} in {self._work_dir}."
                )

        freenrg_rel, repeats_tuple_list = self._calculate_freenrg()

        return (freenrg_rel, repeats_tuple_list)

    def _calculate_freenrg(self):
        """calculate the free energy.

        Returns:
            tuple: ((freenrg, freenrg_err), repeats_tuple_list)
        """

        # list for all the repeats
        repeats_tuple_list = []

        # first check if there are calculated values for each.
        if len(self.bound_calculated) == 0:
            print("there are no results for the bound leg.")
            missing_results = True
        elif len(self.free_calculated) == 0:
            print("there are no results for the free leg.")
            missing_results = True
        else:
            missing_results = False

        if missing_results:
            return (np.nan, np.nan), repeats_tuple_list

        # list of values for each leg
        free_vals = [val.value() for val in self._free_val_dict.values()]
        bound_vals = [val.value() for val in self._bound_val_dict.values()]
        free_avg = _np.mean(free_vals)
        bound_avg = _np.mean(bound_vals)
        # errors depending on how many results
        if len(self.free_calculated) == 1:
            free_err = self._free_err_dict[self.free_calculated[0]].value()
        else:
            free_err = sem(free_vals)
        if len(self.bound_calculated) == 1:
            bound_err = self._bound_err_dict[self.bound_calculated[0]].value()
        else:
            bound_err = sem(bound_vals)
        # freenerg values
        freenrg_val = bound_avg - free_avg
        freenrg_err = _math.sqrt(_math.pow(bound_err, 2) + _math.pow(free_err, 2))
        freenrg_rel = (
            freenrg_val * _Units.Energy.kcal_per_mol,
            freenrg_err * _Units.Energy.kcal_per_mol,
        )

        # set the average bound and free values and their error
        self.free_val = free_avg
        self.free_err = free_err
        self.bound_val = bound_avg
        self.bound_err = bound_err
        self.freenrg_val = freenrg_val
        self.freenrg_err = freenrg_err

        # create tuple list of each repeat that was calculated
        # first check the length of the calculated values and check if this is also the length of the folders found
        if len(self.bound_calculated) != self._no_of_b_repeats:
            print(
                "the number of calculated values for bound does not match the number of bound folders.\n maybe try reanalysing/check errors?"
            )
        if len(self.free_calculated) != self._no_of_f_repeats:
            print(
                "the number of calculated values for free does not match the number of free folders.\n maybe try reanalysing/check errors?"
            )

        # if the numebr of calculated values is the same, match these evenly
        if len(self.bound_calculated) == len(self.free_calculated):
            print(
                f"There are {len(self.bound_calculated)} calculated values for each the bound and the free leg for the folders in {self._work_dir}."
            )
            no_of_repeats = len(self.bound_calculated)

        elif len(self.bound_calculated) != len(self.free_calculated):
            print(
                f"There are {len(self.bound_calculated)} calculated values for the bound and {len(self.free_calculated)} calculated values for the free leg for the folders in {self._work_dir}."
            )
            # use the shorter calculated values as the number of complete repeats
            if len(self.bound_calculated) < len(self.free_calculated):
                no_of_repeats = len(self.bound_calculated)
            else:
                no_of_repeats = len(self.free_calculated)

        # calculate
        print(f"{no_of_repeats} repeats will be calculated.")
        r = 0
        for b, f in zip(self.bound_calculated, self.free_calculated):
            print(f"calculating repeat {r} as {b} and {f}.")
            freenrg_rel_rep = BSS.FreeEnergy.Relative.difference(
                self._bound_pmf_dict[b], self._free_pmf_dict[f]
            )
            freenrg_val = freenrg_rel_rep[0]
            freenrg_err = freenrg_rel_rep[1]
            repeats_tuple_list.append((f"{str(r)}_repeat", freenrg_val, freenrg_err))
            r += 1

        return freenrg_rel, repeats_tuple_list

    def save_pickle(self):
        """save the analysis as a pickle."""

        self._pickle_dir = validate.folder_path(self._pickle_dir, create=True)
        pickle_ext = self.pickle_extension

        if not self.is_analysed:
            warnings.warn(
                "can't save pickle, not all repeats have been analysed. please self.analyse_all_repeats() first!"
            )
        else:
            try:
                # write the pmf as a pickle
                with open(
                    f"{self._pickle_dir}/bound_pmf_{pickle_ext}.pickle", "wb"
                ) as handle:
                    pickle.dump(self._bound_pmf_dict, handle)
                with open(
                    f"{self._pickle_dir}/free_pmf_{pickle_ext}.pickle", "wb"
                ) as handle:
                    pickle.dump(self._free_pmf_dict, handle)
                with open(
                    f"{self._pickle_dir}/bound_matrix_{pickle_ext}.pickle", "wb"
                ) as handle:
                    pickle.dump(self._bound_matrix_dict, handle)
                with open(
                    f"{self._pickle_dir}/free_matrix_{pickle_ext}.pickle", "wb"
                ) as handle:
                    pickle.dump(self._free_matrix_dict, handle)
                with open(
                    f"{self._pickle_dir}/bound_val_{pickle_ext}.pickle", "wb"
                ) as handle:
                    pickle.dump(self._bound_val_dict, handle)
                with open(
                    f"{self._pickle_dir}/free_val_{pickle_ext}.pickle", "wb"
                ) as handle:
                    pickle.dump(self._free_val_dict, handle)
                with open(
                    f"{self._pickle_dir}/bound_err_{pickle_ext}.pickle", "wb"
                ) as handle:
                    pickle.dump(self._bound_err_dict, handle)
                with open(
                    f"{self._pickle_dir}/free_err_{pickle_ext}.pickle", "wb"
                ) as handle:
                    pickle.dump(self._free_err_dict, handle)
            except Exception as e:
                print(e)
                print("could not save pickle :( ")

    def check_overlap(self):
        """check the overlap of the analysed object."""

        if not self.is_analysed:
            warnings.warn(
                "can't check overlap, not all repeats have been analysed. please self.analyse_all_repeats() first!"
            )

            return None

        else:
            # check overlap matrix if okay
            okay = 0
            no_overlaps = 0
            too_smalls = 0
            for b in self._b_repeats:
                try:
                    name = str(b) + "_bound"
                    overlap = self._bound_matrix_dict[name]
                    overlap_okay, too_small = BSS.FreeEnergy.Relative.checkOverlap(
                        overlap, estimator=self.estimator
                    )
                    no_overlaps += 1
                    too_smalls += too_small
                    if overlap_okay:
                        okay += 1
                except Exception as e:
                    print(e)
                    print(f"could not check overlap matrix for {name}")

            for f in self._f_repeats:
                try:
                    name = str(f) + "_free"
                    overlap = self._free_matrix_dict[name]
                    overlap_okay, too_small = BSS.FreeEnergy.Relative.checkOverlap(
                        overlap, estimator=self.estimator
                    )
                    no_overlaps += 1
                    too_smalls += too_small
                    if overlap_okay:
                        okay += 1
                except Exception as e:
                    print(e)
                    print(f"could not check overlap matrix for {name}")

            if no_overlaps:
                percen_okay = (okay / no_overlaps) * 100
                too_smalls_avg = too_smalls / no_overlaps
            else:
                percen_okay = None
                too_smalls_avg = None

            return percen_okay, too_smalls_avg

    def plot_graphs(self):
        """plot (overlap matrix or dHdl) for the analysed run."""

        if not self.is_analysed:
            warnings.warn(
                "can't plot, not all repeats have been analysed. please self.analyse_all_repeats() first!"
            )

        else:
            pass

            if self.estimator == "MBAR":
                for b in self._b_repeats:
                    try:
                        name = str(b) + "_bound"
                        overlap = self._bound_matrix_dict[name]
                        ax = BSS.FreeEnergy.Relative.plot(
                            overlap,
                            work_dir=self._graph_dir,
                            file_name=f"{name}_overlap_MBAR_{self.pickle_extension}",
                        )
                    except Exception as e:
                        print(e)
                        print(f"could not plt overlap matrix for {name}")

                for f in self._f_repeats:
                    try:
                        name = str(f) + "_free"
                        overlap = self._free_matrix_dict[name]
                        ax = BSS.FreeEnergy.Relative.plot(
                            overlap,
                            work_dir=self._graph_dir,
                            file_name=f"{name}_overlap_MBAR_{self.pickle_extension}",
                        )
                    except Exception as e:
                        print(e)
                        print(f"could not plt overlap matrix for {name}")

            elif self.estimator == "TI":
                for b in self._b_repeats:
                    name = str(b) + "_bound"
                    overlap = self._bound_matrix_dict[name]
                    try:
                        ax = BSS.FreeEnergy.Relative.plot(
                            overlap,
                            work_dir=self._graph_dir,
                            file_name=f"{name}_dHdl_TI_{self.pickle_extension}",
                        )
                    except Exception as e:
                        print(e)
                        print(f"could not plt dhdl for {name}")

                for f in self._f_repeats:
                    name = str(f) + "_free"
                    overlap = self._free_matrix_dict[name]
                    try:
                        ax = BSS.FreeEnergy.Relative.plot(
                            overlap,
                            work_dir=self._graph_dir,
                            file_name=f"{name}_dHdl_TI_{self.pickle_extension}",
                        )
                    except Exception as e:
                        print(e)
                        print(f"could not plt dhdl for {name}")

    def calculate_convergence(self):
        if self._try_pickle:
            do_pickle = self._check_convergence_pickle()
        else:
            do_pickle = False

        if not do_pickle:
            # calculate all the truncated data
            sresults_dict, sbound_dict, sfree_dict = analyse._calculate_truncated(
                self._work_dir,
                self.estimator,
                "start",
                self._statistical_inefficiency,
                self._auto_equilibration,
            )
            eresults_dict, ebound_dict, efree_dict = analyse._calculate_truncated(
                self._work_dir,
                self.estimator,
                "end",
                self._statistical_inefficiency,
                self._auto_equilibration,
            )

            self.spert_results_dict = sresults_dict
            self.spert_bound_dict = sbound_dict
            self.spert_free_dict = sfree_dict
            self.epert_results_dict = eresults_dict
            self.epert_bound_dict = ebound_dict
            self.epert_free_dict = efree_dict

        self.save_convergence_pickle()

    def _check_convergence_pickle(self):
        """check if there are all the pickle files present in the pickle folder.

        Returns:
            boolean: whether to try pickle when analysing
        """

        do_pickle = False
        pickle_ext = self.pickle_extension.split("truncate")[0]

        # try loading in if previously calculated
        try:
            print(
                f"trying to locate pickles for convergence in default pickle folder, {self._pickle_dir} for {pickle_ext}..."
            )
            with open(
                f"{self._pickle_dir}/spert_results_dict_{pickle_ext}.pickle", "rb"
            ) as file:
                self.spert_results_dict = pickle.load(file)
            with open(
                f"{self._pickle_dir}/epert_results_dict_{pickle_ext}.pickle", "rb"
            ) as file:
                self.epert_results_dict = pickle.load(file)
            with open(
                f"{self._pickle_dir}/spert_bound_dict_{pickle_ext}.pickle", "rb"
            ) as file:
                self.spert_bound_dict = pickle.load(file)
            with open(
                f"{self._pickle_dir}/epert_bound_dict_{pickle_ext}.pickle", "rb"
            ) as file:
                self.epert_bound_dict = pickle.load(file)
            with open(
                f"{self._pickle_dir}/spert_free_dict_{pickle_ext}.pickle", "rb"
            ) as file:
                self.spert_free_dict = pickle.load(file)
            with open(
                f"{self._pickle_dir}/epert_free_dict_{pickle_ext}.pickle", "rb"
            ) as file:
                self.epert_free_dict = pickle.load(file)
            print("pickles found!")
            do_pickle = True
        except:
            print("loading pickle failed. Calculating normally.")

        return do_pickle

    def save_convergence_pickle(self):
        self._pickle_dir = validate.folder_path(self._pickle_dir, create=True)
        pickle_ext = self.pickle_extension.split("truncate")[0]

        # save pickle
        if self._save_pickle:
            print(
                f"trying to save pickles for convergence in default pickle folder, {self._pickle_dir} for {pickle_ext}..."
            )
            # write the pmf as a pickle
            with open(
                f"{self._pickle_dir}/spert_results_dict_{pickle_ext}.pickle", "wb"
            ) as handle:
                pickle.dump(self.spert_results_dict, handle)
            with open(
                f"{self._pickle_dir}/epert_results_dict_{pickle_ext}.pickle", "wb"
            ) as handle:
                pickle.dump(self.epert_results_dict, handle)
            with open(
                f"{self._pickle_dir}/spert_bound_dict_{pickle_ext}.pickle", "wb"
            ) as handle:
                pickle.dump(self.spert_bound_dict, handle)
            with open(
                f"{self._pickle_dir}/epert_bound_dict_{pickle_ext}.pickle", "wb"
            ) as handle:
                pickle.dump(self.epert_bound_dict, handle)
            with open(
                f"{self._pickle_dir}/spert_free_dict_{pickle_ext}.pickle", "wb"
            ) as handle:
                pickle.dump(self.spert_free_dict, handle)
            with open(
                f"{self._pickle_dir}/epert_free_dict_{pickle_ext}.pickle", "wb"
            ) as handle:
                pickle.dump(self.epert_free_dict, handle)
            print("saved pickles!")

    def plot_convergence(self):
        try:
            sdf = analyse.single_pert_dict_into_df(self.spert_results_dict)
            edf = analyse.single_pert_dict_into_df(self.epert_results_dict)
            # plot individually for perts
            print(
                f"plotting for {self.perturbation}, {self.engine} in {self._graph_dir}/forward_reverse_{self.perturbation}.png..."
            )
            analyse.plot_truncated(
                sdf,
                edf,
                file_path=f"{self._graph_dir}/forward_reverse_{self.perturbation}_{self.file_extension.split('truncate')[0]}.png",
                plot_difference=False,
            )
        except Exception as e:
            print(e)
            print("failed to plot convergence, please check Exception message.")

    @staticmethod
    def _calculate_truncated(
        path_to_dir, estimator, start_end="start", statsineff=False, eq=False
    ):
        start_end = validate.truncate_keep(start_end)

        # analyse the work dir
        analysed_pert = analyse(path_to_dir)

        truncate_percentage = [5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100]

        results_dict = {}
        bound_dict = {}
        free_dict = {}

        for trunc_per in truncate_percentage:
            analysed_pert.set_options(
                {
                    "estimator": estimator,
                    "truncate percentage": trunc_per,
                    "try pickle": True,
                    "truncate keep": start_end,
                    "statistical inefficiency": statsineff,
                    "auto equilibration": eq,
                }
            )

            analysed_pert.analyse_all_repeats()

            results_dict[trunc_per] = (
                analysed_pert.freenrg_val,
                analysed_pert.freenrg_err,
            )
            bound_dict[trunc_per] = (analysed_pert.bound_val, analysed_pert.bound_err)
            free_dict[trunc_per] = (analysed_pert.free_val, analysed_pert.free_err)

        return results_dict, bound_dict, free_dict

    @staticmethod
    def single_pert_dict_into_df(pert_dict):
        df = _pd.DataFrame.from_dict(pert_dict)
        df = df.transpose()
        df.columns = ["avg", "max"]
        df["min"] = df.loc[:, "max"] * -1

        return df

    @staticmethod
    def plot_truncated(
        sdf, edf, file_path=None, plot_error=False, plot_difference=True
    ):
        include_key = True

        plt.rc("font", size=12)
        plt.rcParams["axes.xmargin"] = 0  # plt.margins(x=0)
        fig, ax = plt.subplots(figsize=(10, 10))
        lines = []

        # fill in final value and its error
        plt.axhline(y=sdf["avg"].iloc[-1], color="c")
        plt.fill_between(
            sdf.index,
            sdf["avg"].iloc[-1] + sdf["min"].iloc[-1],
            sdf["avg"].iloc[-1] + sdf["max"].iloc[-1],
            color="paleturquoise",
            alpha=0.8,
        )
        lines += plt.plot(0, 0, c="c", label="final estimate")
        scatterplot = [plt.plot(sdf.index, sdf["avg"], c="lightcoral")]
        plt.fill_between(
            sdf.index,
            sdf["avg"] + sdf["min"],
            sdf["avg"] + sdf["max"],
            color="mistyrose",
            alpha=0.4,
        )
        lines += plt.plot(0, 0, c="lightcoral", label="forward")
        scatterplot = [plt.plot(edf.index, edf["avg"], c="cornflowerblue")]
        plt.fill_between(
            edf.index,
            edf["avg"] + edf["min"],
            edf["avg"] + edf["max"],
            color="lightskyblue",
            alpha=0.3,
        )
        lines += plt.plot(0, 0, c="cornflowerblue", label="reverse")

        if include_key:
            labels = [l.get_label() for l in lines]
            plt.legend(lines, labels, loc="upper left")

        plt.xlabel("Percentage of run used")
        if file_path:
            plt.title(f"{file_path.split('/')[-1].split('.')[0].replace('_',' ')}")
        else:
            pass

        if plot_error:
            if plot_difference:
                plt.ylabel("Computed Error difference to final (kcal/mol)")
            else:
                plt.ylabel("Computed Error (kcal/mol)")
        else:
            if plot_difference:
                plt.ylabel(
                    "difference to final result for computed $\Delta\Delta$G$_{bind}$ (kcal/mol)"
                )
            else:
                plt.ylabel("Computed $\Delta\Delta$G$_{bind}$ (kcal/mol)")

        if file_path:
            plt.savefig(file_path)
