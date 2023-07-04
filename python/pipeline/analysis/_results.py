import os
import itertools as it
import sys
import re

from ..utils import *
from ._network import *
from ._analysis import *
from ._plotting import *
from ._statistics import *
from ._dictionaries import *
from ._convert import *

import cinnabar
from cinnabar import wrangle, plotting, stats

from math import isnan


class analysis_network:
    """class to analyse results files and plot"""

    def __init__(
        self,
        results_directory=None,
        exp_file=None,
        engines=None,
        net_file=None,
        output_folder=None,
        analysis_prot=None,
        method=None,
        extra_options=None,
        verbose=False,
    ):
        """analyses the network for a certain system

        Args:
            results_directory (str, optional): directory where the results are located. Defaults to None.
            exp_file (str, optional): file path to the experimental results file. Defaults to None.
            engines (list, optional): engines to consider (to comapre multiple). Defaults to None.
            net_file (str, optional): file path to the network of perturbations to analyse. Defaults to None.
            output_folder (str, optional): folder path to where to save all the outputs of the analysis. Defaults to None.
            analysis_prot (pipeline.protocol.analysis_protocol, optional): analysis protocol to make file extension to look for for the results. Defaults to None.
            method (str, optional): Method to consider in the method column of the results files. Defaults to None (will consider all results in the file).
            extra_options (dict, optional): extra options (eg temperature). Defaults to None.

        Raises:
            TypeError: analysis ext must be the correct type.
        """

        self.is_verbose(verbose)

        if method:
            self.method = validate.string(method)
        else:
            self.method = None

        # get engines for analysis
        if engines:
            self.engines = validate.engines(engines)
        else:
            self.engines = validate.engines("ALL")

        if not exp_file:
            print(
                "please set an experimental yml file so this can be used, eg using .get_experimental(exp_file). "
            )
            self.exp_file = None
        else:
            self.exp_file = validate.file_path(exp_file)

        if not net_file:
            print(
                "no network file, will use all perturbations found in results files from the results dir."
            )
            self._net_file = None
            self.net_ext = "network"
        else:
            self._net_file = validate.file_path(net_file)
            self.net_ext = validate.string(f"{net_file.split('/')[-1].split('.')[0]}")

        if not analysis_prot:
            self.file_ext = ".+"  # wildcard, all files in the folder included
            self.analysis_options = analysis_protocol(file=None, auto_validate=True)
        else:
            # read in a dictionary or analysis protocol file
            try:
                self.analysis_options = analysis_protocol(
                    analysis_prot, auto_validate=True
                )
                self.file_ext = analyse.file_ext(self.analysis_options)
            except:
                raise TypeError(
                    f"{analysis_prot} analysis protocol must be an analysis protocol file/dictionary"
                )

        if results_directory:
            self._results_directory = validate.folder_path(results_directory)
            # get files from results directory
            self._results_repeat_files = self._get_results_repeat_files()
            self._results_free_repeat_files = self._get_results_repeat_files(leg="free")
            self._results_bound_repeat_files = self._get_results_repeat_files(
                leg="bound"
            )
            self._results_files = self._get_results_files()
        else:
            print(
                "There is no provided results directory. There are no results to analyse. This will probably create an issue for many other functions. please reinstantiate the object with a results directory."
            )
            self._results_directory = None
            self._results_repeat_files = None
            self._results_free_repeat_files = None
            self._results_bound_repeat_files = None
            self._results_files = None

        if not output_folder:
            if self._results_directory:
                print(
                    "no output folder provided, writing all output to the 'results_directory'."
                )
                self.output_folder = f"{self._results_directory}"
                self.graph_dir = validate.folder_path(
                    f"{self._results_directory}/graphs", create=True
                )
            else:
                print(
                    "no output or results directory, so writing files to current folder..."
                )
                self.output_folder = os.getcwd()
        else:
            self.output_folder = validate.folder_path(output_folder, create=True)
            self.graph_dir = validate.folder_path(
                f"{output_folder}/graphs", create=True
            )

        # set defaults
        self.temperature = 300

        # overwrite if in extra options
        if extra_options:
            extra_options = validate.dictionary(extra_options)

            # temperature for converting experimental values
            if "temperature" in extra_options.keys():
                self.temperature = [validate.is_float(extra_options["temperature"])]

        # get info from the network
        self.perturbations = None
        self.ligands = None
        self._set_network()  # get network info

        # read the extra options
        self._set_default_options()
        if extra_options:
            self.set_options(extra_options)

        # as not yet computed, set this to false
        self._is_computed = False

        # set all the dicts for analysis
        # per engine dicts (also used for other results)
        self.calc_pert_dict = {}  # diff from the results repeat files, average
        self.cinnabar_calc_val_dict = {}  # from the cinnabar network analysis
        self.cinnabar_exper_val_dict = (
            {}
        )  # normalised from the cinnabar network analysis
        self.cinnabar_calc_pert_dict = {}  # from cinnabar network edges
        self.cinnabar_exper_pert_dict = {}  # from cinnabar network edges

        # solo dicts for exper
        self.exper_val_dict = None  # yml converted into experimental values, actual, for ligands in object
        self.normalised_exper_val_dict = (
            None  # yml converted into experimental values, then normalised
        )
        self.exper_pert_dict = None  # yml converted into experimental values, actual, for perturbations in object

        # for convergence
        self.spert_results_dict = {}
        self.spert_bound_dict = {}
        self.spert_free_dict = {}
        self.epert_results_dict = {}
        self.epert_bound_dict = {}
        self.epert_free_dict = {}

        # storing the nx digraphs, per engine
        self._cinnabar_networks = {}
        # overall graph
        self.network_graph = None
        # cycles
        self.cycle_dict = {}

        # for other results
        self.other_results_names = []

        # for checking against free energy workflows
        self._fwf_experimental_DDGs = None
        self._fwf_computed_relative_DDGs = {}
        self._fwf_path = None

        # for plotting
        self._plotting_object = None
        # for stats
        self._stats_object = None

    def is_verbose(self, value):
        verbose = validate.boolean(value)
        self._is_verbose = verbose

        return verbose

    def _get_results_repeat_files(self, leg=None):
        """get the files of all the repeats for a specific leg. Used during init to set free and bound repeat files.

        Args:
            leg (str, optional): Which leg to get the repeats for, ie 'free' or 'bound'. Defaults to None.

        Returns:
            dict: dict of engines as keys and their repeat files for the defined leg.
        """
        res_dir = self._results_directory
        all_files = os.listdir(res_dir)

        files_for_dict = []

        if leg:  # leg should be free or bound
            for file in all_files:
                if f"{leg}_" in file:
                    if re.search(self.file_ext, file):
                        files_for_dict.append(f"{res_dir}/{file}")

        else:  # search for the freenrg
            for file in all_files:
                if "freenrg" in file:
                    files_for_dict.append(f"{res_dir}/{file}")

        files_dict = {}

        for eng in self.engines:
            eng_files = []
            for file in files_for_dict:
                if eng in file:
                    if re.search(self.file_ext, file):
                        eng_files.append(file)
            files_dict[eng] = eng_files

        return files_dict

    def _get_results_files(self):
        """get the summary results files

        Returns:
            dict: dict of engines as keys and their summary files as values.
        """
        res_dir = self._results_directory
        all_files = os.listdir(res_dir)
        sum_files = []
        for file in all_files:
            if "summary" in file:
                if re.search(self.file_ext, file):
                    sum_files.append(f"{res_dir}/{file}")

        files_dict = {}

        for eng in self.engines:
            eng_files = []
            for file in sum_files:
                if eng in file:
                    if re.search(self.file_ext, file):
                        eng_files.append(file)
            files_dict[eng] = eng_files

        return files_dict

    def _set_network(self):
        """set the network based on the network file or based on all the found files."""

        # get perturbations and ligands for the network
        # if network file, get from that
        # if not, use all results files for this

        if not self._results_directory:
            if not self._net_file:
                print(
                    "As there is no provided results directory or network file, please set perturbations and ligands manually."
                )
                return
            else:
                file_names = None

        else:
            # get results files from dict into a list, flatten the list
            results_lists = list(self._results_files.values()) + list(
                self._results_repeat_files.values()
            )
            file_names = [res_file for sublist in results_lists for res_file in sublist]

        values = get_info_network(
            results_files=file_names,
            net_file=self._net_file,
            extra_options={"engines": self.engines},
        )

        self.perturbations = values[0]
        self.ligands = values[1]

    def _set_default_options(self):
        """set the default options for the analysis."""

        self._compute_experimental = True
        self._compute_per_ligand = True
        self._compute_cycles = True

        self._save_pickle = True
        self._try_pickle = True

        # TODO set other default self outputs to None

    def set_options(self, options_dict):
        """set the options based on the options dict.

        Args:
            options_dict (dict): options dict to adjust default options.
        """

        options_dict = validate.dictionary(options_dict)

        if "compute_experimental" in options_dict:
            self._compute_experimental = validate.boolean(
                options_dict["compute_experimental"]
            )
        if "compute_per_ligand" in options_dict:
            self._compute_per_ligand = validate.boolean(
                options_dict["compute_per_ligand"]
            )
        if "compute_cycle_closure" in options_dict:
            self._compute_cycles = validate.boolean(options_dict["compute_cycles"])
        if "XXXX" in options_dict:
            self._XXXX = validate.boolean(options_dict["XXXX"])

        if "save_pickle" in options_dict:
            self._save_pickle = validate.boolean(options_dict["save_pickle"])
        if "try_pickle" in options_dict:
            self._try_pickle = validate.boolean(options_dict["try_pickle"])

    def get_experimental(self, exp_file=None):
        """get the experimental value dictionaries from a given yml file.

        Args:
            exp_file (str, optional): file path to the experimental file. Defaults to None.

        Raises:
            ValueError: if the provided file is not in the yml format.

        Returns:
            tuple: (experimental value dictionary, normalised experimental value dictionary)
        """

        if not exp_file:
            if not self.exp_file:
                warnings.warn(
                    "need an experimental file to proceed with most of the calculations. please set using self.get_experimental(file)"
                )
            else:
                exp_file = self.exp_file
        else:
            self.exp_file = validate.file_path(exp_file)

        # check if yml, as this is needed for
        if exp_file.split(".")[-1] != "yml":
            raise ValueError("the provided experimental file should be in yml format")

        # TODO more checks for type of data
        exper_val_dict = convert.yml_into_exper_dict(
            exp_file, temp=self.temperature
        )  # this output is in kcal/mol

        # experimental value dict
        new_exper_val_dict = make_dict._exper_from_ligands(exper_val_dict, self.ligands)
        self.exper_val_dict = new_exper_val_dict

        # normalise the experimental values
        normalised_exper_val_dict = make_dict._exper_from_ligands(
            exper_val_dict, self.ligands, normalise=True
        )
        self.normalised_exper_val_dict = normalised_exper_val_dict

        return new_exper_val_dict, normalised_exper_val_dict

    def get_experimental_pert(self, exper_val_dict=None):
        """calculate the experimental relative values (pert) for the network from the provided experimental values

        Args:
            exper_val_dict (dict, optional): experimental value dictionary. Defaults to None.

        Returns:
            dict: dictionary of perturbations and their experimental values.
        """

        if exper_val_dict:
            exper_val_dict = validate.dictionary(exper_val_dict)
        else:
            exper_val_dict = self.exper_val_dict

        pert_dict = make_dict._exper_from_perturbations(
            exper_val_dict, self.perturbations
        )

        self.exper_pert_dict = pert_dict

        return pert_dict

    def remove_perturbations(self, perts):
        """remove perturbations from the network used.

        Args:
            perts (list): list of perturbations to remove.
        """

        perts = validate.is_list(perts, make_list=True)

        for pert in perts:
            self.perturbations.remove(pert)

        if self._is_computed:
            self._compute_dicts(use_cinnabar=True)
        # remove plotting object as needs to be reintialised with new perturbations
        self._plotting_object = None
        self._stats_object = None

    def remove_ligands(self, ligs):
        """remove ligand and assosciated perturbations from the network used.

        Args:
            ligs (list): list of ligands to remove.
        """

        ligs = validate.is_list(ligs, make_list=True)

        for lig in ligs:
            self.ligands.remove(lig)

            for pert in self.perturbations:
                if lig in pert:
                    self.perturbations.remove(pert)

        if self._is_computed:
            self._compute_dicts(use_cinnabar=True)
        # remove plotting object as needs to be reintialised with new perturbations
        self._plotting_object = None
        self._stats_object = None

    def change_name(self, old_name, new_name):
        """change the name of the data. Can be used for engine or other result names. Will update the self.dicts with this new name.

        Args:
            old_name (str): old name to replace
            new_name (str): new name to replace it with
        """

        dict_list = [self.calc_pert_dict, self.cinnabar_calc_val_dict, self.cinnabar_exper_val_dict,self.cinnabar_calc_pert_dict,
                     self.cinnabar_exper_pert_dict, self._cinnabar_networks, self.spert_results_dict, self.spert_bound_dict,
                     self.spert_free_dict, self.epert_results_dict, self.epert_bound_dict, self.epert_free_dict, self._fwf_computed_relative_DDGs]

        for adict in dict_list:
            try:
                adict[new_name] = adict.pop(old_name)
            except:
                print(f"could not rename one of the dicts, as it does not have this key as one of its keys.")

        if old_name in self.other_results_names:
            self.other_results_names.remove(old_name)
            self.other_results_names.append(new_name)
        elif old_name in self.engines:
            self.engines.remove(old_name)
            self.other_results_names.append(new_name)

        # remove plotting object as needs to be reintialised with new name
        self._plotting_object = None
        self._stats_object = None

    def compute(self, cycle_closure=True, statistics=True, use_cinnabar=True):
        """compute the dictionaries for analysis and those passed to the plotting object.

        Args:
            cycle_closure (bool, optional): whether to compute the cycle closures. Defaults to True.
            statistics (bool, optional): whether to calculate the statistics. Defaults to True.
            use_cinnabar (bool, optional): whether to use cinnabar during the analysis. Defaults to True.
        """

        # get all the dictionaries needed for plotting
        self._compute_dicts(use_cinnabar=use_cinnabar)

        # TODO have these entirely seperate?
        # compute the cycle closure
        if cycle_closure:
            self._initialise_graph_object()
            self._compute_cycle_closures()

        if statistics:
            # TODO get statistics
            self.pert_statistics = {}
            self.val_statistics = {}
            # for eng in self.engines:
            # self.pert_statistics.update({eng: self.compute_statistics(pert_val="pert", engine=eng)})
            # self.val_statistics.update({eng: self.compute_statistics(pert_val="val", engine=eng)})

        # initialise plotting and stats objects
        self._initialise_plotting_object(verbose=self._is_verbose)
        self._initialise_stats_object()

        self._is_computed = True

    def _compute_dicts(self, use_cinnabar=True):
        """calculate the perturbation dicts from the previously passed repeat files.
        If use_cinnabar, calculate the the cinnabar network and the computed values.

        Args:
            use_cinnabar (bool, optional): whether to use cinnabar. Defaults to True.
        """

        use_cinnabar = validate.boolean(use_cinnabar)

        # compute the experimental for perturbations
        self.get_experimental()  # get experimental val dict and normalised dict
        self.get_experimental_pert()  # from cinnabar expeirmental diff ? make_dict class

        # for self plotting of per pert
        for eng in (
            self.engines + self.other_results_names
        ):  # other results will only be added after already computed once in function below
            if not self._results_files[eng]:
                files = self._results_repeat_files[eng]
            else:
                files = self._results_files[eng]

            calc_diff_dict = make_dict.comp_results(
                files, self.perturbations, eng, method=self.method
            )  # older method
            self.calc_pert_dict.update({eng: calc_diff_dict})

            if use_cinnabar:
                self._compute_cinnabar_dict(files, eng, method=self.method)

    def _compute_cinnabar_dict(self, files, eng, method=None):
        """compute cinnabar and get the dictionaries from it."""

        perts, ligs = get_info_network_from_dict(self.calc_pert_dict[eng])

        # get the files into cinnabar format for analysis
        cinnabar_file_name = (
            f"{self.output_folder}/cinnabar_{eng}_{self.file_ext}_{self.net_ext}"
        )

        convert.cinnabar_file(
            files,
            self.exper_val_dict,
            cinnabar_file_name,
            perturbations=perts,
            method=method,
        )

        try:
            # compute the per ligand for the network
            network = wrangle.FEMap(f"{cinnabar_file_name}.csv")
            self._cinnabar_networks.update({eng: network})

            # from cinnabar graph
            self.cinnabar_calc_pert_dict.update(
                {eng: make_dict.from_cinnabar_network_edges(network, "calc", perts)}
            )
            self.cinnabar_exper_pert_dict.update(
                {eng: make_dict.from_cinnabar_network_edges(network, "exp", perts)}
            )

            # for self plotting of per ligand
            self.cinnabar_calc_val_dict.update(
                {eng: make_dict.from_cinnabar_network_node(network, "calc")}
            )
            self.cinnabar_exper_val_dict.update(
                {
                    eng: make_dict.from_cinnabar_network_node(
                        network, "exp", normalise=True
                    )
                }
            )

        except Exception as e:
            print(e)
            print(f"could not create cinnabar network for {eng}")

    def compute_other_results(self, file_names=None, name=None, method=None, use_cinnabar=True):
        """compute other results in a similar manner to the engine results.

        Args:
            file_names (list, optional): list of other results. Defaults to None.
            name (str, optional): name of these other results (for files and graphs and identification). Defaults to None.
            method (str, optional): method in the input files to include only. Defaults to None.
            use_cinnabar (bool, optional): whether to use cinnabar. Defaults to True.
        """

        file_names = validate.is_list(file_names, make_list=True)
        for file in file_names:
            validate.file_path(file)

        # add identifier for the other results
        name = validate.string(name)
        if name in self.other_results_names:
            print(f"{name} is already in the other results. please use a different name!")
            return
        else:
            self.other_results_names.append(name)

        # add files to file list
        self._results_repeat_files[name] = file_names

        new_file_path = f"{file_names[0].replace(file_names[0].split('/')[-1], '')[:-1]}/{name}_results_file"

        # for self plotting of per pert
        calc_diff_dict = make_dict.comp_results(
            file_names,
            perturbations=None,
            engine=None,
            name=name,
            method=method,
            output_file=new_file_path,
        )
        self.calc_pert_dict.update({name: calc_diff_dict})
        self._results_files[name] = f"{new_file_path}.csv"

        use_cinnabar = validate.boolean(use_cinnabar)
        if use_cinnabar:
            self._compute_cinnabar_dict(files=f"{new_file_path}.csv", eng=name)

        # initialise plotting and stats objects again so theyre added
        self._initialise_plotting_object(verbose=self._is_verbose)
        self._initialise_stats_object()

    def compute_convergence(self, main_dir):
        main_dir = validate.folder_path(main_dir)

        for engine in self.engines:
            self.spert_results_dict[engine] = {}
            self.spert_bound_dict[engine] = {}
            self.spert_free_dict[engine] = {}
            self.epert_results_dict[engine] = {}
            self.epert_bound_dict[engine] = {}
            self.epert_free_dict[engine] = {}

            for pert in self.perturbations:
                # find correct path, use extracted if it exists
                if self.method:
                    name = f"_{self.method}"
                else:
                    name = ""
                path_to_dir = (
                    f"{main_dir}/outputs_extracted/{engine}/{pert}{name}/pickle"
                )
                try:
                    validate.folder_path(path_to_dir)
                except:
                    try:
                        path_to_dir = (
                            f"{main_dir}/outputs/{engine}_extracted/{pert}{name}/pickle"
                        )
                        validate.folder_path(path_to_dir)
                    except:
                        try:
                            path_to_dir = (
                                f"{main_dir}/outputs/{engine}/{pert}{name}/pickle"
                            )
                            validate.folder_path(path_to_dir)
                        except:
                            print(path_to_dir)
                            path_to_dir = None
                            print(
                                f"cannot find pickle directory for {pert} in {engine}, where the pickles saved in a 'pickle' dir in that perturbation folder?"
                            )

                try:
                    pickle_ext = analyse.pickle_ext(
                        self.analysis_options.dictionary(), pert, engine
                    ).split("truncate")[0]

                    with open(
                        f"{path_to_dir}/spert_results_dict_{pickle_ext}.pickle", "rb"
                    ) as file:
                        sresults_dict = pickle.load(file)
                    with open(
                        f"{path_to_dir}/epert_results_dict_{pickle_ext}.pickle", "rb"
                    ) as file:
                        eresults_dict = pickle.load(file)
                    with open(
                        f"{path_to_dir}/spert_bound_dict_{pickle_ext}.pickle", "rb"
                    ) as file:
                        sbound_dict = pickle.load(file)
                    with open(
                        f"{path_to_dir}/epert_bound_dict_{pickle_ext}.pickle", "rb"
                    ) as file:
                        ebound_dict = pickle.load(file)
                    with open(
                        f"{path_to_dir}/spert_free_dict_{pickle_ext}.pickle", "rb"
                    ) as file:
                        sfree_dict = pickle.load(file)
                    with open(
                        f"{path_to_dir}/epert_free_dict_{pickle_ext}.pickle", "rb"
                    ) as file:
                        efree_dict = pickle.load(file)

                    self.spert_results_dict[engine][pert] = sresults_dict
                    self.spert_bound_dict[engine][pert] = sbound_dict
                    self.spert_free_dict[engine][pert] = sfree_dict
                    self.epert_results_dict[engine][pert] = eresults_dict
                    self.epert_bound_dict[engine][pert] = ebound_dict
                    self.epert_free_dict[engine][pert] = efree_dict

                except:
                    print(
                        f"could not load pickles for {pert} in {engine}. Was it analysed for convergence?"
                    )

    def successful_runs(self, eng, perts=None):
        """calculate how many successful runs

        Args:
            eng (str): the engine to calc for
            perts (list, optional): The list of perts to use. Defaults to None.

        Returns:
            tuple: (val, percen, perturbations)
        """

        res_dict = self.calc_pert_dict[eng]
        eng = validate.engine(eng)
        if perts:
            perts = validate.is_list(perts)
        else:
            perts = res_dict.keys()

        perturbations = []
        if self._is_computed:
            val = 0
            for key in res_dict.keys():
                if key in perts:
                    if not isnan(res_dict[key][0]):
                        val += 1
                        perturbations.append(key)

            percen = (val / len(perts)) * 100

            print(f"{val} out of {len(perts)} have results, which is {percen} %.")
            return (val, percen, perturbations)

        else:
            print("please compute results from results files first.")
            return None

    def failed_runs(self, eng):
        eng = validate.engine(eng)

        val, percen, perturbations = self.successful_runs(eng)

        failed_perts = []

        for pert in self.perturbations:
            if pert not in perturbations:
                failed_perts.append(pert)

        return failed_perts

    def disconnected_ligands(self, eng):
        eng = validate.engine(eng)

        self._initialise_graph_object(check=True)
        ligs = self.network_graph.disconnected_ligands()

        return ligs

    def draw_failed_perturbations(self, eng):
        eng = validate.engine(eng)

        val, percen, perturbations = self.successful_runs(eng)

        self._initialise_graph_object(check=True)
        for pert in perturbations:
            self.network_graph.draw_perturbation(pert)

    def remove_outliers(
        self, threshold=10, name=None, verbose=False, use_cinnabar=True
    ):
        """remove outliers above a certain difference to the experimental.

        Args:
            threshold (float, optional): difference threshold above which to remove. Defaults to 10.
            name (str, optional): name of the data (engine or other results). Defaults to None.
            verbose (bool, optional): whether to verbose output which pert is removed. Defaults to False.
            use_cinnabar (bool, optional): whether to use cinnabar. This should be consistent with what was used earlier.
        """

        # can get from dict or dataframe
        # probably best from plotting object

        plot_obj = self._initialise_plotting_object(check=True, verbose=verbose)
        threshold = validate.is_float(threshold)

        perts = []

        if name:
            names = [plot_obj._validate_in_names_list(name)]
        else:
            names = self.other_results_names + self.engines

        for name in names:
            freenrg_df_plotting = plot_obj.freenrg_df_dict["experimental"][name][
                "pert"
            ].dropna()
            x = freenrg_df_plotting[f"freenrg_experimental"]
            y = freenrg_df_plotting["freenrg_calc"]
            # get an array of the MUE values comparing experimental and FEP values. Take the absolute values.
            mue_values = abs(x - y)

            # find the n ligand names that are outliers.
            perts = [
                key
                for pert, key in zip(
                    mue_values.gt(threshold), mue_values.gt(threshold).keys()
                )
                if pert
            ]

            for pert in perts:
                del self.calc_pert_dict[name][pert]

                if verbose:
                    print(f"removed {pert} for {name}.")

            self._compute_dicts(use_cinnabar=use_cinnabar, recompute=True)

        # remove plotting object as needs to be reintialised with new perturbations
        self._plotting_object = None
        self._stats_object = None

    def _initialise_graph_object(self, check=False):
        """intialise the graph object

        Args:
            check (bool, optional): whether to check the plotting object. Defaults to False.

        Returns:
            pipeline.analysis.plotting_engines: the plotting object.
        """

        # if not checking, always make
        if not check:
            self.network_graph = net_graph(self.ligands, self.perturbations)

        # if checking, first see if it exists and if not make
        elif check:
            if not self.network_graph:
                self.network_graph = net_graph(self.ligands, self.perturbations)

        return self.network_graph

    def draw_graph(self, output_dir=None, use_cinnabar=False, engine=None):
        """draw the network graph.

        Args:
            output_dir (str, optional): folder to save the image in. Defaults to None.
            use_cinnabar (bool, optional): whether to use the cinnabar data or the self computed data. Defaults to False.
            engine (str, optional): engine to draw the network for. Defaults to None, draws for each engine.
        """

        if use_cinnabar:
            if engine:
                engines = [engine]
            else:
                engines = self.engines

            for eng in engines:
                if output_dir:
                    file_name = f"{output_dir}/cinnabar_network_{eng}_{self.file_ext}_{self.net_ext}.png"
                else:
                    file_name = None
                self._cinnabar_networks[eng].draw_graph(file_name=file_name)

        else:
            self._initialise_graph_object(check=True)

            self.network_graph.draw_graph(file_dir=output_dir)

        # TODO also incl cinnabar graph drawing functionality?
        # TODO eng specific graph drawing or this doesnt matter

    def _compute_cycle_closures(self):
        """compute the cycle closures and their stats for each engine for the network.

        Returns:
            dict: self.cycle_dict (eng: cycles_dict, cycle_vals, np.mean(cycle_vals), np.std(cycle_vals) )
        """

        # cycle closures

        self._initialise_graph_object(check=True)

        network_graph = self.network_graph

        for eng in self.engines:
            pert_dict = self.cinnabar_calc_pert_dict[eng]

            cycle_closures = network_graph.cycle_closures()

            cycles = make_dict.cycle_closures(pert_dict, cycle_closures)

            # print(f"{eng} cycle vals is {cycles[1]}")
            # print(f"{eng} cycle mean is {cycles[2]}")
            # print(f"{eng} cycle deviation is {cycles[3]}")

            self.cycle_dict.update(
                {eng: (cycles[0], cycles[1], cycles[2], cycles[3])}
            )  # the cycles dict

        return self.cycle_dict

    def cycle_closures(self):
        # TODO
        # if is computed
        # take data from cycle closure dict and plot
        pass

    def _initialise_plotting_object(self, check=False, verbose=False):
        """intialise the plotting object

        Args:
            check (bool, optional): whether to check the plotting object. Defaults to False.

        Returns:
            pipeline.analysis.plotting_engines: the plotting object.
        """

        # if not checking, always make
        if not check:
            self._plotting_object = plotting_engines(
                analysis_object=self, verbose=verbose
            )

        # if checking, first see if it exists and if not make
        elif check:
            if not self._plotting_object:
                self._plotting_object = plotting_engines(
                    analysis_object=self, verbose=verbose
                )

        return self._plotting_object

    def plot_bar_pert(self, engine=None, **kwargs):
        """plot the bar plot of the perturbations.

        Args:
            engine (str, optional): engine to plot for. Defaults to None, will use all.
        """
        engine = validate.is_list(engine, make_list=True)
        engine.append("experimental")

        plot_obj = self._initialise_plotting_object(
            check=True, verbose=self._is_verbose
        )
        plot_obj.bar(pert_val="pert", names=engine, **kwargs)

    def plot_bar_lig(self, engine=None, **kwargs):
        """plot the bar plot of the values per ligand.

        Args:
            engine (str, optional): engine to plot for. Defaults to None, will use all.
        """
        engine = validate.is_list(engine, make_list=True)
        engine.append("experimental")

        plot_obj = self._initialise_plotting_object(
            check=True, verbose=self._is_verbose
        )
        plot_obj.bar(pert_val="val", names=engine, **kwargs)

    def plot_scatter_pert(self, engine=None, use_cinnabar=False, **kwargs):
        """plot the scatter plot of the perturbations.

        Args:
            engine (str, optional): engine to plot for. Defaults to None, will use all.
            use_cinnabar (bool, optional): whether to plot via cinnabar. Defaults to False.
        """

        if use_cinnabar:
            if engine:
                try:
                    engine = validate.engine(engine)
                    engines = [engine]
                except:
                    print(
                        "for cinnabar plotting, can only have one engine. Please use the engine keyword to define, or leave blank to autoidentify available engines."
                    )
                    return

            else:
                print(
                    "no engine specified for plotting, will plot seperate graphs for each self.engines "
                )
                engines = self.engines

            for eng in engines:
                plotting.plot_DDGs(
                    self._cinnabar_networks[eng].graph,
                    filename=f"{self.graph_dir}/DDGs_{eng}_{self.file_ext}_{self.net_ext}.png",
                    title=f"DDGs for {eng}, {self.net_ext}",
                    **{"figsize": 5},
                )  # with {self.file_ext}

        else:
            plot_obj = self._initialise_plotting_object(
                check=True, verbose=self._is_verbose
            )
            plot_obj.scatter(pert_val="pert", y_names=engine, **kwargs)

    def plot_scatter_lig(self, engine=None, use_cinnabar=False, **kwargs):
        """plot the scatter plot of the values per ligand.

        Args:
            engine (str, optional): engine to plot for. Defaults to None, will use all.
            use_cinnabar (bool, optional): whether to plot via cinnabar. Defaults to False.
        """

        if use_cinnabar:
            if engine:
                try:
                    engine = validate.engine(engine)
                    engines = [engine]
                except:
                    print(
                        "for cinnabar plotting, can only have one engine. Please use the engine keyword to define, or leave blank to autoidentify available engines."
                    )
                    return

            else:
                print(
                    "no engine specified for plotting, will plot seperate graphs for each self.engines "
                )
                engines = self.engines

            for eng in engines:
                plotting.plot_DGs(
                    self._cinnabar_networks[eng].graph,
                    filename=f"{self.graph_dir}/DGs_{eng}_{self.file_ext}_{self.net_ext}.png",
                    title=f"DGs for {eng}, {self.net_ext}",
                    **{"figsize": 5},
                )

        else:
            plot_obj = self._initialise_plotting_object(
                check=True, verbose=self._is_verbose
            )
            plot_obj.scatter(pert_val="val", y_names=engine, **kwargs)

    def plot_eng_vs_eng(self, engine_a=None, engine_b=None, pert_val="pert"):
        """plot scatter plot of engine_a vs engine_b

        Args:
            engine_a (str, optional): engine_a. Defaults to None.
            engine_b (str, optional): engine_b. Defaults to None.
            pert_val (str, optional): whether perturbations 'pert' or values per ligand 'val'. Defaults to "pert".
        """

        plot_obj = self._initialise_plotting_object(
            check=True, verbose=self._is_verbose
        )

        if pert_val == "pert":
            binding = "$\Delta\Delta$G$_{bind}$ (kcal/mol)"
        elif pert_val == "val":
            binding = "$\Delta$G$_{bind}$ (kcal/mol)"
        plotting_dict = {
            "title": f"{engine_a} vs {engine_b}\n for {self.file_ext}, {self.net_ext}",
            "y label": f"{engine_a} " + binding,
            "x label": f"{engine_b} " + binding,
            "key": False,
        }

        plot_obj.scatter(
            pert_val=pert_val, y_names=engine_a, x_name=engine_b, **plotting_dict
        )

    def plot_other_results(self, name=None, engine=None, pert_val=None, outliers=None, **kwargs):
        """plot any other results as a scatter plot

        Args:
            name (str, optional): name of the other results (that was used when it was added). Defaults to None.
            engine (str, optional): engine to plot against other results. Defaults to None.
            pert_val (str, optional): whether perturbations 'pert' or values per ligand 'val'. Defaults to None.
            outliers (int, optional): number of outliers to identify. Defaults to None.

        Raises:
            NameError: if the name does not exist in the other results
        """

        name = validate.string(name)

        # first, check if name exists in other dict
        if name in self.other_results_names:
            pass
        else:
            raise NameError(f"{name} does not exist as an added other results.")

        plot_obj = self._initialise_plotting_object(
            check=True, verbose=self._is_verbose
        )

        plot_obj.scatter(pert_val=pert_val, y_names=engine, x_name=name, outliers=outliers, **kwargs)

    def plot_outliers(self, engine=None, no_outliers=5, pert_val="pert", **kwargs):
        """plot scatter plot with annotated outliers.

        Args:
            engine (list, optional): engine to plot for. Defaults to None.
            no_outliers (int, optional): number of outliers to annotate. Defaults to 5.
            pert_val (str, optional): whether plotting 'pert' ie perturbations or 'val' ie values (per ligand result). Defaults to None.

        """

        plot_obj = self._initialise_plotting_object(
            check=True, verbose=self._is_verbose
        )
        plot_obj.scatter(pert_val=pert_val, y_names=engine, no_outliers=no_outliers, **kwargs)

    def plot_histogram_runs(self, engine=None, per_repeat=False):
        """plot histograms for the individual runs (ie per repeats)

        Args:
            engine (str, optional): engine to plot histogram for. Defaults to None.
        """

        plot_obj = self._initialise_plotting_object(
            check=True, verbose=self._is_verbose
        )

        # get list of the errors from the repeat files
        error_dict = {}

        for eng in self.engines:
            error_list = make_dict.error_list_from_files(
                self._results_repeat_files[eng]
            )
            error_dict.update({eng: error_list})

        # if per_repeat:

        plot_obj.histogram(
            engines=engine, error_dict=error_dict, file_ext="individual_runs"
        )

    def plot_histogram_sem(self, engine=None, pert_val="val"):
        """plot histograms for the sem of the result (either pert or val).

        Args:
            engine (str, optional): engine to plot for. Defaults to None.
            pert_val (str, optional): whether perturbations 'pert' or values per ligand 'val'. Defaults to "val".
        """

        plot_obj = self._initialise_plotting_object(
            check=True, verbose=self._is_verbose
        )
        if not engine:
            engine = self.engines
        plot_obj.histogram(engines=engine, pert_val=pert_val)

    # TODO histogram plot for the error in the free and the error in the bound legs

    def plot_histogram_legs(self, engine=None, free_bound=None):
        """plot histograms for the errors of a specific leg (and all repeats for that leg)

        Args:
            engine (str, optional): engine to plot for. Defaults to None.
            free_bound (str, optional): whether to plot for the 'free' or the 'bound' leg. Defaults to None.

        Raises:
            ValueError: free_bound must be 'free' or 'bound'.
        """

        free_bound = validate.string(free_bound)
        if free_bound.lower() not in ["free", "bound"]:
            raise ValueError(
                "free_bound must be either 'free' or 'bound' to denote one of the legs of the run."
            )

        plot_obj = self._initialise_plotting_object(
            check=True, verbose=self._is_verbose
        )

        error_dict = {}

        for eng in self.engines:
            if free_bound == "free":
                error_list = make_dict.error_list_from_files(
                    self._results_free_repeat_files[eng]
                )
            if free_bound == "bound":
                error_list = make_dict.error_list_from_files(
                    self._results_bound_repeat_files[eng]
                )

            error_dict.update({eng: error_list})

        plot_obj.histogram(engines=engine, error_dict=error_dict, file_ext=free_bound)

    def plot_convergence(self, engine=None):
        if not self.spert_results_dict:
            raise EnvironmentError(
                f"please run 'calculate_convergence' first with the main_dir set."
            )

        else:
            if not engine:
                engine = self.engines
            else:
                engine = validate.engines(engine)

            plot_obj = self._initialise_plotting_object(
                check=True, verbose=self._is_verbose
            )
            plot_obj.plot_convergence(engines=engine)

    def _initialise_stats_object(self, check=False):
        """intialise the object for statistical analysis.

        Args:
            check (bool, optional): whether to check. Defaults to False.

        Returns:
            pipeline.analysis.stats_engines: statistics object
        """

        # if not checking, always make
        if not check:
            self._stats_object = stats_engines(analysis_object=self)

        # if checking, first see if it exists and if not make
        elif check:
            if not self._stats_object:
                self._stats_object = stats_engines(analysis_object=self)

        return self._stats_object

    def calc_mae(self, pert_val=None):
        """calculate the Mean Absolute Error (MAE) for between all the engines.
        # TODO expand so can also do for against other results. move to stats?

        Args:
            pert_val (str, optional): whether plotting 'pert' ie perturbations or 'val' ie values (per ligand result). Defaults to None.

        Returns:
            tuple: of dataframe of value and error (mae_pert_df, mae_pert_df_err)
        """

        stats_obj = self._initialise_stats_object(check=True)

        pv = validate.pert_val(pert_val)

        values_dict = stats_obj.values_dict
        engines = self.engines

        mae_pert_df = pd.DataFrame(columns=engines, index=engines)
        mae_pert_df_err = pd.DataFrame(columns=engines, index=engines)

        # iterate over all possible combinations
        for combo in it.product(engines, engines):
            eng1 = combo[0]
            eng2 = combo[1]

            values = stats_obj.compute_mue(pv, x=eng1, y=eng2)
            mean_absolute_error = values[0]  # the computed statistic
            mae_err = values[1]  # the stderr from bootstrapping

            mae_pert_df.loc[eng1, eng2] = mean_absolute_error
            mae_pert_df_err.loc[eng1, eng2] = mae_err

        mae_pert_df.to_csv(
            f"{self.output_folder}/mae_pert_{self.file_ext}.csv", sep=" "
        )
        mae_pert_df_err.to_csv(
            f"{self.output_folder}/mae_pert_err_{self.file_ext}.csv", sep=" "
        )

        return mae_pert_df, mae_pert_df_err

    def compute_mue(self, pert_val=None, engines=None):
        # TODO fix

        stats_obj = self._initialise_stats_object(check=True)
        pert_val = validate.pert_val(pert_val)

        if not engines:
            engines = self.engines
        else:
            engines = validate.engines(engines)

        for eng in engines:
            values = stats_obj.compute_mue(pert_val=pert_val, y=eng)  # TODO fix

        return values

    # TODO compute all stats function
    # TODO write an output file for the dictionary

    # freenergworkflows stuff for comparison
    def _add_fwf_path(self, fwf_path):
        # using freenergworkflows
        if not fwf_path:
            raise ValueError("pls incl the path to freenergworkflows")

        fwf_path = validate.folder_path(fwf_path)

        if fwf_path not in sys.path:
            sys.path.insert(1, fwf_path)

        self._fwf_path = fwf_path

    def _get_exp_fwf(self):
        """get experimental values using freenergworkflows

        Returns:
            tuple: (exp_lig_dict, exp_pert_dict)
        """
        if not self._fwf_path:
            raise ValueError("need fwf path added using _add_fwf_path(fwf_path)")
        import experiments

        # first need to convert the yml file into one useable by freenergworkflows
        exp_file_dat = f"{self.exp_file.split('.')[0]}_exp_dat.dat"
        convert.yml_into_freenrgworkflows(self.exp_file, exp_file_dat)
        experiments = experiments.ExperimentalData()
        experiments.compute_affinities(
            exp_file_dat, data_type="IC50", comments="#", delimiter=","
        )
        experimental_DDGs = experiments.freeEnergiesInKcal

        exp_pert_dict, exp_lig_dict = make_dict.experimental_from_freenrgworkflows(
            experimental_DDGs, self.ligands, self.perturbations
        )
        self._fwf_experimental_DDGs = experimental_DDGs

        return exp_lig_dict, exp_pert_dict

    def _get_ana_fwf(self, engine=None):
        """get experimental values using freenergworkflows

        Args:
            engine (str, optional): name of engine. Defaults to None.

        Raises:
            ValueError: need an engine

        Returns:
            dict: freenerg dict of results for that engine
        """
        # using freenergworkflows
        if not self._fwf_path:
            raise ValueError("need fwf path added using _add_fwf_path(fwf_path)")
        import networkanalysis

        if not engine:
            raise ValueError("please incl an engine")

        # using the network analyser
        nA = networkanalysis.NetworkAnalyser()

        first_file = False
        nf = 0
        for file_name in self._results_repeat_files[engine]:
            # rewrite the file to include only lig_0, lig_1, freenrg, error, engine
            new_file_name = f"{self.output_folder}/fwf_{engine}_file_{nf}.csv"
            data = pd.read_csv(file_name, delimiter=",")
            header_data = data[["lig_0", "lig_1", "freenrg", "error", "engine"]]
            clean_data = header_data.replace("kcal/mol", "", regex=True)
            pd.DataFrame.to_csv(clean_data, new_file_name, sep=",", index=False)
            nf += 1

            if first_file is False:
                nA.read_perturbations_pandas(
                    new_file_name, comments="#", source="lig_0", target="lig_1"
                )
                first_file = True
            else:
                # add more replicates to the graph. FreeNrgWorkflows will take care of averaging
                # the free energies as well as propagating the error.
                nA.add_data_to_graph_pandas(
                    new_file_name, comments="#", source="lig_0", target="lig_1"
                )

        # set network analyser graph as graph
        self.fwf_graph = nA._graph

        computed_relative_DDGs = nA.freeEnergyInKcal

        # this is the per ligand results
        freenrg_dict = make_dict.from_freenrgworkflows_network_analyser(
            computed_relative_DDGs
        )
        self._fwf_computed_relative_DDGs.update({engine: computed_relative_DDGs})

        return freenrg_dict

    def _get_stats_fwf(self, engine=None):
        """get stats using freenergworkflows for the ligands

        Args:
            engine (str, optional): name of engine. Defaults to None.

        Raises:
            ValueError: need an engine

        Returns:
            tuple: r_confidence, tau_confidence, mue_confidence
        """
        # using freenergworkflows
        if not self._fwf_path:
            raise ValueError("need fwf path added using _add_fwf_path(fwf_path)")
        import stats

        if not engine:
            raise ValueError("please incl an engine")

        # these are the per ligand results
        computed_relative_DDGs = self._fwf_computed_relative_DDGs[engine]
        experimental_DDGs = self._fwf_experimental_DDGs

        _stats = stats.freeEnergyStats()
        _stats.generate_statistics(
            computed_relative_DDGs, experimental_DDGs, repeats=10000
        )
        r_confidence = _stats.R_confidence
        tau_confidence = _stats.tau_confidence
        mue_confidence = _stats.mue_confidence
        print(
            "R confidence is:   %.2f < %.2f < %.2f"
            % (r_confidence[1], r_confidence[0], r_confidence[2])
        )
        print(
            "MUE confidence is: %.2f < %.2f < %.2f"
            % (mue_confidence[1], mue_confidence[0], mue_confidence[2])
        )
        print(
            "Tau confidence is: %.2f < %.2f < %.2f"
            % (tau_confidence[1], tau_confidence[0], tau_confidence[2])
        )

        return r_confidence, tau_confidence, mue_confidence
