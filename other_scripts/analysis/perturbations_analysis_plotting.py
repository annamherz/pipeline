# %%
# analysis paper
# import libraries
import seaborn as sns
import numpy as np
import scipy.stats as _stats
from functools import reduce
from pipeline.analysis import *
from pipeline.utils import validate
from pipeline import *
import logging
import networkx as nx
import glob
from scipy.stats import sem as sem
import sys
# sys.path.insert(1, "/home/anna/Documents/code/python/pipeline")

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)
# warnings.simplefilter(action='ignore', category=SettingWithCopyWarning)

logging.getLogger().setLevel(logging.ERROR)


print(BSS.__file__)

# %%


def check_normal_dist(values):
    # check normally dist
    if len(values) < 50:
        stat, p = _stats.shapiro(values)
    else:
        stat, p = _stats.kstest(values)
    if p < 0.05:
        return True
    else:
        return False


# %%
# define the analysis method to use
ana_dicts = {"plain": {
    "estimator": "MBAR",
    "method": "alchemlyb",
    "check overlap": True,
    "try pickle": True,
    "save pickle": True,
    "auto equilibration": False,
    "statistical inefficiency": False,
    "truncate lower": 0,
    "truncate upper": 100,
    "name": None,
},
    "subsampling": {
    "estimator": "MBAR",
    "method": "alchemlyb",
    "check overlap": True,
    "try pickle": True,
    "save pickle": True,
    "auto equilibration": False,
    "statistical inefficiency": True,
    "truncate lower": 0,
    "truncate upper": 100,
    "name": None,
},
    "1ns": {
    "estimator": "MBAR",
    "method": "alchemlyb",
    "check overlap": True,
    "try pickle": True,
    "save pickle": True,
    "auto equilibration": False,
    "statistical inefficiency": True,
    "truncate lower": 0,
    "truncate upper": 25,
    "name": None,
},
    "2ns": {
    "estimator": "MBAR",
    "method": "alchemlyb",
    "check overlap": True,
    "try pickle": True,
    "save pickle": True,
    "auto equilibration": False,
    "statistical inefficiency": True,
    "truncate lower": 0,
    "truncate upper": 50,
    "name": None,
},
    "3ns": {
    "estimator": "MBAR",
    "method": "alchemlyb",
    "check overlap": True,
    "try pickle": True,
    "save pickle": True,
    "auto equilibration": False,
    "statistical inefficiency": True,
    "truncate lower": 0,
    "truncate upper": 75,
    "name": None,
},
    "autoeq": {
    "estimator": "MBAR",
    "method": "alchemlyb",
    "check overlap": True,
    "try pickle": True,
    "save pickle": True,
    "auto equilibration": True,
    "statistical inefficiency": True,
    "truncate lower": 0,
    "truncate upper": 100,
    "name": None,
},
    "TI": {
    "estimator": "TI",
    "method": "alchemlyb",
    "check overlap": True,
    "try pickle": True,
    "save pickle": True,
    "auto equilibration": False,
    "statistical inefficiency": True,
    "truncate lower": 0,
    "truncate upper": 100,
    "name": None,
},
    "single": {
    "estimator": "MBAR",
    "method": "alchemlyb",
    "check overlap": True,
    "try pickle": True,
    "save pickle": True,
    "auto equilibration": False,
    "statistical inefficiency": False,
    "truncate lower": 0,
    "truncate upper": 100,
    "name": None,
}
}

# %%
# set the variables
network = "combined"  # lomap rbfenn combined

# all the options
ana_obj_dict = {}

for protein in ["tyk2", "mcl1", "hif2a", "syk", "p38", "cmet"]:

    ana_obj_dict[protein] = {}

    for ana_dict in ana_dicts:
        ana_prot = analysis_protocol(ana_dicts[ana_dict])
        print(protein, ana_dict)

        if protein == "syk" or protein == "cmet":
            main_dir = f"/backup/{protein}/neutral"
        else:
            main_dir = f"/backup/{protein}"

        bench_folder = f"/home/anna/Documents/benchmark"

        # if need size of protein
        try:
            prot = BSS.IO.readMolecules(
                [f"{bench_folder}/inputs/{protein}/{protein}_prep/{protein}.gro", f"{bench_folder}/inputs/{protein}/{protein}_prep/{protein}.top"])[0]
        except:
            prot = BSS.IO.readMolecules(
                [f"{bench_folder}/inputs/{protein}/{protein}_parameterised.prm7", f"{bench_folder}/inputs/{protein}/{protein}_parameterised.rst7"])[0]

        print(f"no of residues in the protein: {prot.nResidues()}")

        # choose location for the files
        if protein == "syk" or protein == "cmet" or protein == "hif2a":
            # the lomap network
            net_file = f"{main_dir}/execution_model/network_all.dat"
        else:
            net_file = f"{main_dir}/execution_model/network_{network}.dat"

        exp_file = f"{bench_folder}/inputs/experimental/{protein}.yml"
        output_folder = f"{main_dir}/outputs_extracted"

        # prot_file = f"{main_dir}/execution_model/protocol.dat" # no protocol used , name added after if needed
        pipeline_prot = pipeline_protocol(auto_validate=True)
        # pipeline_prot.name("")

        # initialise the network object
        all_analysis_object = analysis_network(
            output_folder,
            exp_file=exp_file,
            net_file=net_file,
            analysis_prot=ana_prot,
            method=pipeline_prot.name(),  # if the protocol had a name
            engines=pipeline_prot.engines(),
        )

        # compute
        all_analysis_object.compute_results()

        if ana_dict == "single":
            all_analysis_object.file_ext = all_analysis_object.file_ext + \
                f"_{ana_dict}"

        # add ligands folder
        if os.path.isdir(f"{bench_folder}/inputs/{protein}/ligands"):
            all_analysis_object.add_ligands_folder(
                f"{bench_folder}/inputs/{protein}/ligands")
        else:
            all_analysis_object.add_ligands_folder(
                f"{bench_folder}/inputs/{protein}/ligands_neutral")

        ana_obj_dict[protein][ana_dict] = all_analysis_object

print(ana_obj_dict)

# %%
# make single vs triplicate results
for prot in ana_obj_dict.keys():

    ana_obj = ana_obj_dict[prot]["single"]
    # function for single dicts
    ana_obj.compute_single_repeat_results()
    for eng in ["AMBER", "SOMD", "GROMACS"]:
        print(prot, eng)
        ana_obj.change_name(eng, f"{eng}_old")
        ana_obj.change_name(f"{eng}_single", eng)
        if eng not in ana_obj.engines:
            ana_obj.engines.append(eng)
        if eng in ana_obj.other_results_names:
            ana_obj.other_results_names.remove(eng)
        print(ana_obj.engines + ana_obj.other_results_names)
        # print(ana_obj.calc_pert_dict[eng])

# # error for a perturbation per single run

# uncertainty_dict_single = {}

# for eng in all_analysis_object.engines:
#     uncertainty_dict_single[eng] = {}
#     repeat = 0
#     for file in all_analysis_object._results_repeat_files[eng]:
#         uncertainty_dict_single[eng][repeat] = {}
#         calc_diff_dict = make_dict.comp_results(
#             file, all_analysis_object.perturbations, eng, name=None
#         )

#         for pert in calc_diff_dict.keys():
#             uncertainty_dict_single[eng][repeat][pert] = calc_diff_dict[pert][1]

#         repeat += 1

# %%
# plot convergence. only plot if has been computed (should have run for not stats ineff)

# for prot in ana_obj_dict.keys():

#     try:
#         ana_obj = ana_obj_dict[prot]["plain"]
#         ana_obj.compute_convergence(
#             compute_missing=True
#         )
#         ana_obj.plot_convergence()
#     except Exception as e:
#         print(e)
#         print(f"could not for {prot}")

# %%
# identify any outliers and plot again if needed above
for prot in ana_obj_dict.keys():
    ana_obj = ana_obj_dict[prot]["plain"]
    print(prot)
    for eng in ana_obj.engines:
        print(
            f"failed percentage for {eng}: {100 - ana_obj.successful_perturbations(eng)[1]} ({len(ana_obj.perturbations) - len(ana_obj.successful_perturbations(eng)[2])} / {len(ana_obj.perturbations)})")
        print(f"{eng} failed perturbations: {ana_obj.failed_perturbations(engine=eng)}")
        print(f"{eng} disconnected ligands: {ana_obj.disconnected_ligands(engine=eng)}")


# %%
# outliers
for prot in ana_obj_dict.keys():
    for name in ana_dicts:
        print(prot, name)
        ana_obj = ana_obj_dict[prot][name]

        for eng in ana_obj.engines:
            print(
                f"outliers {eng}: {ana_obj.get_outliers(threshold=10, name=eng)}")

# exclude outliers
threshold = 10
for prot in ana_obj_dict.keys():

    for name in ana_dicts.keys():
        print(prot, name)
        ana_obj = ana_obj_dict[prot][name]

        for eng in ana_obj.engines:
            ana_obj.remove_outliers(threshold=threshold, name=eng)
        ana_obj.file_ext = ana_obj.file_ext + f"_outliers{threshold}removed"
        # print(ana_obj.file_ext)

# %%
# calcualte the differences in SEM
# SEM differences
sem_dict = {}
sem_dict_name = {}

for name in ana_dicts:

    sem_list_name = []

    for prot in ana_obj_dict.keys():

        sem_dict[prot] = {}
        sem_dict[prot][name] = {}

        ana_obj = ana_obj_dict[prot][name]  # subsampling

        for eng in ana_obj.engines:

            sem_dict[prot][name][eng] = {}

            sem_list = []
            sems = [val[1] for val in ana_obj.calc_pert_dict[eng].values()]
            sem_list.append(sems)
            sem_list_name.append(sems)

            sem_list = reduce(lambda xs, ys: xs + ys, sem_list)
            sem_list = [x for x in sem_list if str(x) != 'nan']

            # if not check_normal_dist(sem_list):
            #     print(f"{prot} {name} not normally dist")

            mean = np.mean(sem_list)
            lower_ci, upper_ci = _stats.norm.interval(confidence=0.95,
                                                      loc=np.mean(sem_list),
                                                      scale=_stats.sem(sem_list))
            print(prot, name, eng, mean, lower_ci, upper_ci)
            sem_dict[prot][name][eng] = (mean, _stats.tstd(
                sem_list), (lower_ci, upper_ci), sem_list)

    sem_list_name = reduce(lambda xs, ys: xs + ys, sem_list_name)
    sem_list_name = [x for x in sem_list_name if str(x) != 'nan']
    mean = np.mean(sem_list_name)
    lower_ci, upper_ci = _stats.norm.interval(confidence=0.95,
                                              loc=np.mean(sem_list_name),
                                              scale=_stats.sem(sem_list_name))
    print(name, mean, lower_ci, upper_ci)
    sem_dict_name[name] = (mean, _stats.tstd(
        sem_list_name), (lower_ci, upper_ci), sem_list_name)


# %%
# plot all the ddG
for prot in ana_obj_dict.keys():
    for method in ana_dicts:
        print(prot, method)
        ana_obj = ana_obj_dict[prot][method]

        stats_string_all = ""
        try:
            mae = ana_obj.calc_mae_engines(pert_val="pert")
        except:
            pass

        for eng in ana_obj.engines:

            stats_string = ""
            try:

                stats_string += f"{eng} MAE: {mae[0][eng]['experimental']:.2f} +/- {mae[1][eng]['experimental']:.2f} kcal/mol, "
                if method == "single":
                    errors = [val[1]
                              for val in ana_obj.calc_pert_dict[eng].values()]
                    stats_string += f"error: {np.mean(errors):.2f} +/- {_stats.tstd(errors):.2f} kcal/mol\n"
                else:
                    stats_string += f"SEM: {sem_dict[prot][name][eng][0]:.2f} +/- {sem_dict[prot][name][eng][1]:.2f} kcal/mol\n"
            except Exception as e:
                print(e)
                print(f"could not compute for {prot} {name} {eng}")

            try:
                ana_obj.plot_scatter_ddG(
                    engines=eng, suptitle=f"{prot}, {method}\n", title=f"{stats_string}")
                ana_obj.plot_scatter_ddG(engines=eng, use_cinnabar=True)
            except:
                pass
            stats_string_all += stats_string

        try:
            ana_obj.plot_scatter_ddG(
                title=f"{prot}, {method}\n {stats_string_all}", engines=ana_obj.engines)
        except:
            print(f"could not plot {prot} {method}")
