# analysis per system

# import libraries
from cinnabar import wrangle as _wrangle
import seaborn as sns
import numpy as np
import scipy.stats as _stats
from functools import reduce
from pipeline.analysis import *
from pipeline.utils import *
from pipeline import *
import logging
import networkx as nx
import glob
from scipy.stats import sem as sem
from matplotlib import colormaps
import sys

from matplotlib.ticker import MaxNLocator

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=RuntimeWarning)
# warnings.simplefilter(action='ignore', category=SettingWithCopyWarning)

logging.getLogger().setLevel(logging.ERROR)


print(BSS.__file__)


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


def flatten_comprehension(matrix):
    return [item for row in matrix for item in row]


# define the analysis method to use
ana_dicts = {
    "plain": {
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
}
protein = sys.argv[1]
network_dict = {}
print(protein)

# "rbfenn", "flare", "combined", "lomap-a-optimal", "lomap-d-optimal", "rbfenn-a-optimal", "rbfenn-d-optimal"
for network in ["combined"]:
    # all the options
    ana_obj_dict = {}

    for ana_dict in ana_dicts.items():
        ana_prot = analysis_protocol(ana_dicts[ana_dict])

        bench_folder = f"/home/anna/Documents/benchmark"
        # main_dir = f"{bench_folder}/reruns/{protein}"
        main_dir = f"/backup/{protein}"

        # choose location for the files
        if protein == "syk" or protein == "cmet" or protein == "hif2a":
            net_file = f"{main_dir}/execution_model/network_lomap.dat"
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
            # method=pipeline_prot.name(),  # if the protocol had a name
            # engines=pipeline_prot.engines(),
        )

        # compute
        try:
            all_analysis_object.compute_results()
        except:
            print("failed analysis")

        # add ligands folder
        all_analysis_object.add_ligands_folder(
            f"{bench_folder}/inputs/reruns/{protein}/ligands_intermediates"
        )

        ana_obj_dict[ana_dict[0]] = all_analysis_object

    network_dict[network] = ana_obj_dict

# set the network for the pertubation analysis
network = "combined"
ana_obj = network_dict[network]["plain"]

# autoequilbration
ana_obj.check_convergence(compute_missing=True)

# convergence
ana_obj.compute_convergence(
    compute_missing=True
)

# mbarnet

# # compute all first
# for eng in ana_obj.engines:
#     # try:
#     ana_obj.analyse_mbarnet(
#         compute_missing=True,
#         write_xml=True,
#         run_xml_py=True,
#         use_experimental=True,
#         overwrite=False,
#         engines=[eng],
#         normalise=True,
#     )
#     # except Exception as e:
#     #     print(e)
#     #     print(f"failed for {eng}")

# for eng in ana_obj.engines:
#     dg_list = []
#     print(eng)

#     try:
#         for key in ana_obj._mbarnet_computed_DGs[eng].keys():
#             value = abs(
#                 ana_obj._mbarnet_computed_DGs[eng][key][0]
#                 - ana_obj.normalised_exper_val_dict[key][0]
#             )
#             dg_list.append(value)
#             if value > 5:
#                 print(eng, key, value)
#     except:
#         pass

# print("done.")
