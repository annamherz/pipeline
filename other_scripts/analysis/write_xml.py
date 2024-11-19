from multiprocessing import Process
import random
import sys
from scipy.stats import sem as sem
import glob
import networkx as nx
import logging
from pipeline import *
from pipeline.utils import validate
from pipeline.analysis import *
from functools import reduce
import scipy.stats as _stats
import numpy as np
import seaborn as sns
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

# analysis paper
# import libraries
# sys.path.insert(1, "/home/anna/Documents/code/python/pipeline")


logging.getLogger().setLevel(logging.ERROR)


print(BSS.__file__)


# define the analysis method to use
ana_dicts = {
    "plain": {
        "estimator": "MBAR",
        "method": "alchemlyb",
        "check overlap": False,
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
}

# set the variables

network_dict = {}

for network in ["combined"]:  # lomap rbfenn combined
    # all the options
    ana_obj_dict = {}

    for protein in ["tyk2", "mcl1", "hif2a", "syk", "p38", "cmet"]:
        ana_obj_dict[protein] = {}

        for ana_dict in ana_dicts.items():
            ana_prot = analysis_protocol(ana_dict[1])

            if protein == "syk" or protein == "cmet":
                main_dir = f"/backup/{protein}/neutral"
            else:
                main_dir = f"/backup/{protein}"

            bench_folder = f"/home/anna/Documents/benchmark"

            # if need size of protein
            try:
                prot = BSS.IO.readMolecules(
                    [
                        f"{bench_folder}/inputs/{protein}/{protein}_prep/{protein}.gro",
                        f"{bench_folder}/inputs/{protein}/{protein}_prep/{protein}.top",
                    ]
                )[0]
            except:
                prot = BSS.IO.readMolecules(
                    [
                        f"{bench_folder}/inputs/{protein}/{protein}_parameterised.prm7",
                        f"{bench_folder}/inputs/{protein}/{protein}_parameterised.rst7",
                    ]
                )[0]

            print(f"no of residues in the protein: {prot.nResidues()}")

            # choose location for the files
            if protein == "syk" or protein == "cmet" or protein == "hif2a":
                # the lomap network
                if network == "lomap" or network == "combined":
                    net_file = f"{main_dir}/execution_model/network_all.dat"
                else:
                    ana_obj_dict[protein][ana_dict[0]] = None
                    continue

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

            # add ligands folder
            if os.path.isdir(f"{bench_folder}/inputs/{protein}/ligands"):
                all_analysis_object.add_ligands_folder(
                    f"{bench_folder}/inputs/{protein}/ligands"
                )
            else:
                all_analysis_object.add_ligands_folder(
                    f"{bench_folder}/inputs/{protein}/ligands_neutral"
                )

            ana_obj_dict[protein][ana_dict[0]] = all_analysis_object

    print(ana_obj_dict)

    network_dict[network] = ana_obj_dict

# for prot in ana_obj_dict.keys():
#     print(prot)
#     ana_obj = network_dict["combined"][prot]["subsampling"]

#     for eng in ana_obj.engines:
#         try:
#             def run_ana():
#                 ana_obj.analyse_mbarnet(compute_missing=True,
#                                         write_xml=True, run_xml_py=True,
#                                         use_experimental=True, overwrite=True,
#                                         engines=[eng],
#                                         )
#             p1 = Process(target=run_ana, name='ana mbarnet')
#             p1.start()
#             p1.join(timeout=5)
#             p1.terminate()
#         except Exception as e:
#             print(e)
#             print(f"failed to mbarnet for {prot} {eng}")


df_dict = {}
for prot in ana_obj_dict.keys():
    ana_obj = ana_obj_dict[prot]["subsampling"]
    df = ana_obj.perturbing_atoms_and_overlap()
    df_dict[prot] = df
