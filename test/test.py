import BioSimSpace as BSS
import sys
import numpy as np

try:
    import pipeline
except:
    print("adding code to the pythonpath...")
    code = '/home/anna/Documents/code/python'
    if code not in sys.path:
        sys.path.insert(1, code)
    import pipeline

from pipeline import *


from pipeline.analysis import *
from pipeline.utils import *

from cinnabar import wrangle,plotting

# folders
protein = "tyk2"
file_ext = "MBAR_alchemlyb_benchmark" # for results files

# define all the folder locations
bench_folder = f"/home/anna/Documents/benchmark"
main_folder = f"{bench_folder}/{protein}_benchmark"
out_folder = f"{main_folder}/outputs"
res_folder = f"/home/anna/Documents/code/test/results"
temp_folder = f"{main_folder}/outputs/results/temp"
exp_folder = f"{bench_folder}/inputs/experimental"

# make folders that may not exist
folder_list = [res_folder, temp_folder]
for fold in folder_list:
    validate.folder_path(fold, create=True)

# files
net_file = f"{main_folder}/execution_model/network_lomap.dat"
exp_file = f"{exp_folder}/{protein}.yml"
exp_file_dat = f"{res_folder}/exp_data_{protein}.dat"

# OUTPUT
file_ext_out = "test" # for how files will be written

# files that will get written
comp_pert_file_name = f"computed_perturbations_average_{file_ext_out}"
cinnabar_file = f"cinnabar_format_{file_ext_out}"
# TODO add file names for saving graphs

res_dir = "/home/anna/Documents/benchmark/tyk2_benchmark/outputs"
res_obj = analysis_engines(res_dir, net_file=net_file, exp_file=exp_file, engines="SOMD")
res_obj.compute()

# # fwf exp data
# print("fwf")
# exp_dicts = res_obj._get_exp_fwf(fwf_path='/home/anna/Documents/september_2022_workshops/freenrgworkflows/networkanalysis/')
# for key in exp_dicts[0]:
#     print(f"{key} : {exp_dicts[0][key][0]}")
#     # dG_kj_mol = dG_kj_mol - np.mean(dG_kj_mol) # nomralises the free energy by subtracting the mean from them
#     # in fwf

# print("self normalised")
# # print(res_obj.exper_val_dict)
# values_exp = []
# for val in res_obj.exper_val_dict.values():
#     values_exp.append(float(val[0]))
# avg_exp = np.mean(values_exp)
# for key in res_obj.exper_val_dict:
#     print(f"{key} : {float(res_obj.exper_val_dict[key][0]) - avg_exp}")

# x_data = np.asarray([node[1]["exp_DG"] for node in res_obj._cinnabar_networks["SOMD"].graph.nodes(data=True)])


# for node in res_obj._cinnabar_networks["SOMD"].graph.nodes(data=True):
#     print(node[1]["name"], node[1]["exp_DG"], ( float(node[1]["exp_DG"]) - np.mean(x_data)) )

# x_data = x_data - np.mean(x_data)
# print(x_data)

# print(res_obj.normalised_exper_val_dict["SOMD"])

# need the transdir, sys argv 1
# folder = "/home/anna/Documents/code/test/outputs/GROMACS/lig_ejm31~lig_ejm42"
# # folder format is $MAINDIRECTORY/outputs/$2/$1

# # # simfile header
# # if "SOMD" in folder:
# #     try:
# #         add_header_simfile(folder)
# #     except:
# #         print(f"could not add the header to simfile in {folder}")

# # extract to output folder
# extract_trajectory_frames(folder, traj_lambdas=["0.0000"])

# results_files = ["/home/anna/Documents/code/test/final_summary_AMBER_MBAR_alchemlyb_benchmark.csv",
#                  "/home/anna/Documents/code/test/final_summary_SOMD_MBAR_alchemlyb_benchmark.csv"]
# engine = "SOMD"
# output_folder = "/home/anna/Documents/code/test"
# net_file = "/home/anna/Documents/benchmark/tyk2_benchmark/execution_model/network_lomap.dat"
# weight_file = "/home/anna/Documents/benchmark/tyk2_benchmark/execution_model/network_lomap_scores.dat"

# perturbations, ligands, mod_results_files = get_info_network(results_files, net_file, extra_options={"engine":engine})
# gra = graph(ligands, perturbations)
# gra.add_weight(weight_file)

# work_dir = ["/home/anna/Documents/code/test/AMBER_extracted/lig_ejm31~lig_ejm42",
#             "/home/anna/Documents/code/test/GROMACS_extracted/lig_ejm31~lig_ejm42",
#             "/home/anna/Documents/code/test/SOMD_extracted/lig_ejm31~lig_ejm42"]

# for dir in work_dir:
#     analysis = pipeline.analysis.analyse(dir)

#     analysis_options = {'estimator': "MBAR", "method":"alchemlyb",
#                         "check_overlap":True,
#                         "try_pickle":True, 'save_pickle':True,
#                         "auto_equilibration": True,
#                         "truncate_percentage": 75,
#                         "truncate_keep":"end"}

#     analysis.set_options(analysis_options)
#     analysis.analyse_all_repeats()
#     analysis.plot_graphs()
#     final_results_folder = f"{dir}/results"
#     write_analysis_file(analysis, final_results_folder)

# file = "/home/anna/Documents/benchmark/tyk2_benchmark/execution_model_rbfenn_test/protocol.dat"
# system = BSS.IO.readMolecules(
#         [f"/home/anna/Documents/benchmark/mcl1_benchmark/prep/lig_23_lig_equil_solv.rst7",
#         f"/home/anna/Documents/benchmark/mcl1_benchmark/prep/lig_23_lig_equil_solv.prm7"])

# from pipeline.prep import *
# from pipeline.utils import *

# lig_1 = "lig_ejm31"
# lig_2 = "lig_ejm42"
# engine_query = "AMBER"
# num_lambda = 11

# # files that were set in the run_all script
# main_dir = "/home/anna/Documents/benchmark/tyk2_benchmark"
# prot_file = file = "/home/anna/Documents/benchmark/tyk2_benchmark/execution_model_rbfenn_test/protocol.dat"
# prep_dir = f"{main_dir}/prep"  # define lig prep location
# workdir = f"/home/anna/Documents/code/test/outputs/{engine_query}/{lig_1}~{lig_2}" # pert dir

# # parse protocol file
# protocol = pipeline_protocol(prot_file) # instantiate the protocol as an object
# protocol.validate() # validate all the input
# protocol.rewrite_protocol() # rewrite protocol file
# # add the number of lambdas and engine to the protocol
# protocol.num_lambda = validate.num_lambda(num_lambda)
# protocol.engine = validate.engine(engine_query)


# # create the system for each the free and the bound leg.
# system_free = None
# system_bound = None


# for name, leg in zip(["lig", "sys"], ["free", "bound"]):
#     # Load equilibrated inputs for both ligands
#     system_1 = BSS.IO.readMolecules(
#         [f"{prep_dir}/{lig_1}_{name}_equil_solv.rst7", f"{prep_dir}/{lig_1}_{name}_equil_solv.prm7"])
#     system_2 = BSS.IO.readMolecules(
#         [f"{prep_dir}/{lig_2}_{name}_equil_solv.rst7", f"{prep_dir}/{lig_2}_{name}_equil_solv.prm7"])

#     print(f"Preparing the {leg} leg...")
#     if leg == "free":
#         system_free = merge.merge_system(system_1, system_2, protocol.engine)
#     if leg == "bound":
#         system_bound = merge.merge_system(system_1, system_2, protocol.engine)

# # instantiate each system as a fepprep class with the protocol
# fepprep = fepprep(system_free, system_bound, protocol)
# fepprep.generate_folders(workdir)





# import BioSimSpace as BSS

# process_dict = {"auto_equilibration":False, 
#                 "statistical_inefficiency":False,
#                 "truncate_percentage":25,
#                 "truncate_keep":"start"}

# pmf_bound, overlap_matrix_bound = BSS.FreeEnergy.Relative.analyse(
#                         f'/home/anna/Documents/benchmark/tyk2_benchmark/outputs/SOMD_extracted/lig_ejm45~lig_ejm53/bound_0',
#                         estimator="MBAR",
#                         method="alchemlyb",
#                         **process_dict
#                         )

# pmf_free, overlap_matrix_free = BSS.FreeEnergy.Relative.analyse(
#                         f'/home/anna/Documents/benchmark/tyk2_benchmark/outputs/SOMD_extracted/lig_ejm45~lig_ejm53/free_0',
#                         estimator="MBAR",
#                         method="alchemlyb",
#                         **process_dict
#                         )

# freenrg_rel = BSS.FreeEnergy.Relative.difference(pmf_bound, pmf_free)
# print(f"start 25 is {freenrg_rel}")

# process_dict = {"auto_equilibration":False, 
#                 "statistical_inefficiency":False,
#                 "truncate_percentage":50,
#                 "truncate_keep":"start"}

# pmf_bound, overlap_matrix_bound = BSS.FreeEnergy.Relative.analyse(
#                         f'/home/anna/Documents/benchmark/tyk2_benchmark/outputs/SOMD_extracted/lig_ejm45~lig_ejm53/bound_0',
#                         estimator="MBAR",
#                         method="alchemlyb",
#                         **process_dict
#                         )

# pmf_free, overlap_matrix_free = BSS.FreeEnergy.Relative.analyse(
#                         f'/home/anna/Documents/benchmark/tyk2_benchmark/outputs/SOMD_extracted/lig_ejm45~lig_ejm53/free_0',
#                         estimator="MBAR",
#                         method="alchemlyb",
#                         **process_dict
#                         )

# freenrg_rel = BSS.FreeEnergy.Relative.difference(pmf_bound, pmf_free)
# print(f"start 50 is {freenrg_rel}")

# process_dict = {"auto_equilibration":False, 
#                 "statistical_inefficiency":False,
#                 "truncate_percentage":75,
#                 "truncate_keep":"start"}

# pmf_bound, overlap_matrix_bound = BSS.FreeEnergy.Relative.analyse(
#                         f'/home/anna/Documents/benchmark/tyk2_benchmark/outputs/SOMD_extracted/lig_ejm45~lig_ejm53/bound_0',
#                         estimator="MBAR",
#                         method="alchemlyb",
#                         **process_dict
#                         )

# pmf_free, overlap_matrix_free = BSS.FreeEnergy.Relative.analyse(
#                         f'/home/anna/Documents/benchmark/tyk2_benchmark/outputs/SOMD_extracted/lig_ejm45~lig_ejm53/free_0',
#                         estimator="MBAR",
#                         method="alchemlyb",
#                         **process_dict
#                         )

# freenrg_rel = BSS.FreeEnergy.Relative.difference(pmf_bound, pmf_free)
# print(f"start 75 is {freenrg_rel}")


# process_dict = {"auto_equilibration":False, 
#                 "statistical_inefficiency":False,
#                 "truncate_percentage":25,
#                 "truncate_keep":"start"}

# pmf_bound, overlap_matrix_bound = BSS.FreeEnergy.Relative.analyse(
#                         f'/home/anna/Documents/benchmark/tyk2_benchmark/outputs/SOMD_extracted/lig_ejm45~lig_ejm53/bound_0',
#                         estimator="MBAR",
#                         method="alchemlyb",
#                         **process_dict
#                         )

# pmf_free, overlap_matrix_free = BSS.FreeEnergy.Relative.analyse(
#                         f'/home/anna/Documents/benchmark/tyk2_benchmark/outputs/SOMD_extracted/lig_ejm45~lig_ejm53/free_0',
#                         estimator="MBAR",
#                         method="alchemlyb",
#                         **process_dict
#                         )

# freenrg_rel = BSS.FreeEnergy.Relative.difference(pmf_bound, pmf_free)
# print(f"end 25 is {freenrg_rel}")

# process_dict = {"auto_equilibration":False, 
#                 "statistical_inefficiency":False,
#                 "truncate_percentage":50,
#                 "truncate_keep":"start"}

# pmf_bound, overlap_matrix_bound = BSS.FreeEnergy.Relative.analyse(
#                         f'/home/anna/Documents/benchmark/tyk2_benchmark/outputs/SOMD_extracted/lig_ejm45~lig_ejm53/bound_0',
#                         estimator="MBAR",
#                         method="alchemlyb",
#                         **process_dict
#                         )

# pmf_free, overlap_matrix_free = BSS.FreeEnergy.Relative.analyse(
#                         f'/home/anna/Documents/benchmark/tyk2_benchmark/outputs/SOMD_extracted/lig_ejm45~lig_ejm53/free_0',
#                         estimator="MBAR",
#                         method="alchemlyb",
#                         **process_dict
#                         )

# freenrg_rel = BSS.FreeEnergy.Relative.difference(pmf_bound, pmf_free)
# print(f"end 50 is {freenrg_rel}")


# process_dict = {"auto_equilibration":False, 
#                 "statistical_inefficiency":False,
#                 "truncate_percentage":75,
#                 "truncate_keep":"start"}

# pmf_bound, overlap_matrix_bound = BSS.FreeEnergy.Relative.analyse(
#                         f'/home/anna/Documents/benchmark/tyk2_benchmark/outputs/SOMD_extracted/lig_ejm45~lig_ejm53/bound_0',
#                         estimator="MBAR",
#                         method="alchemlyb",
#                         **process_dict
#                         )

# pmf_free, overlap_matrix_free = BSS.FreeEnergy.Relative.analyse(
#                         f'/home/anna/Documents/benchmark/tyk2_benchmark/outputs/SOMD_extracted/lig_ejm45~lig_ejm53/free_0',
#                         estimator="MBAR",
#                         method="alchemlyb",
#                         **process_dict
#                         )

# freenrg_rel = BSS.FreeEnergy.Relative.difference(pmf_bound, pmf_free)
# print(f"end 75 is {freenrg_rel}")

# process_dict = {"auto_equilibration":False, 
#                 "statistical_inefficiency":False,
#                 "truncate_percentage":0}

# pmf_bound, overlap_matrix_bound = BSS.FreeEnergy.Relative.analyse(
#                         f'/home/anna/Documents/benchmark/tyk2_benchmark/outputs/SOMD_extracted/lig_ejm45~lig_ejm53/bound_0',
#                         estimator="MBAR",
#                         method="alchemlyb",
#                         **process_dict
#                         )

# pmf_free, overlap_matrix_free = BSS.FreeEnergy.Relative.analyse(
#                         f'/home/anna/Documents/benchmark/tyk2_benchmark/outputs/SOMD_extracted/lig_ejm45~lig_ejm53/free_0',
#                         estimator="MBAR",
#                         method="alchemlyb",
#                         **process_dict
#                         )

# freenrg_rel = BSS.FreeEnergy.Relative.difference(pmf_bound, pmf_free)
# print(f"all is {freenrg_rel}")