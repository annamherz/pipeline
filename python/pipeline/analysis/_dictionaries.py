
# import libraries
import BioSimSpace as BSS
import os
import glob
import csv
import numpy as np
import math
import pandas as pd
import networkx as nx
import yaml
from scipy.stats import sem as sem
from scipy.stats import bootstrap
from sklearn.metrics import mean_absolute_error as mae
import pickle
import tempfile
import itertools

import matplotlib.pyplot as plt
import csv
import numpy as np
import pandas as pd 

from ..utils import *

# functions
# TODO clean up and make sure have description at start


# make dictionary of values of results
def make_dict_comp_results(results_files=None, perturbations=None, file_name=None, engine=None):

    # TODO some way to get engine from results files?
    # make a dictionary with the results of the files
    comp_dict_list = {}
    comp_err_dict_list = {}

    # append for results file
    for res_file in results_files:
        res_df = pd.read_csv(res_file)
        for index,row in res_df.iterrows():
            lig_0 = row[0]
            lig_1 = row[1]
            pert = f"{lig_0}~{lig_1}"
            if not isinstance(row[2], float):
                ddG = BSS.Types.Energy(float(row[2].split()[0]),row[2].split()[-1])
            else:
                ddG = row[2]
            if not isinstance(row[3], float):
                ddG_err = BSS.Types.Energy(float(row[3].split()[0]),row[3].split()[-1])
            else:
                ddG_err = row[3]
                
            if pert in comp_dict_list:
                # Key exist in dict, check if is a list
                if not isinstance(comp_dict_list[pert], list):
                    # If type is not list then make it list
                    comp_dict_list[pert] = [comp_dict_list[pert]]
                if not isinstance(comp_err_dict_list[pert], list):
                    # If type is not list then make it list
                    comp_err_dict_list[pert] = [comp_err_dict_list[pert]]
                # Append the value in list
                comp_dict_list[pert].append(ddG)
                comp_err_dict_list[pert].append(ddG_err)
            else:
                # As key is not in dict,
                # so, add key-value pair
                comp_dict_list[pert] = ddG
                comp_err_dict_list[pert] = ddG_err

    # now calculate all the avg and SEM for the network perturbations
    # put these into a dictionary
    comp_diff_dict = {}

    # TODO if none, don't write to csv file
    # write these to a csv file
    with open(f"{file_name}.csv", "w") as comp_pert_file:
        writer = csv.writer(comp_pert_file, delimiter=",")
        writer.writerow(["lig_1","lig_2","freenrg","error","engine"])
        for pert in perturbations:
            lig_0 = pert.split("~")[0]
            lig_1 = pert.split("~")[1]
            
            # check if the perturbations calculated are also those in the network file and if any are missing
            try:
                # find the values in the dictionary
                ddGs = comp_dict_list[pert]
                ddGs_error = comp_err_dict_list[pert]
                # calculate the average and the error
                comp_ddG = np.average(ddGs)
                # comp_ddG = np.average([ddG.value() for ddG in ddGs])
                if len(ddGs) == 1:
                    comp_err = ddGs_error.value()
                else:
                    comp_err = sem(ddGs)
                    # comp_err = sem([ddG.value() for ddG in ddGs])
        
            # if unable to calculate one of the perturbations, this is a None value.
            except:
                comp_ddG = None
                comp_err = None
            
            # TODO some way to incl units

            #update the dictionary for plotting later
            comp_diff_dict.update({pert:(comp_ddG, comp_err)})

            writer.writerow([lig_0, lig_1, comp_ddG, comp_err, engine])
    
    return comp_diff_dict

def freenrgworkflows_into_dict(experimental_DDGs, ligands, perturbations):
    # create a dictionary for the experimental values
    exper_val_dict = {}

    # convert the list of dicitonaries from freenrgworkflows into a single dictionary
    for lig_dict in experimental_DDGs:
        lig_name = list(lig_dict.keys())[0]
        exper = lig_dict[lig_name]
        exper_err = lig_dict["error"]
        exper_val_dict.update({lig_name:(exper, exper_err)})

    # add any ligands that are in the ligands file but dont have experimental values for
    for lig_name in ligands:
        if lig_name in exper_val_dict:
            pass
        else:
            exper_val_dict.update({lig_name:(None, None)})

    # now that we have our dictionary, 
    # we can also create a dictionary with all the experimental values for the perturbations
    exper_diff_dict = {}

    # calculate the experimental RBFEs
    # write these to a csv file
    with open("experimental_perturbations.csv", "w") as exp_pert_file:
        writer = csv.writer(exp_pert_file, delimiter=",")
        writer.writerow(["lig_1","lig_2","freenrg","error","engine"])

        for pert in perturbations:
            lig_0 = pert.split("~")[0]
            lig_1 = pert.split("~")[1]
            # exclude from calculating if one of the ligands is not available
            if exper_val_dict[lig_0][0] is None or exper_val_dict[lig_1][0] is None:
                exper_ddG =None
                exper_err = None
                exper_diff_dict.update({pert:(None, None)})
            # if experimental data is available, calculate experimental perturbation and propagate
            else:
                exper_ddG = exper_val_dict[lig_1][0] - exper_val_dict[lig_0][0]
                exper_err = math.sqrt(math.pow(exper_val_dict[lig_0][1], 2) + math.pow(exper_val_dict[lig_1][1], 2))
                exper_diff_dict.update({pert:(exper_ddG, exper_err)})

            writer.writerow([lig_0, lig_1, exper_ddG, exper_err, "experimental"])

    return exper_diff_dict, exper_val_dict



def calc_mae(values_dict=None, perts_ligs = None):
    # calc mae for a provided dictionary in the format wanted

    values_dict_list = []
    for key in values_dict.keys():
        values_dict_list.append(key)

    mae_pert_df = pd.DataFrame(columns=values_dict_list,index=values_dict_list)
    mae_pert_df_err = pd.DataFrame(columns=values_dict_list,index=values_dict_list)

    # iterate over all possible combinations
    for combo in itertools.product(values_dict_list,values_dict_list):
        eng1 = combo[0]
        eng2 = combo[1]

        eng1_vals = []
        eng2_vals = []

        # first create df of values
        if perts_ligs == "perts":
            for pert in values_dict[eng1]["pert_results"]:
                if pert in values_dict[eng2]["pert_results"]:
                    if values_dict[eng1]["pert_results"][pert][0] != None:
                        if values_dict[eng2]["pert_results"][pert][0] != None:
                            eng1_vals.append(values_dict[eng1]["pert_results"][pert][0])
                            eng2_vals.append(values_dict[eng2]["pert_results"][pert][0])
        elif perts_ligs == "ligs":
            for pert in values_dict[eng1]["val_results"]:
                if pert in values_dict[eng2]["val_results"]:
                    if values_dict[eng1]["val_results"][pert][0] != None:
                        if values_dict[eng2]["val_results"][pert][0] != None:
                            eng1_vals.append(values_dict[eng1]["val_results"][pert][0])
                            eng2_vals.append(values_dict[eng2]["val_results"][pert][0])
        else:
            raise ValueError("'perts_ligs' must be either 'perts' or 'ligs'!")


        if len(eng1_vals) == len(eng2_vals):
            mean_absolute_error = mae(eng1_vals,eng2_vals)  
            data_for_df = {"eng1":eng1_vals,"eng2":eng2_vals}
            data_df= pd.DataFrame(data_for_df)
        else:
            print("cant calc")

        boots = []
        n_boots = 10000

        for n in range(n_boots):
            sample_df = data_df.sample(n=len(eng1_vals), replace=True)
            mae_sample = (abs(sample_df['eng1'] - sample_df['eng2']).sum())/len(eng1_vals)
            boots.append(mae_sample)
            mae_err = (np.std(boots))

        mae_pert_df.loc[eng1,eng2]=mean_absolute_error
        mae_pert_df_err.loc[eng1,eng2]=mae_err

    mae_pert_df.to_csv(f"{res_folder}/mae_pert_{file_ext_out}.csv", sep=" ")
    mae_pert_df_err.to_csv(f"{res_folder}/mae_pert_err_{file_ext_out}.csv", sep=" ")

    return mae_pert_df, mae_pert_df_err