#!/usr/bin/python3

import sys
from argparse import ArgumentParser
import glob

print("adding code to the pythonpath...")
code = "/home/anna/Documents/code/python"
if code not in sys.path:
    sys.path.insert(1, code)
import pipeline

from pipeline.prep import *
from pipeline.utils import *
from pipeline.analysis import *

import subprocess


def search_string_with_grep(file_path, search_string):
    try:
        output = subprocess.check_output(["grep", search_string, file_path])
        return str(output.strip())
    except subprocess.CalledProcessError:
        return False


def check_arguments(args):
    # pass the checks to the other check functions

    if args.main_folder:
        main_folder = args.main_folder
    else:
        main_folder = str(
            input(
                "what is the main folder where the pipeline was run (should include slurm_logs folder in it)?: "
            )
        ).strip()

    if args.protocol_file:
        protocol_file = args.protocol_file
    else:
        protocol_file = str(
            input(
                "what is the path to the protocol file? (default if left blank: execution_model/protocol.dat): "
            ).strip()
        )
        if not protocol_file:
            protocol_file = f"{main_folder}/execution_model/protocol.dat"

    # if args.network_file:
    #     network_file = args.network_file
    # else:
    #     try:
    #         network_file = validate.file_path(
    #             f"{main_folder}/execution_model/network_combined.dat"
    #         )
    #     except:
    #         network_file = str(input("what is the path to the network file?: ").strip())

    network_file = None

    return main_folder, protocol_file, network_file


def main():
    # accept all options as arguments
    parser = ArgumentParser(description="featurise the perturbations")
    parser.add_argument(
        "-mf",
        "--main_folder",
        type=str,
        default=None,
        help="main folder for the run",
    )

    parser.add_argument(
        "-p", "--protocol_file", type=str, default=None, help="path to protocol file"
    )

    parser.add_argument(
        "-n", "--network_file", type=str, default=None, help="path to network file"
    )
    args = parser.parse_args()

    # check arguments
    print("checking the provided command line arguments...")
    main_folder, protocol_file, network_file = check_arguments(args)
    prot = pipeline_protocol(protocol_file, auto_validate=True)

    # check all the perturbations in the network
    perts, ligs = get_info_network(net_file=network_file)

    slurm_log_files = glob.glob(f"{main_folder}/slurm_logs/prod*")

    minutes_list = []

    for eng in prot.engines():
        for file in slurm_log_files:
            has_engine = search_string_with_grep(file, eng)

            if has_engine:
                output = search_string_with_grep(file, "runtime")
                if output:
                    minutes = int(output.split()[-2])
                    minutes_list.append(minutes)
                else:
                    pass

    print(f"average per production log is: {np.mean(minutes_list)} ")
    print(
        f"divided by the number of repeats from the protocol: {np.mean(minutes_list) / prot.repeats()} "
    )


if __name__ == "__main__":
    main()
