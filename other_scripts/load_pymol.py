# from finlay

from pymol import cmd
import subprocess


def set_colours() -> None:
    """Set up preffered colours and display settings"""
    cmd.set_color("white", [1.0, 1.0, 1.0])
    # cmd.bg_color("white")
    cmd.set("antialias", 1)
    cmd.set("orthoscopic", 1)
    cmd.set("gamma", 1.15)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_fancy_sheets", 1)
    cmd.set_color("wred", [0.788, 0.000, 0.140])
    cmd.set_color("wblue", [0.31, 0.506, 0.686])
    cmd.set_color("wgold", [0.855, 0.647, 0.125])
    cmd.set_color("wgreen", [0.134, 0.545, 0.134])
    cmd.set_color("wgray", [0.800, 0.800, 0.800])
    cmd.set_color("wrose", [0.65, 0.47, 0.55])
    cmd.set_color("wpurple", [0.37, 0.31, 0.62])
    cmd.set_color("mpurple", [0.75, 0.57, 0.80])
    cmd.set_color("mpgrey", [0.73, 0.68, 0.82])
    cmd.set("ray_shadows", 0)
    cmd.set("ray_trace_fog", 1)
    cmd.set("ray_trace_mode", 1)
    cmd.set("shininess", 1000)


def load_complex(
    top_file: str = "somd.parm7",
    coord_file: str = "somd.rst7",
    traj_file: str = "traj000000001.dcd",
) -> None:
    """Load the complex and a trajectory"""
    # Load protein
    cmd.load(top_file, "sys")
    cmd.load(coord_file, "sys")

    # Ensure secondary structure is shown
    cmd.dss()

    # Hide boring stuff
    cmd.hide("everything")
    cmd.show("cartoon", "polymer")
    cmd.show("licorice", "organic")

    # Show binding site residues
    cmd.select("lig", "(sys and organic)")
    cmd.select("bs_res", "(byres polymer within 5 of lig)")
    cmd.show("licorice", "bs_res")
    cmd.select("bs_waters", "(solvent within 5 of lig)")
    cmd.show("licorice", "bs_waters")

    # Center camera on ligand
    cmd.center("lig")

    # Load trajectory and align to the protein
    if traj_file:
        cmd.load_traj(traj_file)
        cmd.intra_fit("polymer")
        cmd.smooth()

    # Uncomment below to update selection of waters each frame. This is very slow.
    # Define the movie range based on the number of frames in the trajectory
    # n_frames = cmd.count_states("sys")  # Get the number of frames in the trajectory
    # movie_range = f"1-{n_frames}"  # Set the movie range

    # Set the movie range
    # cmd.mset(movie_range)

    # Loop through the trajectory and update the selection
    # for f_no in range(1, n_frames + 1):
    # cmd.mdo(
    # f_no,
    # "hide licorice, bs_waters; select bs_waters, (solvent within 5 of lig), state="
    # + str(f_no)
    # + "; show licorice, bs_waters",
    # )

    # Play the movie
    # cmd.mplay()


def unwrap_gmx_traj(
    tpr_file: str = "gromacs.gro",
    traj_file: str = "gromacs.xtc",
) -> str:
    """
    Unwrap a GROMACS trajectory using trjconv. Returns the name of the
    unwrapped trajectory file.
    """
    cmd = f"gmx trjconv -f {traj_file} -s {tpr_file} -o unwrapped.xtc -pbc mol"
    subprocess.run(cmd, shell=True, check=True)
    return "unwrapped.xtc"


def main():
    import pymol
    import argparse

    # Parse command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--folder", help="Folder for the files", default=None, type=str)

    parser.add_argument("--top", help="Topology file", default=None, type=str)
    parser.add_argument(
        "--coord", help="Coordinate file", default=None, type=str
    )
    parser.add_argument(
        "--traj",
        help="Trajectory file. Specifying --traj '' will cause the trajectory to be ignored.",
        default=None,
        type=str,
    )
    parser.add_argument(
        "--unwrap",
        help="Unwrap trajectory. This is only required for GROMACS trajectories.",
        action="store_true",
    )
    # Also read in tpr file if --unwrap is specified
    parser.add_argument(
        "--tpr",
        help="GROMACS tpr file. Only required if --unwrap is specified.",
        default=None,
        type=str,
    )

    args = parser.parse_args()

    if args.folder:
        if not args.top:
            top_file = f"{args.folder}/somd.prm7"
        else:
            top_file = args.top
        if not args.coord:
            coord_file = f"{args.folder}/somd.rst7"
        else:
            coord_file = args.coord
        if not args.traj:
            traj_file = f"{args.folder}/traj000000001.dcd"
        else:
            coord_file = args.traj
        if not args.tpr:
            tpr_file = f"{args.folder}/gromacs.tpr"
        else:
            coord_file = args.tpr
    else:
        if args.top:
            top_file = args.top
        else:
            top_file = "somd.prm7"
        if args.coord:
            coord_file = args.coord
        else:
            coord_file = "somd.rst7"
        if args.traj:
            traj_file = args.traj
        else:
            traj_file = "traj000000001.dcd"
        if args.tpr:
            tpr_file = args.tpr
        else:
            tpr_file = "gromacs.tpr"

    # Soft link to the prm7 file if supplied
    if top_file[-4:] == "prm7":
        if args.folder:
            subprocess.run(
                f"ln -sf {top_file} {args.folder}/somd.parm7", shell=True, check=True)
            top_file = f"{args.folder}/somd.parm7"
        else:
            subprocess.run(
                f"ln -sf {top_file} somd.parm7", shell=True, check=True)
            top_file = "somd.parm7"

    # Unwrap trajectory if required
    if args.unwrap:
        if not traj_file or not tpr_file:
            raise ValueError(
                "Must specify both --traj and --tpr if --unwrap is specified."
            )
        traj_file = unwrap_gmx_traj(tpr_file, traj_file)

    # Run
    pymol.finish_launching()
    set_colours()
    load_complex(top_file, coord_file, traj_file)


if __name__ == "__main__":
    main()
