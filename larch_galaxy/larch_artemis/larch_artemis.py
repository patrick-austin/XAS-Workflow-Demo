import os
from pathlib import Path
import sys

import json
from typing import List

from larch import Interpreter
from larch.io import create_athena, read_ascii
from larch.symboltable import Group
from larch.utils import group2dict, dict2group
from larch.xafs import pre_edge

import matplotlib.pyplot as plt

from lib import atoms_feff, handle_csv, manage_athena, manage_fit


def calc_with_defaults(xafs_group: Group) -> Group:
    """Calculate pre_edge and background with default arguments"""
    # TODO not sure what the implications of this change are?
    # Didn't have background in our data so had to remove it from the pre_edge
    pre_edge(xafs_group)
    manage_athena.autobk(xafs_group)
    manage_athena.xftf(xafs_group)
    return xafs_group


def plot_normalised_all(plot_path: str, xafs_group: List[Group]):
    for i, group in enumerate(xafs_group):
        plt.subplot(len(xafs_group), 1, i + 1)
        plt.plot(group.energy, group.mu, label=group.filename)
        plt.grid(color="r", linestyle=":", linewidth=1)
        plt.xlabel("Energy (eV)")
        plt.ylabel("x$\mu$(E)")  # noqa: W605
        plt.title("pre-edge and post_edge fitting to $\mu$")  # noqa: W605
        plt.legend()
    plt.savefig(plot_path, format="png")


def run_feff(input_files: List[str]) -> List[str]:
    """Copies each input file and runs FEFF6l calculation"""
    # TODO check implications of changes
    # Removed the calls to perl which convert to FEFF inp format
    # TODO can we remove the implicit reliance on structure names in filepaths?
    feff_dir_list = []
    for inp_file in input_files:
        crystal_f = Path(inp_file)
        # use the name of the input file to define the
        # names of the feff directory and inp file
        feff_dir = crystal_f.name[:-4] + "_feff"
        feff_inp = crystal_f.name[:-4] + "_feff.inp"
        path = Path(feff_dir, feff_inp)
        atoms_ok = atoms_feff.copy_to_feff_dir(crystal_f, path)
        if atoms_ok:
            # run feff to generate the scattering paths
            atoms_feff.feff6l(folder=feff_dir, feffinp=feff_inp)
            feff_dir_list.append(feff_dir)
    return feff_dir_list


def main(
    prj_file: str,
    structure_files: List[str],
    gds_file: str,
    sp_file: str,
    fit_vars: dict,
    plot_graph: bool,
):
    session = Interpreter()
    # TODO don't know how to produce this, but the repo implies we should
    # fpj_path = "out.fpj"
    athena_project = manage_athena.read_project(project_file=prj_file)
    athena_groups = manage_athena.get_groups(athena_project=athena_project)
    data_group = calc_with_defaults(athena_groups[0])  # TODO hardcoded 0?

    # Generate the FEFF dirs for the path selection
    run_feff(structure_files)

    # TODO In principle might want to edit this rather than loading from file
    gds = manage_fit.read_gds(gds_file, session)

    # TODO In principle might want to edit this rather than loading from file
    # TODO hardcoded 0
    # path_sheet = manage_fit.show_feff_paths(structure_files[0])
    # sp_sheet = manage_fit.show_selected_paths(path_sheet)
    # selected_paths = manage_fit.build_selected_paths_list(sp_sheet, session)
    # manage_fit.save_selected_paths_list(sp_sheet, spl_file)
    selected_paths = manage_fit.read_selected_paths_list(sp_file, session)

    # TODO removed code for showing the spreadsheet - doesn't work without an ipynb
    # TODO consider replacing with a pandas table or similar?

    # TODO not doing anything with trans?
    trans, dset, out = manage_fit.run_fit(
        data_group, gds, selected_paths, fit_vars, session
    )

    fit_report = manage_fit.get_fit_report(out, session)
    with open("fit_report.txt", "w") as fit_report_file:
        fit_report_file.write(fit_report)

    if plot_graph:
        plot_normalised_all(plot_path="norm.png", xafs_group=athena_groups)
        rmr_p = manage_fit.plot_rmr(dset, fit_vars["rmin"], fit_vars["rmax"])
        rmr_p.savefig("rmr.png", format="png")
        chikr_p = manage_fit.plot_chikr(
            dset, fit_vars["rmin"], fit_vars["rmax"], fit_vars["kmin"], fit_vars["kmax"]
        )
        chikr_p.savefig("chikr.png", format="png")


if __name__ == "__main__":
    prj_file = sys.argv[1]
    structure_files = sys.argv[2].split(",")[:-1]  # Account for trailing ","
    gds_file = sys.argv[3]
    sp_file = sys.argv[4]
    inputs = json.load(open(sys.argv[5], "r", encoding="utf-8"))
    fit_vars = inputs["fit_vars"]
    plot_graph = inputs["plot_graph"]
    main(prj_file, structure_files, gds_file, sp_file, fit_vars, plot_graph)
