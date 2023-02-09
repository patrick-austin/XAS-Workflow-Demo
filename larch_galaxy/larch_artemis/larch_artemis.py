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


def rename_cols(xafs_group: Group):
    mu_e = xafs_group.xmu
    # get a dictionary from te group
    xafs_dict = group2dict(xafs_group)
    # add mu and energy to the dictionary
    xafs_dict["mu"] = mu_e
    xafs_group = dict2group(xafs_dict)
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


def main(prj_file: str, structure_file: str, gds_file: str, sp_file: str,
        fit_vars: dict, plot_graph: bool):
    plot_path = "out.png"
    athena_project = manage_athena.read_project(project_file=prj_file)
    athena_groups = manage_athena.get_groups(athena_project=athena_project)

    group_names = [group.label for group in athena_groups]
    print(f"group names: {group_names}")

    data_group = manage_athena.calc_with_defaults(athena_groups[0])

    if plot_graph:
        plot_normalised_all(plot_path=plot_path, xafs_group=athena_groups)

    # TODO type? if structure file is multiple is this a string or a list?
    feff_dirs = atoms_feff.run_feff(structure_file)

    # TODO
    # At this point we're supposed to use the fit manager but can't because it uses the larch_plugins import
    return
    path_sheet = manage_fit.show_feff_paths(structure_file[0])
    display(path_sheet)  # This isn't a function, even in the namespace of the notebook - check some of the other material?



if __name__ == "__main__":
    prj_file = sys.argv[1]
    structure_file = sys.argv[2]
    gds_file = sys.argv[3]
    sp_file = sys.argv[4]
    inputs = json.load(open(sys.argv[5], "r", encoding="utf-8"))
    fit_vars = inputs["fit_vars"]
    plot_graph = inputs["true"]
    main(prj_file, structure_file, gds_file, sp_file, fit_vars, plot_graph)
