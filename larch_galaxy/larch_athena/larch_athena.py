import sys

from larch.io import create_athena, read_ascii
from larch.symboltable import Group
from larch.utils import group2dict, dict2group
from larch.xafs import pre_edge

import matplotlib.pyplot as plt


def rename_cols(xafs_group: Group):
    mu_e = xafs_group.xmu
    # get a dictionary from te group
    xafs_dict = group2dict(xafs_group)
    # add mu and energy to the dictionary
    xafs_dict["mu"] = mu_e
    print(f"XAFS data keys: {xafs_dict.keys()}")
    xafs_group = dict2group(xafs_dict)
    return xafs_group


def plot_normalised(plot_path: str, xafs_group: Group):
    plt.plot(xafs_group.energy, xafs_group.pre_edge, "g", label="pre-edge")
    plt.plot(xafs_group.energy, xafs_group.post_edge, "r", label="post-edge")
    plt.plot(xafs_group.energy, xafs_group.mu, "b", label=xafs_group.filename)
    plt.grid(color="r", linestyle=":", linewidth=1)
    plt.xlabel("Energy (eV)")
    plt.ylabel("x$\mu$(E)")  # noqa: W605
    plt.title("pre-edge and post_edge fitting to $\mu$")  # noqa: W605
    plt.legend()
    plt.savefig(plot_path, format="png")


def main(dat_file: str, plot_graph: bool):
    prj_path = "out.prj"
    plot_path = "out.png"

    xas_data: Group = read_ascii(dat_file)
    xas_data = rename_cols(xas_data)
    xas_data.filename = prj_path

    pre_edge(energy=xas_data.energy, mu=xas_data.mu, group=xas_data)

    if plot_graph:
        plot_normalised(plot_path, xas_data)

    xas_project = create_athena(prj_path)
    xas_project.add_group(xas_data)
    xas_project.save()


if __name__ == "__main__":
    dat_file = sys.argv[1]
    plot_graph = sys.argv[2] == "true"
    main(dat_file, plot_graph)
