import gc
import os
import sys
from zipfile import ZipFile

from larch.io import create_athena, read_ascii, set_array_labels
from larch.symboltable import Group
from larch.utils import group2dict, dict2group
from larch.xafs import pre_edge

import matplotlib
import matplotlib.pyplot as plt


def rename_cols(xafs_group: Group):
    labels = xafs_group.array_labels
    new_labels = []
    for label in labels:
        if label == "col1":
            new_labels.append("energy")
        elif label == "col2" or label == "xmu":
            new_labels.append("mu")
        else:
            new_labels.append(label)

    if new_labels != labels:
        return set_array_labels(xafs_group, new_labels)
    else:
        return xafs_group


def plot_edge_fits(plot_path: str, xafs_group: Group):
    plt.figure()
    plt.plot(xafs_group.energy, xafs_group.pre_edge, "g", label="pre-edge")
    plt.plot(xafs_group.energy, xafs_group.post_edge, "r", label="post-edge")
    plt.plot(xafs_group.energy, xafs_group.mu, "b", label=xafs_group.filename)
    plt.grid(color="r", linestyle=":", linewidth=1)
    plt.xlabel("Energy (eV)")
    plt.ylabel("x$\mu$(E)")  # noqa: W605
    plt.title("pre-edge and post_edge fitting to $\mu$")  # noqa: W605
    plt.legend()
    plt.savefig(plot_path, format="png")
    plt.close("all")


def plot_flattened(plot_path: str, xafs_group: Group):
    plt.figure()
    plt.plot(xafs_group.energy, xafs_group.flat, label=xafs_group.filename)
    plt.grid(color="r", linestyle=":", linewidth=1)
    plt.xlabel("Energy (eV)")
    plt.ylabel("normalised x$\mu$(E)")  # noqa: W605
    plt.legend()
    plt.savefig(plot_path, format="png")
    plt.close("all")


def main(
    dat_file: str,
    plot_graph: bool,
    prj_path: str = "out.prj",
    edge_plot_path: str = "edge.png",
    flat_plot_path: str = "flat.png",
):
    xas_data: Group = read_ascii(dat_file)
    xas_data = rename_cols(xas_data)
    xas_data.filename = prj_path

    pre_edge(energy=xas_data.energy, mu=xas_data.mu, group=xas_data)

    if plot_graph:
        plot_edge_fits(edge_plot_path, xas_data)
        plot_flattened(flat_plot_path, xas_data)

    xas_project = create_athena(prj_path)
    xas_project.add_group(xas_data)
    xas_project.save()

    # Ensure that we do not run out of memory when running on large zips
    gc.collect()


if __name__ == "__main__":
    # larch imports set this to an interactive backend, so need to change it
    matplotlib.use("Agg")

    dat_file = sys.argv[1]
    plot_graph = sys.argv[2] == "true"
    extension = sys.argv[3]
    if extension == "zip":
        i = 0
        ZipFile(dat_file).extractall("dat_files")
        for dirpath, _, filenames in os.walk("dat_files"):
            for filename in filenames:
                main(
                    os.path.join(dirpath, filename),
                    plot_graph,
                    f"out/{i}.binary",
                    f"edge/{i}.png",
                    f"flat/{i}.png",
                )
                i += 1
    else:
        main(dat_file, plot_graph)
