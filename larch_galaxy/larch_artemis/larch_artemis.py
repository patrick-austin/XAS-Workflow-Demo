from pathlib import Path
import sys

import json
from typing import List

from larch import Interpreter
from larch.fitting import guess, param, param_group
from larch.symboltable import Group
from larch.xafs import (
    feffit,
    feffit_report,
    FeffitDataSet,
    FeffPathGroup,
    pre_edge,
    TransformGroup,
)

import matplotlib
import matplotlib.pyplot as plt

from lib import atoms_feff, handle_csv, manage_athena


def calc_with_defaults(xafs_group: Group) -> Group:
    """Calculate pre_edge and background with default arguments"""
    # TODO not sure what the implications of this change are?
    # Didn't have background in our data so had to remove it from the pre_edge
    pre_edge(xafs_group)
    manage_athena.autobk(xafs_group)
    manage_athena.xftf(xafs_group)
    return xafs_group


def dict_to_gds(data_dict, session):
    dgs_group = param_group(_larch=session)
    for par_idx in data_dict:
        # gds file structure:
        gds_name = data_dict[par_idx]["name"]
        gds_val = 0.0
        gds_expr = ""
        try:
            gds_val = float(data_dict[par_idx]["value"])
        except ValueError:
            gds_val = 0.00
        gds_expr = data_dict[par_idx]["expr"]
        gds_vary = (
            True
            if str(data_dict[par_idx]["vary"]).strip().capitalize() == "True"
            else False
        )
        one_par = None
        if gds_vary:
            # equivalent to a guess parameter in Demeter
            one_par = guess(name=gds_name, value=gds_val, vary=gds_vary, expr=gds_expr)
        else:
            # equivalent to a defined parameter in Demeter
            one_par = param(name=gds_name, value=gds_val, vary=gds_vary, expr=gds_expr)
        if one_par is not None:
            dgs_group.__setattr__(gds_name, one_par)
    return dgs_group


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


def plot_rmr(data_set, rmin, rmax):
    plt.figure()
    plt.plot(data_set.data.r, data_set.data.chir_mag, color="b")
    plt.plot(data_set.data.r, data_set.data.chir_re, color="b", label="expt.")
    plt.plot(data_set.model.r, data_set.model.chir_mag, color="r")
    plt.plot(data_set.model.r, data_set.model.chir_re, color="r", label="fit")
    plt.ylabel("Magnitude of Fourier Transform of $k^2 \cdot \chi$/$\mathrm{\AA}^{-3}$")
    plt.xlabel("Radial distance/$\mathrm{\AA}$")
    plt.xlim(0, 5)

    plt.fill([rmin, rmin, rmax, rmax], [-rmax, rmax, rmax, -rmax], color="g", alpha=0.1)
    plt.text(rmax - 0.65, -rmax + 0.5, "fit range")
    plt.legend()
    return plt


def plot_chikr(data_set, rmin, rmax, kmin, kmax):
    fig = plt.figure(figsize=(16, 4))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.plot(
        data_set.data.k,
        data_set.data.chi * data_set.data.k**2,
        color="b",
        label="expt.",
    )
    ax1.plot(
        data_set.model.k,
        data_set.model.chi * data_set.data.k**2,
        color="r",
        label="fit",
    )
    ax1.set_xlim(0, 15)
    ax1.set_xlabel("$k (\mathrm{\AA})^{-1}$")
    ax1.set_ylabel("$k^2$ $\chi (k)(\mathrm{\AA})^{-2}$")

    ax1.fill([kmin, kmin, kmax, kmax], [-rmax, rmax, rmax, -rmax], color="g", alpha=0.1)
    ax1.text(kmax - 1.65, -rmax + 0.5, "fit range")
    ax1.legend()

    ax2.plot(data_set.data.r, data_set.data.chir_mag, color="b", label="expt.")
    ax2.plot(data_set.model.r, data_set.model.chir_mag, color="r", label="fit")
    ax2.set_xlim(0, 5)
    ax2.set_xlabel("$R(\mathrm{\AA})$")
    ax2.set_ylabel("$|\chi(R)|(\mathrm{\AA}^{-3})$")
    ax2.legend(loc="upper right")

    ax2.fill([rmin, rmin, rmax, rmax], [-rmax, rmax, rmax, -rmax], color="g", alpha=0.1)
    ax2.text(rmax - 0.65, -rmax + 0.5, "fit range")
    return plt


def read_gds(gds_file, session):
    gds_pars, _ = handle_csv.read_csv_data(gds_file)
    dgs_group = dict_to_gds(gds_pars, session)
    return dgs_group


def read_selected_paths_list(file_name, session):
    sp_dict, _ = handle_csv.read_csv_data(file_name)
    sp_list = []
    for path_id in sp_dict:
        new_path = FeffPathGroup(
            filename=sp_dict[path_id]["filename"],
            label=sp_dict[path_id]["label"],
            s02=sp_dict[path_id]["s02"],
            e0=sp_dict[path_id]["e0"],
            sigma2=sp_dict[path_id]["sigma2"],
            deltar=sp_dict[path_id]["deltar"],
            _larch=session,
        )
        sp_list.append(new_path)
    return sp_list


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


def run_fit(data_group, gds, selected_paths, fv, session):
    # create the transform grup (prepare the fit space).
    trans = TransformGroup(
        fitspace=fv["fitspace"],
        kmin=fv["kmin"],
        kmax=fv["kmax"],
        kw=fv["kw"],
        dk=fv["dk"],
        window=fv["window"],
        rmin=fv["rmin"],
        rmax=fv["rmax"],
        _larch=session,
    )

    dset = FeffitDataSet(
        data=data_group, pathlist=selected_paths, transform=trans, _larch=session
    )

    out = feffit(gds, dset, _larch=session)
    return trans, dset, out


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
    gds = read_gds(gds_file, session)

    # TODO In principle might want to edit this rather than loading from file
    # TODO hardcoded 0
    # path_sheet = manage_fit.show_feff_paths(structure_files[0])
    # sp_sheet = manage_fit.show_selected_paths(path_sheet)
    # selected_paths = manage_fit.build_selected_paths_list(sp_sheet, session)
    # manage_fit.save_selected_paths_list(sp_sheet, spl_file)
    selected_paths = read_selected_paths_list(sp_file, session)

    # TODO removed code for showing the spreadsheet - doesn't work without an ipynb
    # TODO consider replacing with a pandas table or similar?

    # TODO not doing anything with trans?
    trans, dset, out = run_fit(data_group, gds, selected_paths, fit_vars, session)

    fit_report = feffit_report(out, _larch=session)
    with open("fit_report.txt", "w") as fit_report_file:
        fit_report_file.write(fit_report)

    if plot_graph:
        plot_normalised_all(plot_path="norm.png", xafs_group=athena_groups)
        rmr_p = plot_rmr(dset, fit_vars["rmin"], fit_vars["rmax"])
        rmr_p.savefig("rmr.png", format="png")
        chikr_p = plot_chikr(
            dset, fit_vars["rmin"], fit_vars["rmax"], fit_vars["kmin"], fit_vars["kmax"]
        )
        chikr_p.savefig("chikr.png", format="png")


if __name__ == "__main__":
    # larch imports set this to an interactive backend, so need to change it
    matplotlib.use("Agg")

    prj_file = sys.argv[1]
    structure_files = sys.argv[2].split(",")[:-1]  # Account for trailing ","
    gds_file = sys.argv[3]
    sp_file = sys.argv[4]
    inputs = json.load(open(sys.argv[5], "r", encoding="utf-8"))
    fit_vars = inputs["fit_vars"]
    plot_graph = inputs["plot_graph"]
    main(prj_file, structure_files, gds_file, sp_file, fit_vars, plot_graph)
