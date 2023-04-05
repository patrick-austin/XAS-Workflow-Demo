import csv
import sys

import json

from larch import Interpreter
from larch.fitting import guess, param, param_group
from larch.io import read_athena, extract_athenagroup
from larch.symboltable import Group
from larch.xafs import (
    autobk,
    feffit,
    feffit_report,
    FeffitDataSet,
    FeffPathGroup,
    pre_edge,
    TransformGroup,
    xftf,
)

import matplotlib
import matplotlib.pyplot as plt


def get_groups(athena_project):
    athena_groups = []
    group_keys = list(athena_project._athena_groups.keys())
    for group_key in group_keys:
        gr_0 = extract_athenagroup(athena_project._athena_groups[group_key])
        athena_groups.append(gr_0)
    return athena_groups


def read_csv_data(input_file, id_field="id"):
    csv_data = {}
    try:
        with open(input_file, encoding="utf8") as csvfile:
            reader = csv.DictReader(csvfile, skipinitialspace=True)
            for row in reader:
                csv_data[int(row[id_field])] = row
    except FileNotFoundError:
        print("The specified file does not exist")
    return csv_data


def calc_with_defaults(xafs_group: Group) -> Group:
    """Calculate pre_edge and background with default arguments"""
    pre_edge(xafs_group)
    autobk(xafs_group)
    xftf(xafs_group)
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
    plt.savefig("rmr.png", format="png")


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
    fig.savefig("chikr.png", format="png")


def read_gds(gds_file, session):
    gds_pars = read_csv_data(gds_file)
    dgs_group = dict_to_gds(gds_pars, session)
    return dgs_group


def read_selected_paths_list(file_name, session):
    sp_dict = read_csv_data(file_name)
    sp_list = []
    for path_id in sp_dict:
        new_path = FeffPathGroup(
            filename="feff/" + sp_dict[path_id]["filename"],
            label=sp_dict[path_id]["label"],
            s02=sp_dict[path_id]["s02"],
            e0=sp_dict[path_id]["e0"],
            sigma2=sp_dict[path_id]["sigma2"],
            deltar=sp_dict[path_id]["deltar"],
            _larch=session,
        )
        sp_list.append(new_path)
    return sp_list


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
    gds_file: str,
    sp_file: str,
    fit_vars: dict,
    plot_graph: bool,
):
    session = Interpreter()
    athena_project = read_athena(prj_file)
    athena_groups = get_groups(athena_project=athena_project)
    data_group = calc_with_defaults(athena_groups[0])

    gds = read_gds(gds_file, session)
    selected_paths = read_selected_paths_list(sp_file, session)
    trans, dset, out = run_fit(data_group, gds, selected_paths, fit_vars, session)

    fit_report = feffit_report(out, _larch=session)
    with open("fit_report.txt", "w") as fit_report_file:
        fit_report_file.write(fit_report)

    if plot_graph:
        plot_rmr(dset, fit_vars["rmin"], fit_vars["rmax"])
        plot_chikr(
            dset, fit_vars["rmin"], fit_vars["rmax"], fit_vars["kmin"], fit_vars["kmax"]
        )


if __name__ == "__main__":
    # larch imports set this to an interactive backend, so need to change it
    matplotlib.use("Agg")

    prj_file = sys.argv[1]
    gds_file = sys.argv[2]
    sp_file = sys.argv[3]
    input_values = json.load(open(sys.argv[4], "r", encoding="utf-8"))
    fit_vars = input_values["fit_vars"]
    plot_graph = input_values["plot_graph"]
    main(prj_file, gds_file, sp_file, fit_vars, plot_graph)
