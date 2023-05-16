import csv
import json
import re
import sys


SP_DATA = [
    f"{'id':>4s}, {'filename':>12s}, {'label':>24s}, {'s02':>3s}, {'e0':>4s}, {'sigma2':>24s}, {'deltar':>10s}\n"
]
GDS_DATA = [f"{'id':>4s}, {'name':>24s}, {'value':>5s}, {'expr':>4s}, {'vary':>4s}\n"]
SP_ROW_ID = 1
GDS_ROW_ID = 1
AMP = "amp"
ENOT = "enot"
ALPHA = "alpha"
ALPHA_REFF = "alpha*reff"


def write_selected_path(
    row: "list[str]",
    s02: str = AMP,
    e0: str = ENOT,
    sigma2: str = "",
    deltar: str = ALPHA_REFF,
):
    global SP_ROW_ID
    filename = row[0].strip()
    label = row[-2].strip()
    if not sigma2:
        sigma2 = "s" + label.replace(".", "").lower()
        write_gds(sigma2)
    SP_DATA.append(
        f"{SP_ROW_ID:>4d}, {filename:>12s}, {label:>24s}, {s02:>3s}, {e0:>4s}, {sigma2:>24s}, {deltar:>10s}\n"
    )
    SP_ROW_ID += 1


def write_gds(name: str, value: float = 0.003, expr: str = None, vary: bool = True):
    global GDS_ROW_ID
    if not expr:
        expr = "    "
    GDS_DATA.append(
        f"{GDS_ROW_ID:4d}, {name:>24s}, {str(value):>5s}, {expr:>4s}, {str(vary):>4s}\n"
    )
    GDS_ROW_ID += 1


def main(
    paths_file: str,
    select_all: bool,
    path_values: list,
    amp_values: dict,
    enot_values: dict,
    alpha_values: dict,
    gds_values: list,
):
    path_values_ids = [path_value["id"] for path_value in path_values]
    gds_names = [gds_value["name"] for gds_value in gds_values]

    write_gds(AMP, **amp_values)
    write_gds(ENOT, **enot_values)
    write_gds(ALPHA, **alpha_values)
    for gds_value in gds_values:
        write_gds(**gds_value)

    with open(paths_file) as file:
        reader = csv.reader(file)
        for row in reader:
            id_match = re.search(r"\d+", row[0])
            if id_match:
                path_id = int(id_match.group())
                if path_id in path_values_ids:
                    path_value = path_values[path_values_ids.index(path_id)]
                    s02 = path_value["s02"]
                    e0 = path_value["e0"]
                    sigma2 = path_value["sigma2"]
                    if sigma2 and sigma2 not in gds_names:
                        write_gds(sigma2)
                    deltar = path_value["deltar"]
                    write_selected_path(row, s02, e0, sigma2, deltar)
                elif select_all or int(row[-1]):
                    write_selected_path(row)

    with open("sp.csv", "w") as out:
        out.writelines(SP_DATA)

    with open("gds.csv", "w") as out:
        out.writelines(GDS_DATA)


if __name__ == "__main__":
    paths_file = sys.argv[1]
    input_values = json.load(open(sys.argv[2], "r", encoding="utf-8"))
    main(
        paths_file,
        input_values["select_all"],
        input_values["paths"],
        input_values["amp"],
        input_values["enot"],
        input_values["alpha"],
        input_values["gds"],
    )
