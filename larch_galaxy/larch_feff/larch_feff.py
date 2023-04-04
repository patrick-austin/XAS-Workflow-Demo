import re
import shutil
import sys
from pathlib import Path

from larch.xafs.feffrunner import feff6l


def get_path_labels(paths_file):
    is_meta = True
    count = 0
    a_path = {}
    all_paths = {}
    with open(paths_file) as datfile:
        dat_lines = datfile.readlines()
        for a_line in dat_lines:
            count += 1
            if re.match("-{15}", a_line.strip()) is not None:
                is_meta = False
            elif not is_meta:
                if re.match("\s*\d*\s{4}\d*\s{3}", a_line) is not None:
                    if a_path:
                        all_paths[a_path["index"]] = a_path
                    line_data = a_line.split()
                    a_path = {
                        "index": line_data[0],
                        "nleg": line_data[1],
                        "degeneracy": line_data[2],
                    }
                elif (
                    re.match("\s{6}x\s{11}y\s{5}", a_line) is None
                ):  # ignore the intermediate headings
                    line_data = a_line.split()
                    if "label" not in a_path:
                        a_path["label"] = line_data[4].replace("'", "")
                    else:
                        a_path["label"] += "." + line_data[4].replace("'", "")
    if a_path and "index" in a_path:
        all_paths[a_path["index"]] = a_path
    return all_paths


def main(structure_file: str):
    crystal_f = Path(structure_file)
    feff_dir = "feff"
    feff_inp = "feff.inp"
    path = Path(feff_dir, feff_inp)
    try:
        path.parent.mkdir(parents=True, exist_ok=True) 
        print(f"Copying {crystal_f.name} to {path}")
        shutil.copy(crystal_f, path)
        feff6l(folder=feff_dir, feffinp=feff_inp)
    except Exception as exception:
        raise ValueError("Unable to run feff") from exception

    feff_files = "files.dat"
    input_file = Path(feff_dir, feff_files)
    # the .dat data is stored in fixed width strings
    comma_positions = [13, 21, 32, 41, 48, 61]
    paths_data = []
    # get the list of paths info to assign labels to paths
    paths_info = get_path_labels(Path(feff_dir, "paths.dat"))
    print("Reading from: " + str(input_file))
    with open(input_file) as datfile:
        # Read until we find the line at the end of the header
        line = datfile.readline()
        while not re.match("-+", line.strip()):
            line = datfile.readline()

        # Build headers
        line = datfile.readline()
        header = ""
        start = 0
        for end in comma_positions:
            header += line[start:end] + ","
            start = end
        header += f" {'label':12s}, {'select':6d}"
        paths_data.append(header)

        # Read data
        line = datfile.readline()
        while line:
            data = ""
            start = 0
            for end in comma_positions[:-1]:
                data += line[start:end] + ","
                start = end

            # Last column needs padding to align
            data += line[start:] + "     ,"

            # Add label and select column
            path_id = data[5:9]
            label = paths_info[path_id]["label"] + f".{int(path_id):d}"
            data += f" {label:12s}, {0:6d}"
            paths_data.append(data)
            line = datfile.readline()

    with open("out.csv", "w") as out:
        out.writelines(paths_data)


if __name__ == "__main__":
    structure_file = sys.argv[1]
    main(structure_file)
