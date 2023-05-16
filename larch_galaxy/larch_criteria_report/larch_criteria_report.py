import csv
import json
import sys
from typing import Iterable

import matplotlib.pyplot as plt


def plot(variable: str, column: Iterable[float]):
    variable_stripped = variable.strip()
    path = f"plots/{variable_stripped}.png"
    plt.figure(figsize=(8, 4))
    plt.plot(column)
    plt.xticks(range(len(column)))
    plt.xlabel("Dataset number")
    plt.ylabel(variable_stripped)
    plt.savefig(path, format="png")


def load(filepath: str) -> "list[list[str|float]]":
    with open(filepath) as f:
        reader = csv.reader(f)
        header = next(reader)
        columns = [[h] for h in header]
        for row in reader:
            for i, value in enumerate(row):
                columns[i].append(float(value))

    return columns


def parse_reports(input_data: str) -> "dict[str, list[float]]":
    # Need to scrape variables from individual files
    report_criteria = input_values["format"]["report_criteria"]
    data = {c["variable"]: [] for c in report_criteria}
    headers = list(data.keys())
    with open("criteria_report.csv", "w") as f_out:
        writer = csv.writer(f_out)
        writer.writerow([f"{h:>12s}" for h in headers])
        for input_file in input_data.split(","):
            row = parse_row(data, headers, input_file)
            writer.writerow(row)

    return data


def parse_row(
    data: "dict[str, list[float]]", headers: "list[str]", input_file: str
) -> "list[str]":
    row = [None] * len(headers)
    with open(input_file) as f_in:
        line = f_in.readline()
        while line:
            words = line.split()
            try:
                variable = words[0]
                value = words[2]
                if variable in headers:
                    row[headers.index(variable)] = f"{value:>12s}"
                    data[variable].append(float(value))
                    if all(row):
                        return row
            except IndexError:
                # Not all lines will have potential variables/values, so just pass
                pass

            line = f_in.readline()

    # Only reach here if we run out of lines without finding a value for each variable
    raise RuntimeError(
        f"One or more criteria missing, was looking for {headers} but found {row}"
    )


if __name__ == "__main__":
    input_data = sys.argv[1]
    input_values = json.load(open(sys.argv[2], "r", encoding="utf-8"))

    if "report_criteria" in input_values["format"]:
        data = parse_reports(input_data)
        for variable, column in data.items():
            plot(variable, column)
    else:
        columns = load(input_data)
        for column in columns:
            plot(column[0], column[1:])
