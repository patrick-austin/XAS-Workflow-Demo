import csv
import sys


SP_HEADER = "id,filename,label,s02,e0,sigma2,deltar"
GDS_HEADER = "id,name,value,expr,vary"


def main(paths_file: str):
    # with open(paths_file) as file:
    #     file.readline()
    #     line = file.readline()
    #     while line:
    #         values = line.split(",")
    gds_data = []
    sp_data = []
    row_id = 1

    with open(paths_file) as file:
        reader = csv.reader(file)
        header = True
        for row in reader:
            if header:
                header = False
                continue

            selected = int(row[-1])
            if selected:
                filename = row[0].strip()
                # path_id = file[-8:-4]
                label = row[-2].strip()
                sigma2 = 's' + label.lower()
                sp_data.append(f"{row_id:d},{filename:s},{label:s},amp,enot,{sigma2:s},alpha*reff")
                gds_data.append(f"{row_id:d},{sigma2:s},0.003,,True")
                row_id += 1

    gds_data.append(f"{row_id:d},alpha,1e-07,,True")
    gds_data.append(f"{row_id + 1:d},amp,1.0,,True")
    gds_data.append(f"{row_id + 2:d},enot,1e-07,,True")

    with open("sp.csv", "w") as out:
        out.writelines(sp_data)

    with open("gds.csv", "w") as out:
        out.writelines(gds_data)


if __name__ == "__main__":
    paths_file = sys.argv[1]
    main(paths_file)
