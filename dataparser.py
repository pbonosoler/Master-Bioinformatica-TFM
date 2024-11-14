directories_files = {}

with open("/mnt/Drive-D/DATA/Sc_RNAseq_evol_data_23.10.05.csv", "r") as file:
    for line in file:
        line = line.strip('\n').split('\t')

        if line[4] == "YPD" and line[6] == "Glycerol" and line[5] == "200":
            directory = line[0]
            if directory in directories_files.keys():
                directories_files[directory].append(line[1])
            else:
                directories_files[directory] = []
                directories_files[directory].append(line[1])

        elif line[4] == "Glycerol" and line[5] == "200":
            directory = line[0]
            if directory in directories_files.keys():
                directories_files[directory].append(line[1])
            else:
                directories_files[directory] = []
                directories_files[directory].append(line[1])

        elif line[4] == "YPD" and line[5] == "200" and line[6] == "YPD":
            directory = line[0]
            if directory in directories_files.keys():
                directories_files[directory].append(line[1])
            else:
                directories_files[directory] = []
                directories_files[directory].append(line[1])

output = "/mnt/Drive-E/seq/Sc_glycerol_t200/files_list.txt"
with open(output, 'w') as out:
        for directory, filenames in directories_files.items():
            for filename in filenames:
                out.write("/mnt/Drive-D/DATA/" + directory + "/Raw_data/" + filename+ "\n")

