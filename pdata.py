directories_files = []

with open("/mnt/Drive-D/DATA/Sc_RNAseq_evol_data_23.10.05.csv", "r") as file:
    for line in file:
        line = line.strip('\n').split('\t')

        directory = line[0]
        sline = line[2]
        medium = line[4]
        time = line[5]
        challenge = line[6]
        filename = line[1]

        if (medium == "YPD" and time == "200" and challenge == "Glycerol") or (medium == "Glycerol" and time == "200") or (medium == "YPD" and time == "200" and challenge == "YPD"):
            directories_files.append((filename, sline, medium, time, challenge))

output = "/mnt/Drive-E/seq/Sc_glycerol_t200/pdata.txt"
with open(output, 'w') as out:
    # Write column names
    out.write("line\tevolved\ttime\tchallenge\n")
    # Write each row
    for filename, line, medium, time, challenge in directories_files:
        out.write(f"{filename}\t{line}\t{medium}\t{time}\t{challenge}\n")

