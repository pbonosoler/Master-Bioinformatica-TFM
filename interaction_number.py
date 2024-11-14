gene_n_interactions = {}


# Read and process input
with open("/home/pbonosoler/TFM/glycerol/counts/GeneInteraction.txt", "r") as gi:
    for line in gi:
        if line.startswith("Query"):
            continue
        else:
            line = line.split("\t")
            name = line[0]
            name2 = line[2]

            if name in gene_n_interactions.keys():
                gene_n_interactions[name] += 1
            else:
                gene_n_interactions[name] = 1

            if name2 in gene_n_interactions.keys():
                gene_n_interactions[name2] += 1
            else:
                gene_n_interactions[name2] = 1


# Write the dictionary to the new file
with open("/home/pbonosoler/TFM/glycerol/counts/gene_interactions.tsv", "w") as output_file:        
    for gene, interactions in gene_n_interactions.items():
        output_file.write(f"{gene}\t{interactions}\n")
