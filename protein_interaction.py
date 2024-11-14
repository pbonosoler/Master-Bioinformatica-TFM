protein_interaction = {}
unique_interactions = set()

with open("/home/pbonosoler/TFM/glycerol/counts/ProteinInteraction.txt", "r") as file:
	for line in file:
		line = line.strip('\n').split('\t')
		protein1 = line[0]
		protein2 = line[1]
		if (protein1,protein2) in unique_interactions or (protein2,protein1) in unique_interactions:
			continue
		else:
			unique_interactions.add((protein1,protein2))
        	
			if protein1 in protein_interaction.keys():
				protein_interaction[protein1] += 1
			else:
				protein_interaction[protein1] = 1

			if protein2 in protein_interaction.keys():
				protein_interaction[protein2] += 1
			else:
				protein_interaction[protein2] = 1


# Write the dictionary to the new file
with open("/home/pbonosoler/TFM/glycerol/counts/protein_interactions.tsv", "w") as output_file:        
    for gene, interactions in protein_interaction.items():
        output_file.write(f"{gene}\t{interactions}\n")
