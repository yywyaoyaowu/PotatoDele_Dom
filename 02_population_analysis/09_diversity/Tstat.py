
with open("gene_in_top0.05_2", "r") as f:
	data = f.readlines()

region_dict = {}
for i in range(len(data)):
	line = data[i].split()
	if line[0] not in region_dict:
		region_dict[line[0]] = []

for i in range(len(data)):
	line = data[i].split()
	region_dict[line[0]].append(line[7])

with open("gene_in_top0.05_merge.txt", "w") as f:
	for key in region_dict:
		gene_string = ','.join(map(str, region_dict[key]))
		f.writelines(str(key) + "\t" + gene_string + "\t" + str(len(region_dict[key])) + "\n")



