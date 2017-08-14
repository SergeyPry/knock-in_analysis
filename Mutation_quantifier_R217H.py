import re
import sys

filename = sys.argv[-1]

# The program is for quantifying the different types of CRISPR-based ssODN knock-in 
# mutation outcomes

# 0. Define some important sequence variables for matching patterns in the SAM alignment files

# The following patterns are for the tp53 R143H knock-in
HDR_pattern =  "CATGAACCACC"
WT_pattern =   "GATGAACCGCA"
HDR_sgRNA_pattern =  "CTGCATGGGGGGCATGAACCACC"
WT_sgRNA_pattern =   "CTGCATGGGGGGGATGAACCGCA"
current_ref = "R217_reference"

# categories to which to assign the alignments
CRISPR_cats = {}

CRISPR_cats["wt"] = 0
CRISPR_cats["del noHDR"] = 0
CRISPR_cats["ins noHDR"] = 0
CRISPR_cats["del HDR"] = 0
CRISPR_cats["ins HDR"] = 0
CRISPR_cats["correct HDR"] = 0
CRISPR_cats["unmapped HDR"] = 0 
CRISPR_cats["unmapped wt"] = 0
CRISPR_cats["unmapped other"] = 0


# Open files for writing the reads for alignments
del_HDR_handle = str(filename) + "_R217H_del_HDR.fa"
del_HDR_file = open(del_HDR_handle, "a")

ins_HDR_handle = str(filename) + "_R217H_ins_HDR.fa"
ins_HDR_file = open(ins_HDR_handle, "a")

corrHDR_handle = str(filename) + "_R217H_correct_HDR.fa"
correct_HDR_file = open(corrHDR_handle, "a")


Counter = 0


# 1. Reading each line of the SAM file
with open(filename) as f:
    for line in f:
		Counter += 1
		line_data = line.split("\t")
		
		# assign relevant parts of the data to variables
		ref = line_data[2]
		CIGAR = line_data[5]
		cigar_matches = re.findall(r'(\d+)([A-Z]{1})', CIGAR)
		seq = line_data[9]
		
		## Determining the category

		# 1. Divide the alignments into mapped and unmapped
		
		if ref == current_ref: # most of the reads will have this property
			# case 1: the read is wild-type in the relevant part of the read
			# meaning that both the relevant targeted sites of the read were not modified
			# and no mutations were generated at the site of sgRNA cutting
			if re.search(WT_sgRNA_pattern, seq):
				CRISPR_cats["wt"] += 1;
			else:
			# case 2: the read contains an HDR event
				if re.search(HDR_pattern, seq):
					# case 3: a fully correct HDR
					if re.search(HDR_sgRNA_pattern, seq):
						CRISPR_cats["correct HDR"] += 1
						# writing this read out to a file
						correct_HDR_file.write(">" + "correct_HDR_" + CIGAR + "-" + str(Counter) + "\n" + seq + "\n")
						
					# case 4: HDR events with indel mutations
					else:
						# find insertion or deletion in the CIGAR string
						pos = 0                 # current location
						netIndel = 0            # net indel length
						
						# quantify the total net insertions or deletions in the target region 
						# 
						for m in cigar_matches:
							pos += int(m[0])
							
							if pos >= 100:
								break
							else:
								if pos > 30 and m[1] == 'I':
									netIndel += int(m[0])
								elif pos > 30 and m[1] == 'D':
									netIndel -= int(m[0])
						
						# classify the corresponding indel events:
						if netIndel > 0:
							CRISPR_cats["ins HDR"] += 1
							# writing this read out to a file
							ins_HDR_file.write(">" + "ins_HDR_" + CIGAR + "-" + str(Counter) + "\n" + seq + "\n")
							
						elif netIndel < 0:
							CRISPR_cats["del HDR"] += 1
							# writing this read out to a file
							del_HDR_file.write(">" + "del_HDR_" + CIGAR + "-" + str(Counter) + "\n" + seq + "\n")
					
				# case 5: any other indel events
				else: 

				# quantify the total net insertions or deletions in the target region 
					pos = 0                 # current location
					netIndel = 0            # net indel length
					
					for m in cigar_matches:
						pos += int(m[0])
							
						if pos >= 100:
							break
						else:
							if pos > 30 and m[1] == 'I':
								netIndel += int(m[0])
							elif pos > 30 and m[1] == 'D':
								netIndel -= int(m[0])
						
					# classify the corresponding indel events:
					if netIndel > 0:
						CRISPR_cats["ins noHDR"] += 1
						
					elif netIndel < 0:
						CRISPR_cats["del noHDR"] += 1

				
		elif ref == "*": # these reads are unmapped
			if re.search(HDR_pattern, seq):
				# assign to a category
				CRISPR_cats["unmapped HDR"] += 1
				
								
			elif re.search(WT_pattern, seq):
				# assign to a category
				CRISPR_cats["unmapped wt"] += 1
								
			else: 
				CRISPR_cats["unmapped other"] += 1

# close all files
del_HDR_file.close()
ins_HDR_file.close()
correct_HDR_file.close()
				
# sum of all events
sum = 0

for k in CRISPR_cats:
	sum += CRISPR_cats[k]

# open a file for writing the summary results
 
outfile = str(filename) + "_summary.txt"
results = open(outfile, "w")

header = "type" + "\t" + "count" + "\n"
results.write(header)

results.write("del - no HDR" + "\t")
results.write(str(CRISPR_cats["del noHDR"]) + "\n")

results.write("ins - no HDR" + "\t")
results.write(str(CRISPR_cats["ins noHDR"]) + "\n")

results.write("correct HDR" + "\t")
results.write(str(CRISPR_cats["correct HDR"]) + "\n")

results.write("del - HDR" + "\t")
results.write(str(CRISPR_cats["del HDR"]) + "\n")

results.write("ins - HDR" + "\t")
results.write(str(CRISPR_cats["ins HDR"]) + "\n")

results.write("unmapped - HDR" + "\t")
results.write(str(CRISPR_cats["unmapped HDR"]) + "\n")

results.write("unmapped - wt" + "\t")
results.write(str(CRISPR_cats["unmapped wt"])+ "\n")

results.write("unmapped - other" + "\t")
results.write(str(CRISPR_cats["unmapped other"]) + "\n")

results.write("wild-type" + "\t")
results.write(str(CRISPR_cats["wt"]) + "\n")

				
# summary and results output
print("The number of wt events is " + str(CRISPR_cats["wt"]))
print("The number of deletion events without HDR is " + str(CRISPR_cats["del noHDR"]))
print("The number of insertion events without HDR is " + str(CRISPR_cats["ins noHDR"]))
print("The number of correct HDR events is " + str(CRISPR_cats["correct HDR"]))
print("The number of HDR events with deletions is " + str(CRISPR_cats["del HDR"]))
print("The number of HDR events with insertions is " + str(CRISPR_cats["ins HDR"]))
print("The number of unmapped HDR events is " + str(CRISPR_cats["unmapped HDR"]))
print("The number of unmapped wt events is " + str(CRISPR_cats["unmapped wt"]))
print("The number of other unmapped events is " + str(CRISPR_cats["unmapped other"]))

print("The total number of events in all categories is "  + str(sum))

		


