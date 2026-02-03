#!/bin/bash

### match PL number from plasmid seq to fasta file and copy to relevant folder
### run like this
### match_plasmid_fasta_refs.bash <PL_to_fasta.tsv> folder_name_fasta_refs

fastaRefs=$2

while IFS=$'\t' read -r -a myArray
do

	# match destination folder that contains the relevant PL number
	d=$(find . -type d -iname "${myArray[0]}*" -exec readlink -f {} \;)
	f="${myArray[1]}"

	# check if column 2 has an apropriate fasta file ending in either .fa or .fasta
	if [ "${f##*.}" = "fa" ] || [ "${f##*.}" = "fasta" ]; then

		# find the location of the fasta file
		r=$(find "${fastaRefs}" -type f -iname "${f}*")

		echo "${r} ${d}"

		# if user misspelled the fasta file and therefore it cannot be found, add note and touch na
		if [ "${r##*.}" = "fa" ] || [ "${r##*.}" = "fasta" ]; then
			# copy fasta file to PL destination folder
			cp "${r}" "${d}"
		else
			touch "${d}/na"
			touch "${d}/FASTA_REF_PROB_MISSPELLED"
		fi
	else
	# If no appropriate fasta file was found, add "na" empty file to the folder for the next step
	touch "${d}/na"
	fi	

done < "$1"