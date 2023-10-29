#!/bin/bash


#Joins the Tco sample count files for each group vertically

for dir in ~/analysis/groups/*dir; do
        group=$(basename -s _dir $dir) #Extracts the group name of each directory e.g Clone1_24_Induced. -s passes the basename of the file minus the  _counts suffix
        cd $dir
	cat $dir/* | cut -f 4,5,6 | sort -k1,1 >  "$group"_counts
done



#Calculates the mean expression for each gene in each treatment group

for file in ~/analysis/groups/*dir/*counts; do
	group=$(basename -s _counts $file)
	cat $file | awk -F"\t" '{genetotal[$1] += $3
				count[$1]++
				genedesc[$1] = $2;} #Sum the total gene expression for each gene and increment count of records (i.e. number of replicates) for each gene
				END {
					for (gene in genetotal) #Prints the average for each gene along with the gene name and description
						print gene"\t"genedesc[gene]"\t"genetotal[gene]/count[gene];}' |
	sort -k1,1 > "$group"_avg
	echo "Mean gene expression values calculated for ${group} treatment group"
done


#Calculates the fold change in gene expression over time for each Tco culture e.g. Clone1 and of each induction state e.g. Induced

for file in ~/analysis/groups/*dir/*avg; do
        culture=$(basename $file | awk -F"_" '{print $1}') #Extracts the culture group of each file
	type=$(basename $file | awk -F"_" '{print $3}') #Extracts the induction state of each file
        paste ~/analysis/groups/*dir/"$culture"_0_Uninduced_avg ~/analysis/groups/*dir/"$culture"_24_"$type"_avg ~/analysis/groups/*dir/"$culture"_48_"$type"_avg |
	cut -f 1,2,3,6,9  > ~/analysis/foldcalc/"$culture"_"$type"_avg #Horizontally joins the 3 files containing gene expression averages at different time points
	cat ~/analysis/foldcalc/"$culture"_"$type"_avg | awk -F"\t" '{if ($3 != 0 && $4 != 0 && $5 != "0") {print $0}}' | #Only includes averages not equal to 0
	awk -F"\t" '{
			nohrs[$1] += $3 #The average for each gene at 0hrs 
			twofourhrs[$1] += $4 #The average for each gene at 24hrs
			foureighthrs[$1] += $5 #The average for each gene at 48hrs
			genedesc[$1] = $2;}
			END {
				for (gene in nohrs)
					print gene"\t"twofourhrs[gene]/nohrs[gene]"\t"foureighthrs[gene]/twofourhrs[gene]"\t"foureighthrs[gene]/nohrs[gene]"\t"genedesc[gene];}' |
	sort -k4,4nr | #Calculates fold changes for each gene and sorts 0-48hrs fold changes in order of decreasing magnitude
	awk -F"\t" 'BEGIN{print "gene\t0-24hrs\t24-48hrs\t0-48hrs\tdescription"}1' > ~/analysis/foldcalc/"$culture"_"$type"_foldchange #Adds a header row for easy viewing
	cat ~/analysis/foldcalc/"$culture"_"$type"_avg | awk -F"\t" '{if ($3 == 0 || $4 == 0 || $5 == 0) {print $0}}' > ~/analysis/foldcalc/"$culture"_"$type"_avg_0
	#Outputs rows containing an average value of 0 into a separate folder - genes were excluded from fold change calculations
	echo "fold changes in gene expression over time for ${type} ${culture} Tco samples generated - see ${culture}_${type}_foldchange in the culture directory"
	echo "NOTE: ${type} ${culture} Tco genes with an average expression of 0 at any time point were omitted - these data can be viewed in ${culture}_${type}_avg_0"
done



