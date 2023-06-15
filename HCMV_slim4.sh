#!/bin/bash

for i in {35..40}
do
mkdir replicate_"$i"
cd replicate_"$i" 

	#generate all .ms and .fix files according to the list of param combinati9ons
	while IFS="," read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12
        	do
	        slim -d DFE="$col1" -d f0="$col2" -d f1="$col3" -d f2="$col4" -d f3="$col5" -d recombrate="$col6" -d recomb="$col7" -d murate="$col8" -d mu="$col9" -d progeny="$col10" -d burnin="$col11" -d gr="$col12" ../HCMV_congenital_final.slim 
	done < ../params_congenital.txt

	#after the initial files are generated  run them through an initial stats calcuator to see how many sites pass the 2% filter
	for x in *1000.output.ms
	do
	base=$(basename $x .output.ms)
	echo $x
	python ../sc2-full2_stats-v0.2_1.py $base m 23500
	done

	for x in *100.output.ms
	do
	base=$(basename $x .output.ms)
	echo $x
	python ../sc2-full2_stats-v0.2_1_for100.py $base m 23500
	done
	
	#if there are no seg sites in an ms file the sc2.py script will throw an error and not make a csv - but you need it to line up all the files that pass
	ls congenital_plasma*100.*csv > congenital_plasma_100.check
	ls congenital_plasma*1000.*csv > congenital_plasma_1000.check
	ls congenital_urine*100.*csv > congenital_urine_100.check
	ls congenital_urine*1000.*csv > congenital_urine_1000.check
	#compare the number of csvs youve actually output to what should be there
	comm -3 congenital_plasma_100.check ../congenital_plasma_100_full.check > missing.list
	comm -3 congenital_plasma_1000.check ../congenital_plasma_1000_full.check >> missing.list
	comm -3 congenital_urine_100.check ../congenital_urine_100_full.check >> missing.list
	comm -3 congenital_urine_1000.check ../congenital_urine_1000_full.check >> missing.list

	#create an 0 SNPS csv for each missing csv
	while read -r line
	do
	base=$(basename $line .output.csv)
	yes "${base},full,0.02,1,23500,filtSNPs,0" | head -n 1700 > $line
	done < missing.list

	#create the 4 column file that the R script will use to check which param combos have all 4 output files +/- 100snps of the B103 numbers
	for x in congenital_plasma*100.output.csv
	do
	sed -n 1626p $x >> congenital.plasma.100.check
	done 
	for x in congenital_plasma*1000.output.csv
        do
        sed -n 1630p $x >> congenital.plasma.1000.check
        done
	for x in congenital_urine*100.output.csv
        do
        sed -n 1626p $x >> congenital.urine.100.check
        done
	for x in congenital_urine*1000.output.csv
        do
        sed -n 1630p $x >> congenital.urine.1000.check
        done

	paste -d "," congenital.plasma.100.check congenital.plasma.1000.check congenital.urine.100.check congenital.urine.1000.check > checkcombos.csv

	#get out the list of parameter combinations you want to go to the trouble of filtering their ms files and re-running their stats using only the 2% passing 
	Rscript ../twopercentfilter.R 

	#do that on the good list of parameter combinations - need a diff sc2.py for 100s
	while read -r line
	do
	Rscript ../HCMV_filtering_forterbotscriptinput.R "congenital_plasma.${line}.100.output.ms" "congenital_plasma.${line}.100.filter.output.ms"
	cp congenital_plasma.${line}.100.output.fix congenital_plasma.${line}.100.filter.output.fix
	python ../sc2-full2_stats-v0.2_1_for100.py congenital_plasma.${line}.100.filter m 23500
	done < goodcombos2.csv

	while read -r line
        do
        Rscript ../HCMV_filtering_forterbotscriptinput_for1000.R "congenital_plasma.${line}.1000.output.ms" "congenital_plasma.${line}.1000.filter.output.ms"
        cp congenital_plasma.${line}.1000.output.fix congenital_plasma.${line}.1000.filter.output.fix
        python ../sc2-full2_stats-v0.2_1.py congenital_plasma.${line}.1000.filter m 23500
        done < goodcombos2.csv
	
	while read -r line
        do
        Rscript ../HCMV_filtering_forterbotscriptinput.R "congenital_urine.${line}.100.output.ms" "congenital_urine.${line}.100.filter.output.ms"
        cp congenital_urine.${line}.100.output.fix congenital_urine.${line}.100.filter.output.fix
        python ../sc2-full2_stats-v0.2_1_for100.py congenital_urine.${line}.100.filter m 23500
        done < goodcombos2.csv

	while read -r line
        do
        Rscript ../HCMV_filtering_forterbotscriptinput_for1000.R "congenital_urine.${line}.1000.output.ms" "congenital_urine.${line}.1000.filter.output.ms"
        cp congenital_urine.${line}.1000.output.fix congenital_urine.${line}.1000.filter.output.fix
        python ../sc2-full2_stats-v0.2_1.py congenital_urine.${line}.1000.filter m 23500
        done < goodcombos2.csv
	
	 mkdir good_combos
        while read -r line
        do
        mv *${line}.* good_combos
        done < goodcombos2.csv

        rm congenital*
	
	cd ../
done 
