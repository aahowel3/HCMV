for i in {2..10}
do 
(
#mkdir replicate_"$i"
cd replicate_${i}_04

        for x in congenital*.output.ms
        do
        base=$(basename $x .output.ms)
        Rscript ../99.7filter.R $x ${base}.sub.output.ms
        done

        for x in congenital*.output.fix
        do
        base=$(basename $x .output.fix)
        mv $x ${base}.sub.output.fix
        done

        #after the initial files are generated  run them through an initial stats calcuator to see how many sites pass the 2% filter
        for x in *1000.sub.output.ms
        do
        base=$(basename $x .output.ms)
        python ../sc2-full2_stats-v0.2_1.py $base m 23430
        done

        for x in *100.sub.output.ms
        do
        base=$(basename $x .output.ms)
        echo $x
        python ../sc2-full2_stats-v0.2_1_for100.py $base m 23430
        done

        #if there are no seg sites in an ms file the sc2.py script will throw an error and not make a csv - but you need it to line up all>
        ls congenital_plasma*100.*csv > congenital_plasma_100.check
        ls congenital_plasma*1000.*csv > congenital_plasma_1000.check
        ls congenital_urine*100.*csv > congenital_urine_100.check
        ls congenital_urine*1000.*csv > congenital_urine_1000.check
        #compare the number of csvs youve actually output to what should be there
        #full2 is renamed so the output is .sub.output.csv to match new subsampled ms files
        comm -3 congenital_plasma_100.check ../congenital_plasma_100_full2.check > missing.list
        comm -3 congenital_plasma_1000.check ../congenital_plasma_1000_full2.check >> missing.list
        comm -3 congenital_urine_100.check ../congenital_urine_100_full2.check >> missing.list
        comm -3 congenital_urine_1000.check ../congenital_urine_1000_full2.check >> missing.list

	while read -r line
        do
        base=$(basename $line .output.csv)
        yes "${base},full,0.02,1,23500,filtSNPs,0" | head -n 1700 > $line
        done < missing.list

        #create the 4 column file that the R script will use to check which param combos have all 4 output files +/- 100snps of the B103 numbers
        for x in congenital_plasma*100.sub.output.csv
        do
        sed -n 1623p $x >> congenital.plasma.100.check
        done
        for x in congenital_plasma*1000.sub.output.csv
        do
        sed -n 1627p $x >> congenital.plasma.1000.check
        done
        for x in congenital_urine*100.sub.output.csv
        do
        sed -n 1623p $x >> congenital.urine.100.check
        done
        for x in congenital_urine*1000.sub.output.csv
        do
        sed -n 1627p $x >> congenital.urine.1000.check
        done

###        paste congenital.plasma.100.check congenital.plasma.1000.check congenital.urine.100.check congenital.urine.1000.check | column -s "," > checkcombos.csv

        #get out the list of parameter combinations you want to go to the trouble of filtering their ms files and re-running their stats using only the 2% passing
	Rscript ../twopercentfilter.R

        #do that on the good list of parameter combinations - need a diff sc2.py for 100s
        while read -r line
        do
        Rscript ../HCMV_filtering_forterbotscriptinput.R "congenital_plasma.${line}.100.sub.output.ms" "congenital_plasma.${line}.100.filter.output.ms"
        cp congenital_plasma.${line}.100.sub.output.fix congenital_plasma.${line}.100.filter.output.fix
        python ../sc2-full2_stats-v0.2_1_for100.py congenital_plasma.${line}.100.filter m 23500
        done < goodcombos2.csv

        while read -r line
        do
        Rscript ../HCMV_filtering_forterbotscriptinput_for1000.R "congenital_plasma.${line}.1000.sub.output.ms" "congenital_plasma.${line}.1000.filter.output.ms"
        cp congenital_plasma.${line}.1000.sub.output.fix congenital_plasma.${line}.1000.filter.output.fix
        python ../sc2-full2_stats-v0.2_1.py congenital_plasma.${line}.1000.filter m 23500
        done < goodcombos2.csv

        while read -r line
        do
        Rscript ../HCMV_filtering_forterbotscriptinput.R "congenital_urine.${line}.100.sub.output.ms" "congenital_urine.${line}.100.filter.output.ms"
        cp congenital_urine.${line}.100.sub.output.fix congenital_urine.${line}.100.filter.output.fix
        python ../sc2-full2_stats-v0.2_1_for100.py congenital_urine.${line}.100.filter m 23500
        done < goodcombos2.csv

        while read -r line
        do
        Rscript ../HCMV_filtering_forterbotscriptinput_for1000.R "congenital_urine.${line}.1000.sub.output.ms" "congenital_urine.${line}.1000.filter.output.ms"
        cp congenital_urine.${line}.1000.sub.output.fix congenital_urine.${line}.1000.filter.output.fix
        python ../sc2-full2_stats-v0.2_1.py congenital_urine.${line}.1000.filter m 23500
        done < goodcombos2.csv
) &
done
