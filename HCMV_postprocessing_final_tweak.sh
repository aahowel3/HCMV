for i in {1..25}
do
cd replicate_${i}_02

for file in congenital_plasma*.100.output.ms
do
base=$(basename $file .100.output.ms)
strip="${base#*.}"

Rscript ../HCMV_postprocessing_final_tweak.R $file congenital_plasma.${strip}.1000.output.ms congenital_urine.${strip}.100.output.ms congenital_urine.${strip}.1000.output.ms  >> goodgoodcombos2_tweak.csv
#Rscript HCMV_postprocessing_final.R 
done
cd ../
done
