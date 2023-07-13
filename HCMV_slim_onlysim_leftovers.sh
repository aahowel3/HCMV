for i in {1..25}
do
(
#mkdir replicate_"$i"
cd replicate_"$i"
        #generate all .ms and .fix files according to the list of param combinati9ons
        while IFS="," read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12
                do
        coll1=$(sed -e 's/^"//' -e 's/"$//' <<<"$col1")
	coll11=$(sed -e 's/^"//' -e 's/"$//' <<<"$col11")
	coll7=$(sed 's/./&0/5' <<<"$col7") 
	repls=2.0
	coll9="${col9/2/$repls}"
	colll9=$(sed 's/./&0/5' <<<"$coll9")
	coll10=$(sed -e 's/[0]*$//g' <<<"$col10")
	fname="congenital_plasma."${coll1}".${colll9}.${coll7}.${coll10}.${col12}.${coll11}.100.output.ms"
        if [[ -f "$fname" ]]
        then
                continue
        fi 
        	slim -d DFE="$col1" -d f0="$col2" -d f1="$col3" -d f2="$col4" -d f3="$col5" -d recombrate="$col6" -d recomb="$col7" -d murate="$col8" -d mu="$col9" -d progeny="$col10" -d burnin="$col11" -d gr="$col12" ../HCMV_congenital_final.slim
        done < ../xaa_nor.txt

) &
done
