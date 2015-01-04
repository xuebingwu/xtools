
for f in /lab/bartel3_ata/agarwal/papers/human/analysis_files/targetsites_genelevel/GSM*
do 
 fname=`basename $f`
 uniqLines -i $f -o $fname.uniq -c 1
 awk '{print $5$6$7"\t"$2}' $fname.uniq > $fname.PKA.input.txt
 PKA $fname.PKA.input.txt -weighted -o $fname
done

