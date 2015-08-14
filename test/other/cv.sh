	for k in 3 
	do 
	    for shift in 1
	    do
	    echo $k, $shift
	    rm log output.* data.*; cross_validation -data input.weighted.txt  -train_cmd "PKA data.train.@ -p 2 -FDR -o output.@ -max_k $k -max_shift $shift -weighted 2>> log " -test_cmd "PKA data.test.@ -predict output.@ -weighted 2>> log" -nfold 5
        rm tmp
	    cat output.*.score | grep -v inf > tmp
	    Rscript cor.r
	    done 
	done
