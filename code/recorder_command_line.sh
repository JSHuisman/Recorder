#!/bin/bash

small_number_of_cores=3

for conj_rate in {1..37}
do
    for mig_rate in {1..37}
    do

        if [ ! -f "../simulations/${conj_rate}_${mig_rate}.csv" ]; then
            # Control will enter here if $File doesn't exist.
            #echo $conj_rate $mig_rate
            Rscript --vanilla --slave recorder_command_line.R $small_number_of_cores $conj_rate $mig_rate
        fi

    done
done