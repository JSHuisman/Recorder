#!/bin/bash

small_number_of_cores=3

for carrying_cap in {1..29}
do

        if [ ! -f "../simulations_cc/${conj_rate}_${mig_rate}.csv" ]; then
            # Control will enter here if $File doesn't exist.
            #echo $conj_rate $mig_rate
            Rscript --vanilla --slave carrying_cap_command_line.R $small_number_of_cores $carrying_cap
        fi

done