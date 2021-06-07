#!/bin/bash
pimc_bin="/users/e/c/ecasiano/pimc/pimc/Build_VACC/pimc/pigsl.e";

pimc_args="-D 1 -L 4 -N 4 -l 2 -U 3.3 --mu 6.0 --subgeometry square --measurement-frequency 1 --beta 1.0 --rng boost_mt19937 --num-replicas 2 --bin-size 1000 --bins-wanted 1000"

command_file="pimc_commands";
for i in {0..499}
do
    cmd="${pimc_bin} ${pimc_args} --seed ${i}";
    echo "${cmd}" >> $command_file ;
done
