#!/bin/bash
pimc_bin="/users/e/c/ecasiano/pimc/pimc/Build_VACC/pimc/pigsl.e";

pimc_args="-D 1 -L 16 -N 16 -l 8 -U 0.500000 --mu 6.0 --subgeometry square --measurement-frequency 1 --beta 6.0 --rng boost_mt19937 --num-replicas 2 --bin-size 10000 --bins-wanted 100"

command_file="pimc_commands";
for i in {0..164}
do
    cmd="${pimc_bin} ${pimc_args} --seed ${i}";
    echo "${cmd}" >> $command_file ;
done

pimc_args="-D 1 -L 16 -N 16 -l 8 -U 0.730000 --mu 6.0 --subgeometry square --measurement-frequency 1 --beta 6.0 --rng boost_mt19937 --num-replicas 2 --bin-size 10000 --bins-wanted 100"

command_file="pimc_commands";
for i in {0..164}
do
    cmd="${pimc_bin} ${pimc_args} --seed ${i}";
    echo "${cmd}" >> $command_file ;
done

pimc_args="-D 1 -L 16 -N 16 -l 8 -U 1.065800 --mu 6.0 --subgeometry square --measurement-frequency 1 --beta 6.0 --rng boost_mt19937 --num-replicas 2 --bin-size 10000 --bins-wanted 100"

command_file="pimc_commands";
for i in {0..164}
do
    cmd="${pimc_bin} ${pimc_args} --seed ${i}";
    echo "${cmd}" >> $command_file ;
done

pimc_args="-D 1 -L 16 -N 16 -l 8 -U 1.556100 --mu 6.0 --subgeometry square --measurement-frequency 1 --beta 6.0 --rng boost_mt19937 --num-replicas 2 --bin-size 10000 --bins-wanted 100"

command_file="pimc_commands";
for i in {0..164}
do
    cmd="${pimc_bin} ${pimc_args} --seed ${i}";
    echo "${cmd}" >> $command_file ;
done

pimc_args="-D 1 -L 16 -N 16 -l 8 -U 2.272000 --mu 6.0 --subgeometry square --measurement-frequency 1 --beta 6.0 --rng boost_mt19937 --num-replicas 2 --bin-size 10000 --bins-wanted 100"

command_file="pimc_commands";
for i in {0..164}
do
    cmd="${pimc_bin} ${pimc_args} --seed ${i}";
    echo "${cmd}" >> $command_file ;
done

pimc_args="-D 1 -L 16 -N 16 -l 8 -U 3.300000 --mu 6.0 --subgeometry square --measurement-frequency 1 --beta 6.0 --rng boost_mt19937 --num-replicas 2 --bin-size 10000 --bins-wanted 100"

command_file="pimc_commands";
for i in {0..164}
do
    cmd="${pimc_bin} ${pimc_args} --seed ${i}";
    echo "${cmd}" >> $command_file ;
done
