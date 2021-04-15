
#!/bin/bash
pimc_bin="$HOME/pimc/pimc/pimc.e";

pimc_args="-D 1 -L 4 -N 4 -l 2 --bin-size 10000 --bins-wanted 1000 -U 3.3 --mu 6.0 --subgeometry square --num-replicas 1 --sweeps 50000000 --measurement-frequency 25 --beta 4.0 --seed 2"

command_file="pimc_commands";
for i in {0..499}
do
    cmd="${pimc_bin} ${pimc_args} --seed ${i}";
    echo "${cmd}" >> $command_file ;
done
