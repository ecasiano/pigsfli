
#!/bin/bash
pimc_bin="$HOME/PATH/TO/pimc.e";

pimc_args="-U 100.000 -D 1 -L 10 -N 10 -l 5 --beta 4.0 --bin-size 1000 --bins-wanted 1000 --mu 6 --subgeometry  square";

command_file="pimc_commands";
for i in {0..499}
do
    cmd="${pimc_bin} ${pimc_args} --seed ${i}";
    echo "${cmd}" >> $command_file ;
done
