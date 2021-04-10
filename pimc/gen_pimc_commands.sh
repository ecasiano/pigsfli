
#!/bin/bash
pimc_bin="$HOME/pimc/pimc/pimc.e";

pimc_args="-U 3.3 -D 1 -L 16 -N 16 -l 8 --beta 7 --bin-size 1000 --bins-wanted 1000 --mu 6 --subgeometry  square";

command_file="pimc_commands";
for i in {0..499}
do
    cmd="${pimc_bin} ${pimc_args} --seed ${i}";
    echo "${cmd}" >> $command_file ;
done
