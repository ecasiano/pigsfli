# Documentation  

## Introduction 

This webpage contains the details of a worm algorithm path integral quantum Monte Carlo (WA-PIMC) code actively developed in c++ since 2020 in the [Del Maestro group](http://delmaestro.org/adrian) based on:

- T = 0: [N.V. Prokof'ev, B.V. Svistunov and I.S. Tupitsyn, J. Exp. Theor <b>113</b>, 1366 (1997)](https://arxiv.org/pdf/cond-mat/9703200.pdf)

It can be used to simulate the [Bose-Hubbard Model](https://en.wikipedia.org/wiki/Bose–Hubbard_model) in one, two and three spatial dimensions. In it's current version, it's main application is the calculation of Renyi Entanglement Entropies to quantify spatial entanglement between two spatial subregions of a lattice.

If you have questions, please contact me at ecasiano@vols.@utk.edu.

The development and maintenance of this code base has been supported in part by the National Science Foundation under Awar Nos. [DMR-1553991](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1553991).

[<img width="100" src="https://www.nsf.gov/images/logos/NSF_4-Color_bitmap_Logo.png">](http://www.nsf.gov/awardsearch/showAward?AWD_ID=1553991)


## Installation

This program has been successfully compiled and run on both Intel and AMD systems using clang, g++, pathscale and icpc. Before installing, one needs to ensure that all dependencies are met.  We recommend that the required libraries (boost) are installed in a `local` folder inside your home directory: `$HOME/local`.

## Dependencies 

## Path Integral Monte Carlo

After successfully installing the dependencies (which I think is only boost and the header librarires) you are now ready to compile the
main pimc program on your system.
PIMC uses CMake for build, test and installation automation. For details on using CMake consult https://cmake.org/documentation/. In short, the following steps should work on UNIX-like systems:

  ```
  mkdir build
  cd build
  cmake ../src
  make
  sudo make install
  ```

On Windows try:

  ```
  md build
  cd build
  cmake ../src
  cmake --build . --config Release
  cmake --build . --target install
  ```

As above, and with further details below, but you should consider using the following CMake options with the appropriate value instead of xxx :

- `-D NDIM=1|2|3` the number of spatial dimensions
- `-D CMAKE_C_COMPILER=xxx` equal to the name of the C99 Compiler you wish to use (or the environment variable `CC`)
- `-D CMAKE_CXX_COMPILER=xxx` equal to the name of the C++17 compiler you wish to use (or the environment variable `CXX`)
- `-D CMAKE_PREFIX_PATH=xxx` to add a non-standard location for CMake to search for libraries, headers or programs
- `-D CMAKE_INSTALL_PREFIX=xxx` to install pimc to a non-standard location
- `-D BOOST_ROOT=xxx` to add non-standard location for Boost install
- `-D STATIC=1` to enable a static build
- `-D CMAKE_BUILD_TYPE=Debug` to build pigsl in debug mode
- `-D CMAKE_BUILD_TYPE=PIMC` to build pimcl
- `-D CMAKE_BUILD_TYPE=PIMCDebug` to build pimcl in debug mode
- `-E env CXXFLAGS="xxx"` add additional compiler flags
- `-E env LDFLAGS="xxx"` add additional linker flags

Executables will be installed to `CMAKE_INSTALL_PREFIX` location or if the install is skiped will be located in `build/pimc`.
Executables produced are `pigsl.e`, `pigsld.e`, `pimcl.e`, and `pimcld.e` for `CMAKE_BUILD_TYPE=Release|Debug|PIMC|PIMCDebug` respectively.

If you run into problems, failures with linking etc., common errors may include
not properly setting your `LD_LIBRARY_PATH` or not starting from a clean build
directory (issue `make clean` or `rm -rf ./*` inside the build directory).

## Usage

FIXME!!!!!

In order to get a quick idea of the options which the code accepts type:
```bash
pigsl.e --help
```

The code requires various combinations of these options to run, and the help message should give
you an idea about which ones are mandatory.

### Quick Start 

If you want to perform a quick test-run for a one-dimensional chain of 4 sites and 4 particles, you can run the following command:
```bash
./pigsl.e -D 1 -L 4 -N 4 -l 2 -U 3.300000 --mu 6.0 --subgeometry square --measurement-frequency 1 --beta 1.0 --rng boost_mt19937 --num-replicas 1 --bin-size 10001 --bins-wanted 100 --seed 0
```

In order for this to work, you will need a folder named `OUTPUT` in the directory where you type the command as it will produce output files in `OUTPUT` that contain all the results of the code.  Each run of the code is associated with a unique identifying integer: the `PIMCID`.  The options used in this demo include a subset of all the possible options:

| Code Option | Description |
| :-----------: | ----------- |
|`D`     |  Dimension of hypercubic lattice |
|`L`     |  Linear size of hypercube |
|`N`     |  Total number of particles |
|`U`     |  Interaction potential |
|`l`     |  Linear size of hypercubic subregion |
|`sweeps`     |  Number of Monte Carlo sweeps |
|`beta`     |  Set length of imaginary time |
|`mu`     |  Chemical potential |
|`t`     |  Tunneling parameter |
|`canonical`     |  set to false for grand canonical simulation |
|`seed`     |  Random seed value |
|`sweeps-pre`     |  Number pre-equilibration sweeps |
|`bin-size` |  Number of measurements per bin |
|`subgeometry`     |  Shape of subregion: square OR strip|
|`num-replicas`     |  Number of replicas|
|`measurement-frequency`     |  Measurements will be performed every other this amount |
|`rng`     |  Random Number Generator type|
|`restart`     |  continue simulation from a loaded rng state |
|`no-accessible`     |  do not calculate accessible entanglement entropies|

All options, including lists of possible values and default values can be seen
by using the `--help flag`.

The output of the above command should yield something like:
```bash
uuid: 758c7ce2-8924-4982-9dd6-330e6fd4297f
sub-sites: 0 1
U: 3.3

  _           _   _   _           ______ _____ _____  _____
 | |         | | | | (_)          | ___ \_   _|  __ \/  ___|
 | |     __ _| |_| |_ _  ___ ___  | |_/ / | | | |  \/\ `--.
 | |    / _` | __| __| |/ __/ _ \ |  __/  | | | | __  `--. \
 | |___| (_| | |_| |_| | (_|  __/ | |    _| |_| |_\ \/\__/ /
 \_____/\__,_|\__|\__|_|\___\___| \_|    \___/ \____/\____/


Stage (1/3): Determining mu and eta...

mu: 6 eta: 0.5 Z-frac: 38.1003%
N     P(N)
4     *
5     ****************************************************************************************************
<N>: 4.99561

mu: 0.575612 eta: 0.207591 Z-frac: 33.9793%
N     P(N)
4     **********************************************
5     *******************************************************
<N>: 4.54745

mu: 0.385238 eta: 0.217824 Z-frac: 30.9912%
N     P(N)
3     ***
4     ***************************************************
5     ***********************************************
<N>: 4.44279

mu: -1.09811 eta: 0.220017 Z-frac: 32.5504%
N     P(N)
3     *********************************************************
4     *****************************************
5     ***
<N>: 3.46202

mu: 0.397318 eta: 0.230872 Z-frac: 24.2013%
N     P(N)
3     ******************
4     ******************************************************************
5     *****************
<N>: 3.98683

Fine tuning eta... (Want: 10% < Z-frac < 15%)

mu: 0.397318 eta: 0.220616 Z-frac: 30.7196%
N     P(N)
4     ****************************************************
5     *************************************************
<N>: 4.48507

mu: 0.457067 eta: 0.110308 Z-frac: 52.8473%
N     P(N)
4     **************************************************
5     ***************************************************
<N>: 4.50258

mu: 0.446738 eta: 0.159947 Z-frac: 32.3957%
N     P(N)
3     *****************************
4     *******************************************************************
5     *****
<N>: 3.76762

mu: 0.446738 eta: 0.0799734 Z-frac: 59.7352%
N     P(N)
4     *******************************************************
5     **********************************************
<N>: 4.45898

mu: 0.611192 eta: 0.115961 Z-frac: 46.8267%
N     P(N)
3     ***********
4     ******************************************************
5     *************************************
<N>: 4.26167

mu: 0.611192 eta: 0.168144 Z-frac: 34.4293%
N     P(N)
3     ************
4     ************************************************************
5     *****************************
<N>: 4.17189

mu: 0.611192 eta: 0.084072 Z-frac: 60.8987%
N     P(N)
3     *
4     **********************************************
5     ******************************************************
<N>: 4.52409

mu: -1.48159 eta: 0.121904 Z-frac: 42.0507%
N     P(N)
3     *********************************************
4     ***************************************************
5     ******
<N>: 3.61595

mu: -0.458574 eta: 0.121904 Z-frac: 47.8991%
N     P(N)
3     ********************************************
4     *****************************************************
5     ****
<N>: 3.59397

mu: 0.843493 eta: 0.176761 Z-frac: 26.578%
N     P(N)
3     *********************
4     *****************************************************************************
5     ****
<N>: 3.82353

mu: 0.843493 eta: 0.0883807 Z-frac: 43.504%
N     P(N)
3     **********************
4     *******************************************************************************
<N>: 3.78818

Stage (2/3): Equilibrating...

Stage (3/3): Main Monte Carlo loop...

-------- Detailed Balance --------

Insert Worm: 8715/678452
Delete Worm: 7888/8459

Insert Anti: 8346/526521
Delete Anti: 7561/13666

InsertZero Worm: 126201/1421550
DeleteZero Worm: 127192/136187

InsertZero Anti: 125879/1041305
DeleteZero Anti: 126870/225080

InsertBeta Worm: 127788/1417684
DeleteBeta Worm: 128409/138049

InsertBeta Anti: 125852/1045108
DeleteBeta Anti: 126473/225221

Advance Head: 325136/328752
Recede  Head: 325192/325705

Advance Tail: 325756/326279
Recede  Tail: 325976/329340

IKBH: 292313/757246
DKBH: 291279/328708

IKAH: 175833/756466
DKAH: 176165/208877

IKBT: 176453/582403
DKBT: 176670/209532

IKAT: 291004/756383
DKAT: 291485/329540

SWAP: 0/0
UNSWAP: 0/0

SWAP Advance Head: 0/0
SWAP Recede Head: 0/0

SWAP Advance Tail: 0/0
SWAP Recede Tail: 0/0

beta: 1

sweeps: 1e+07
Z_ctr: 3101
Z_frac: 58.2785% (3101/5321)

<N>: 4.53563

Elapsed time: 7.34362 seconds
```

The above shows the process in which $\mu$ is calibrated such that we get an average number of particles close to our chosen $N$. In stage 2, $\eta$ is calibrated such that we obtain a diagonal fraction of about 45%. And finally, stage 3 is composed of a sub-stage where the system runs for a while for equilibration, then after that, measurements start to be taken until the number of bins desired are obtained.

To analyze the results of the simulation, you will need to use a number of python programs located in the `scripts` directory which can be obtained via:

### Output

The results of running the code are a number of data and state files that live in the directory where the code is ran.  If the code is run for more than one replica, there will be an additional files related to the computation of Renyi entanglement entropies.  The generic output files are:

| Output File | Description |
| ----------- | ----------- |
|`D_L_N_l_U_t_beta_bin-V_seed_subgeometry_num-replicas.dat` | potential energy |
|`D_L_N_l_U_t_beta_bin-K_seed_subgeometry_num-replicas.dat` | kinetic energy |
|`D_L_N_l_U_t_beta_bin-SWAP_seed_subgeometry_num-replicas.dat` | times that different number of SWAP kinks were measured |
|`D_L_N_l_U_t_beta_bin-SWAP_seed_subgeometry_num-replicas.dat` | times that different number of SWAP kinks were measured |
|`D_L_N_l_U_t_beta_bin-SWAP_seed_subgeometry_num-replicas.dat` | times that different number of SWAP kinks were measured |
|`D_L_N_l_U_t_beta_bin-size_system-state_seed_subgeometry_num-replicas.dat` |  state of the system at which restarted simulation will begin |
|`D_L_N_l_U_t_beta_bin-rng-state_seed_subgeometry_num-replicas.dat` |  state of the rng at which restarted simulation will begin |

Each line in either the scalar or vector estimator files contains a bin which is the average of some measurement over a certain number of Monte Carlo steps.  By averaging bins, one can get the final result along with its uncertainty via the variance.

## General Description

A full understanding of this path integral Monte Carlo code requires an
understanding of the WA-PIMC algorithm alluded to in the introduction.  In this
section, we describe the large building blocks of the code.  The actual
specific details of the implementation can be understood by reading the doxygen
documentation included here as well as reading the actual source code.

Any Monte Carlo simulation whether quantum or classical shares a number of
features in common.  Some type of simulation cell is created with a set of
parameters that describe its physical environment.  The initial state of the
system is guessed, and a series of Moves are performed on the constituents
of the system in such a way that detailed balance is maintained.  After some
suitable equilibration period, measurements are made and their results are
stored to disk.

As discussed above, the driver file for this PIMC program is called pdrive.cpp.
It takes a series of command line options, which are used by the Setup class to
initialize ConstantParameters, Container, LookupTable and Communicator objects.
Next, a Potential object is created which describes the potential environment
(any walls etc.) and the interactions between bosons. A Path object is then
instantiated which holds all the details of the actual world lines of the
quantum particles. An Action object is created based on the Potential which
holds an approximation of the action to be discretized in the path integral
decomposition of the partition function. Finally, the main operating object of
the program, of type PathIntegralMonteCarlo is created, which requires both the
Path and the [Action](@ref ActionBase).  This object performs the actual
simulation via a series of [Moves](@ref MoveBase), all of which generate trial
world line configurations that exactly sample the kinetic part of the density
matrix.  All measurements are made via specific [Estimators](@ref EstimatorBase)
with the results being output to disk.

The main kernel of this program should remain relatively untouched, as it has
been extensively tested and optimized.  Generality can come from modifying just
a few things.  For example, in order to implement a new type of measurement,
one would need to write a derived [Estimator](@ref EstimatorBase) class along
with modifying the Communicator class to define an output path.  New types of
particles and external environments can be added by adding new
[Potential](@ref PotentialBase) then updating Setup to allow for their
specification at the command line.  Finally, radically different systems can be
studied by modifying the [Container](@ref Container) class.
