# Documentation  

## Introduction 

This webpage contains the details of a ground state (i.e, zero-temperature) lattice worm algorithm path-integral quantum Monte Carlo (WA-PIMC) code actively developed in c++ since 2020 based on:

- [N.V. Prokof'ev, B.V. Svistunov, I.S. Tupitsyn: Exact, Complete, and Universal Continuous-Time Worldline Monte Carlo Approach to the Statistics of Discrete Quantum Systems](https://arxiv.org/abs/cond-mat/9703200)
- [N. Prokof'ev Lattice Path-Integral Monte Carlo Lecture Notes](http://mcwa.csi.cuny.edu/umass/lectures/part5.pdf)

In it's current version, it can be used to simulate indistinguishable bosons in the one,two, and three-dimensional hypercubic Bose-Hubbard Model. As written, it takes a large number of command line options and allows for the measurement of the system energy and the Rényi Entanglement Entropy between spatial bipartitions of the lattice.

If you have questions on code usage or bug reports, please contact me at ecasiano@vols.utk.edu.

The development and maintenance of this code base has been supported in part by the National Science Foundation under Award No. [DMR-2041995](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2041995&HistoricalAwards=false) and Award No. [DMR-1553991](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1553991)

[<img width="100" src="https://www.nsf.gov/images/logos/NSF_4-Color_bitmap_Logo.png">](http://www.nsf.gov/awardsearch/showAward?AWD_ID=1553991)


## Installation

This program has been successfully compiled and runs on Intel systems using clang and g++. Before installing, one needs to ensure that all dependencies are met.  We recommend that the required libraries (boost) are installed in a `local` folder inside your home directory: `$HOME/local`.

```
git clone https://github.com/ecasiano/pigsfli.git
cd pigsfli
mkdir build
cd build
cmake ../src
make
sudo make install
```

## Dependencies 

The code is written in c++ and makes use of <a href="http://www.boost.org/">boost</a> libraries and, for handling of command-line arguments, the <a href="https://github.com/jarro2783/cxxopts">cxxopts</a> header-only library. For generation of random numbers and the ability to save the state of an RNG for simulation restarts, we use <a href="https://github.com/ajibadeshd/RNG_CLASS">RNG_CLASS</a>.

## Path Integral Monte Carlo

After dependencies are satisfied, you are now ready to compile the
main pigsfli program on your system.
pigsfli uses CMake for build, test and installation automation. For details on using CMake consult https://cmake.org/documentation/. In short, the following steps should work on UNIX-like systems:

  ```
  mkdir build
  cd build
  cmake ../src
  make
  sudo make install
  ```

As above, and with further details below, but you should consider using the following CMake options with the appropriate value instead of xxx :

- `-D NDIM=1|2|3` the number of spatial dimensions
- `-D CMAKE_C_COMPILER=xxx` equal to the name of the C99 Compiler you wish to use (or the environment variable `CC`)
- `-D CMAKE_CXX_COMPILER=xxx` equal to the name of the C++17 compiler you wish to use (or the environment variable `CXX`)
- `-D CMAKE_PREFIX_PATH=xxx` to add a non-standard location for CMake to search for libraries, headers or programs
- `-D CMAKE_INSTALL_PREFIX=xxx` to install pigsfli to a non-standard location
- `-D BOOST_ROOT=xxx` to add non-standard location for Boost install
- `-D STATIC=1` to enable a static build
- `-D CMAKE_BUILD_TYPE=Debug` to build pigsfli in debug mode
- `-D CMAKE_BUILD_TYPE=PIMC` to build pigsfli
- `-D CMAKE_BUILD_TYPE=PIMCDebug` to build pigsfli in debug mode
- `-E env CXXFLAGS="xxx"` add additional compiler flags
- `-E env LDFLAGS="xxx"` add additional linker flags

Executables will be installed to `CMAKE_INSTALL_PREFIX` location or if the install is skiped will be located in `build/pigsfli`.
Executables produced are `pigsl.e`, `pigsld.e`, `pimcl.e`, and `pimcld.e` for `CMAKE_BUILD_TYPE=Release|Debug|PIMC|PIMCDebug` respectively.
^ FIX ME. What executables will be produced now?

If you run into problems, failures with linking etc., common errors may include
not properly setting your `LD_LIBRARY_PATH` or not starting from a clean build
directory (issue `make clean` or `rm -rf ./*` inside the build directory).

## Usage

FIXME!!!!!

NOTE: NEED TO ADD EXAMPLE OF COMPILATION

In order to get a quick idea of the options which the code accepts type:
```bash
pigsfli.e --help
```

The code requires various combinations of these options to run, and the help message should give
you an idea about which ones are mandatory.

### Quick Start 

If you want to perform a quick test-run for a small one-dimensional Bose-Hubbard lattice you could try something like:
```bash
./pigsfli.e -D 1 -L 4 -N 4 -l 2 -U 100.0 --mu 6.0 --subgeometry square --measurement-frequency 1 --beta 2.0 --rng boost_mt19937 --num-replicas 1 --seed 17 --bin-size 2 --bins-wanted 5
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
|`sweeps-pre` |  sweeps-pre |
|`bin-size`     |  Number of measurements per bin|
|`bins-wanted`     |  Number of bins desired in data file|
|`subgeometry`     |  Shape of subregion: square OR strip|
|`num-replicas`     |  Number of replicas|
|`measurement-frequency`     |  Measurements will be performed every other this amount|
|`rng`     |  Random Number Generator type|
|`restart`     |  continue simulation from a loaded rng state|
|`no-accessible`     |  do not calculate accessible entanglement entropies|

All options, including lists of possible values and default values can be seen
by using the `--help flag`.

The output of the above command should yield something like:
```bash
sub-sites: 0 1
U: 100

  _           _   _   _           ______ _____ _____  _____
 | |         | | | | (_)          | ___ \_   _|  __ \/  ___|
 | |     __ _| |_| |_ _  ___ ___  | |_/ / | | | |  \/\ `--.
 | |    / _` | __| __| |/ __/ _ \ |  __/  | | | | __  `--. \
 | |___| (_| | |_| |_| | (_|  __/ | |    _| |_| |_\ \/\__/ /
 \_____/\__,_|\__|\__|_|\___\___| \_|    \___/ \____/\____/


Stage (1/3): Determining mu and eta...

mu: 6 eta: 0.5 Z-frac: 35.9215%
N     P(N)
4     *****************************************************************************************************
<N>: 4

mu: 6 eta: 0.373263 Z-frac: 45.8905%
N     P(N)
3     *
4     ****************************************************************************************************
<N>: 3.99975

Fine tuning eta... (Want: 10% < Z-frac < 15%)

mu: 6 eta: 0.376703 Z-frac: 45.1917%
N     P(N)
3     *
4     ****************************************************************************************************
<N>: 3.99886

mu: 6 eta: 0.546219 Z-frac: 32.8821%
N     P(N)
3     *
4     ****************************************************************************************************
<N>: 3.99916

mu: 6 eta: 0.273109 Z-frac: 55.6037%
N     P(N)
4     *****************************************************************************************************
<N>: 4

mu: 6 eta: 0.396009 Z-frac: 43.3001%
N     P(N)
3     *
4     ****************************************************************************************************
<N>: 3.99993

Stage (2/3): Equilibrating...

Stage (3/3): Main Monte Carlo loop...

-------- Detailed Balance --------

Insert Worm: 29800/1160621
Delete Worm: 29469/116950

Insert Anti: 301805/1084255
Delete Anti: 288628/418619

InsertZero Worm: 100333/2156306
DeleteZero Worm: 106928/1447180

InsertZero Anti: 405936/2101899
DeleteZero Anti: 412531/885665

InsertBeta Worm: 100714/2162849
DeleteBeta Worm: 107627/1444361

InsertBeta Anti: 404204/2092847
DeleteBeta Anti: 411117/887977

Advance Head: 754399/754399
Recede  Head: 754909/755013

Advance Tail: 751144/751287
Recede  Tail: 750145/750145

IKBH: 51408/2053874
DKBH: 51350/217882

IKAH: 558472/2049219
DKAH: 558245/593583

IKBT: 555257/1997134
DKBT: 555590/590990

IKAT: 51200/2043258
DKAT: 51152/215920

SWAP: 0/0
UNSWAP: 0/0

SWAP Advance Head: 0/0
SWAP Recede Head: 0/0

SWAP Advance Tail: 0/0
SWAP Recede Tail: 0/0

beta: 2

sweeps: 1e+07
Z_ctr: 10
Z_frac: 83.3333% (10/12)

<N>: 4

Elapsed time: 6.75028 seconds
```

After this has been completed, you can analyze the results of your run using the scripts in the https://github.com/DelMaestroGroup/papers-code-pigsfli.
