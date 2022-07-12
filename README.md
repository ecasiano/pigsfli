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
```
## Dependencies 

The code is written in c++ and makes use of <a href="http://www.boost.org/">boost</a> libraries and, for handling of command-line arguments, the <a href="https://github.com/jarro2783/cxxopts">cxxopts</a> header-only library. For generation of random numbers and the ability to save the state of an RNG for simulation restarts, we use <a href="https://github.com/ajibadeshd/RNG_CLASS">RNG_CLASS</a>.

## Compilation

After dependencies are satisfied, you are now ready to compile the
main pigsfli program on your system.
pigsfli uses CMake for build, test and installation automation. For details on using CMake consult https://cmake.org/documentation/. In short, the following steps should work on UNIX-like systems:

```
cd pigsfli
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
- `-E env CXXFLAGS="xxx"` add additional compiler flags
- `-E env LDFLAGS="xxx"` add additional linker flags

Executables will be installed to `CMAKE_INSTALL_PREFIX` location or if the install is skiped will be located in `build/pigsfli`.
Executables produced are `pigsfli.e`, `pigsflid.e` for `CMAKE_BUILD_TYPE=Release|Debug` respectively.

If you run into problems, failures with linking etc., common errors may include
not properly setting your `LD_LIBRARY_PATH` or not starting from a clean build
directory (issue `make clean` or `rm -rf ./*` inside the build directory).

## Usage

In order to get a quick idea of the options which the code accepts type:
```bash
pigsfli.e --help
```

The code requires various combinations of these options to run, and the help message should give
you an idea about which ones are mandatory.

### Quick Start 

If you want to perform a quick test-run for a small one-dimensional Bose-Hubbard lattice you could try something like:
```bash
./pigsfli.e -D 1 -L 4 -N 4 -l 2 -U 1.995 --mu 1.998 --sweeps 2001 --seed 1968 --measurement-frequency 1 --rng boost_mt19937 --bin-size 10 --bins-wanted 1000 --num-replicas 2 --beta 1.2
```

In order for this to work, you will need a folder named `OUTPUT` in the directory where you type the command as it will produce output files in `OUTPUT` that contain all the results of the code.  Each run of the code is associated with a unique identifying integer: the `PIMCID`.  The options used in this demo include a subset of all the possible options:

| Code Option | Description |
| :-----------: | ----------- |
|`D`     |  Dimension of hypercubic lattice |
|`L`     |  Linear size of hypercube |
|`N`     |  Total number of particles |
|`U`     |  Interaction potential |
|`l`     |  Linear size of hypercubic subregion |
|`sweeps`     |  Number of sweeps before attempting measurements |
|`beta`     |  Set length of imaginary time |
|`mu`     |  Chemical potential |
|`t`     |  Tunneling parameter |
|`canonical`     |  set to false for grand canonical simulation |
|`seed`     |  Random seed value |
|`sweeps-pre` |  Number sweeps for each pre-equilibration step |
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
        _            __ _ _
       (_)          / _| (_)
  _ __  _  __ _ ___| |_| |_
 | '_ \| |/ _` / __|  _| | |
 | |_) | | (_| \__ \ | | | |
 | .__/|_|\__, |___/_| |_|_|
 | |       __/ |
 |_|      |___/

 Path-Integral Ground State (Monte Carlo) For Lattice Implementations


Stage (1/3): Determining mu and eta...

mu: 1.998 eta: 0.5 Z-frac: 1.285%
N     P(N)
3     ***
4     **************************************************************************************************
<N>: 3.97285

mu: -0.984349 eta: 0.190673 Z-frac: 27.9514%
N     P(N)
3     **********************************************
4     *******************************************************
<N>: 3.54578

mu: -1.13738 eta: 0.208524 Z-frac: 21.5446%
N     P(N)
4     *******************************************************************
5     **********************************
<N>: 4.33464

mu: -0.564649 eta: 0.205204 Z-frac: 21.4041%
N     P(N)
3     **********************************
4     *******************************************************************
<N>: 3.66301

mu: -1.12861 eta: 0.204957 Z-frac: 24.5636%
N     P(N)
3     ***********************
4     ***********************************************************
5     ********************
<N>: 3.97396

Fine tuning eta... (Want: 10% < Z-frac < 15%)

mu: -1.12861 eta: 0.206175 Z-frac: 24.4551%
N     P(N)
3     **************************************
4     ******************************************************
5     **********
<N>: 3.72151

mu: -1.12861 eta: 0.103088 Z-frac: 40.3038%
N     P(N)
3     *******
4     **********************************************************************
5     ************************
<N>: 4.17121

Stage (2/3): Equilibrating...

Stage (3/3): Main Monte Carlo loop...

-------- Detailed Balance --------

Insert Worm: 3721/163694
Delete Worm: 3438/3994

Insert Anti: 3525/123099
Delete Anti: 3232/3604

InsertZero Worm: 38086/378100
DeleteZero Worm: 38496/45527

InsertZero Anti: 40086/265865
DeleteZero Anti: 40496/52085

InsertBeta Worm: 37640/375699
DeleteBeta Worm: 37806/45155

InsertBeta Anti: 39576/269448
DeleteBeta Anti: 39742/50686

Advance Head: 96559/97385
Recede  Head: 96107/98277

Advance Tail: 98285/100513
Recede  Tail: 98454/99213

IKBH: 87069/234971
DKBH: 87130/95065

IKAH: 53877/234475
DKAH: 53999/61909

IKBT: 54638/165410
DKBT: 54503/62090

IKAT: 88504/238488
DKAT: 88439/96651

SWAP: 347255/861002
UNSWAP: 347253/934478

SWAP Advance Head: 28600/29543
SWAP Recede Head: 29080/31833

SWAP Advance Tail: 28552/31357
SWAP Recede Tail: 28457/29492

Elapsed time: 1.74691 seconds
```

After this has been completed, you can analyze the results of your run using the scripts in the https://github.com/DelMaestroGroup/papers-code-pigsfli.
