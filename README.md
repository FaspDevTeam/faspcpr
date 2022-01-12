# FASPCPR: README

## Introduction

An efficient parallel Constrained Pressure Residual (CPR) preconditioner with an adaptive "setup phase" (denoted as ASCPR) is developed  for the black oil model. The ASCPR preconditioner is designed based on [Fast Auxiliary Space Preconditioners](http://www.multigrid.org/fasp/)  (FASP) framework, denoted as FASPCPR package. 

FASPCPR is based on the preconditioner method described in the following article:

> Li Zhao and Chunsheng Feng and Chensong Zhang and Shi Shu,
>
> [Parallel Multi-Stage Preconditioners with Adaptive Setup for the Black Oil Model](https://arxiv.org/abs/2201.01970),
>
> 2022, preprint."

## Directory Structure

- data : This folder contains data files 
- include : This folder contains header files
- lib : This folder contains library file
- main : This folder contains main function of the ASCPR method
- src  :  This folder contains sources code 
- util :  This folder contains tools - automatically generate header files

## Build

FASPCPR has the following dependent package:

- FASP,  we recommend version: https://github.com/zhaoli0321/faspsolver.git

  

To build the FASPCPR , first download FASP from the link above and *faspcpr.tar.gz* and *faspsolver.tar.gz* are in the same directory, e.g.: 

```
~> ls
faspcpr.tar.gz  faspsolver.tar.gz
```

Build *FASP*:

```makefile
~> tar -zxvf faspsolver.tar.gz
~> cd faspsolver
~> mkdir Build; cd Build; 
~> cmake -DUSE_OPENMP=ON .. # or cmake -DUSE_PARDISO=ON -DUSE_OPENMP=ON ..
~> make -j 8
~> make install
```

Build *FASPCPR*:

```makefile
~> tar -zxvf faspcpr.tar.gz
~> cd faspcpr
~> make
```

## Running

[**The data**](https://pan.baidu.com/s/1JUHI1y6uSpPjCNRawHOFMw) (Extraction code: fasp) needs to be downloaded and moved into the faspcpr/data directory.

```makefile
export OMP_NUM_THREADS=16  # Specify the number of threads
./test_ascpr.ex
```

## License

This software is free software distributed under the Lesser General Public License or LGPL, version 3.0 or any later versions. This software distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program. If not, see http://www.gnu.org/licenses/.

