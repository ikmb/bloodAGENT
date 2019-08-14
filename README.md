# deep blood group typing of long read sequencng data 
There are some dependencies that I solve in a special way. So I do not implement everything by submodules but clone repositories next to each other. 

## Installing the tool at the ikmbhead.uni-kiel.de

Setup environment first

`module load IKMB`

`module load GCC/6.1`

To install this pipeline, simply clone the repository to a location:

`git clone git@git.ikmb.uni-kiel.de:m.wittig/BfxCppClasses.git MyTools`

`cd MyTools`

`make all`

`cd ..`

MyTools libs created

`git clone https://github.com/samtools/htslib.git`

`cd htslib`

`make`

`cd ..`

htslib created

`git clone https://github.com/dpryan79/libBigWig.git`

`cd libBigWig`

`make`

`cd ..`

libBigWig created

`git clone --recursive git@git.ikmb.uni-kiel.de:m.wittig/deepBlood.git`

`cd deepBlood`

`make all`

`cd ..`

ngsblood created


# setup environment and run

`export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<PATH to libhts.so.2> #/ifs/data/nfs_share/sukko545/software/htslib`

`./dist/Release/GNU-Linux/deepBlood`

run blood typing example:


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
Components of this project are licensed under different licenses



