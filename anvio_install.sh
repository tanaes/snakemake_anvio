# Create conda env for anvio2

conda create -n anvio2 python=2 numpy h5py scipy scikit-learn six requests 

source activate anvio2

condaenv=~/miniconda/envs/anvio2
condabin=${condaenv}/bin

# ETE2
conda install -c etetoolkit ete2

# Prodigal
conda install -c biocore prodigal

pip install bottle

# Cython
pip install Cython


cd ${condabin}

# install HMMER
wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz
tar -xzvf hmmer-3.1b2-linux-intel-x86_64.tar.gz
ln -s hmmer-3.1b2-linux-intel-x86_64/binaries/* .


export C_INCLUDE_PATH=$C_INCLUDE_PATH:${condaenv}/include/
export LIBRARY_PATH=$LIBRARY_PATH:${condaenv}/lib/

# GSL
cd ${condaenv}
wget ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
tar -xzvf gsl-latest.tar.gz
cd gsl-2.1/
./configure --prefix=${condaenv}
make
make install


# Centrifuge
export CENTRIFUGE_BASE=${condaenv}/centrifuge
mkdir $CENTRIFUGE_BASE
cd $CENTRIFUGE_BASE
git clone https://github.com/infphilo/centrifuge
cd centrifuge && make

cd $CENTRIFUGE_BASE
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/b+h+v.tar.gz
tar -zxvf b+h+v.tar.gz && rm -rf b+h+v.tar.gz


# install anvio
pip install anvio

source deactivate


# Set paths for Conda env
mkdir -p ${condaenv}/etc/conda/activate.d
mkdir -p ${condaenv}/etc/conda/deactivate.d
touch ${condaenv}/etc/conda/activate.d/env_vars.sh
touch ${condaenv}/etc/conda/deactivate.d/env_vars.sh

# store old env vars
echo '#!/bin/sh/' >> ${condaenv}/etc/conda/activate.d/env_vars.sh
echo "export OLD_PATH=\$PATH" >> ${condaenv}/etc/conda/activate.d/env_vars.sh

# set new env vars
echo "export CENTRIFUGE_BASE=${CENTRIFUGE_BASE}" >> ${condaenv}/etc/conda/activate.d/env_vars.sh
echo "export PATH=\$PATH:${CENTRIFUGE_BASE}/centrifuge" >> ${condaenv}/etc/conda/activate.d/env_vars.sh

# unset env vars
echo '#!/bin/sh/' >> ${condaenv}/etc/conda/deactivate.d/env_vars.sh
echo "unset CENTRIFUGE_BASE" >> ${condaenv}/etc/conda/deactivate.d/env_vars.sh
echo "export PATH=\$OLD_PATH" >> ${condaenv}/etc/conda/deactivate.d/env_vars.sh
