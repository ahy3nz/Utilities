#!/bin/bash

cd $HOME

echo "Git configure ..."
git config --global user.email "ayang41@gmail.com"
git config --global user.name "Alex Yang"

echo "Installing miniconda..."
curl -o ~/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ~/miniconda.sh -b -p $HOME/miniconda

echo "Creating conda env..."
conda create -n md37 python=3.7 mamba -c conda-forge -y
conda activate md37
mamba install -c conda-forge -c mosdef -c omnia mbuild foyer parmed rdkit openbabel numpy pandas scikit-learn dask matplotlib seaborn networkx -y

echo "Pip installing..."
python -m pip install -e "git+https://github.com/mdtraj/mdtraj.git#egg=mdtraj" --src $HOME/software/
python -m pip install -e "git+https://github.com/mosdef-hub/mbuild.git#egg=mbuild" --src $HOME/software/
python -m pip install -e "git+https://github.com/mosdef-hub/foyer.git#egg=foyer" --src $HOME/software/

echo "Personal libraries..."
cd $HOME/software
git clone https://github.com/ahy3nz/ahy3nz.github.io.git
curl -o ~/.vimrc https://raw.githubusercontent.com/ahy3nz/Utilities/master/.vimrc


