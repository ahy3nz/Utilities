LS_COLORS=$LS_COLORS:'di=1;36:fi=1;35:ex=1;37:'; export LS_COLORS
alias ls='ls -h -G'
export WHITBY='yangah@whitby.vuse.vanderbilt.edu'
export PERCUS='ahy3nz@percus.vuse.vanderbilt.edu'
export FILEY='ahy3nz@filey.vuse.vanderbilt.edu'
export CORI='ahy3nz@cori.nersc.gov'
export EDISON='ahy3nz@edison.nersc.gov'
export TITAN='ahy3nz@titan.ccs.ornl.gov'
export RAHMAN='ahy3nz@rahman.vuse.vanderbilt.edu'
export FILEY='ahy3nz@filey.vuse.vanderbilt.edu'
export ACCRE="yangah@login.accre.vanderbilt.edu"
#export PATH="/raid6/homes/ahy3nz/Programs/packmol:$PATH"
export PYTHONPATH="/raid6/homes/ahy3nz/Programs/":$PYTHONPATH
export PYTHONPATH="/raid6/homes/ahy3nz/Programs/Utilities":$PYTHONPATH
export PYTHONPATH="/raid6/homes/ahy3nz/Programs/groupy2_alex/groupy2":$PYTHONPATH
export PYTHONPATH="/raid6/homes/ahy3nz/Programs/permeability/permeability":$PYTHONPATH
export PYTHONPATH="/raid6/homes/ahy3nz/Programs/Analysis/Bilayers":$PYTHONPATH
export PYTHONPATH="/raid6/homes/ahy3nz/Programs/cg_mapping":$PYTHONPATH
#export PYTHONPATH="/raid6/homes/ahy3nz/Programs/setup/Bilayer/":$PYTHONPATH
#export PYTHONPATH="/raid6/software/hoomd/hoomd_1.3.3/lib/hoomd/python-module/":$PYTHONPATH
#export PYTHONPATH="/raid6/homes/ahy3nz/Programs/setup/CoarseGraining/":$PYTHONPATH
#export PYTHONPATH="/raid6/homes/ahy3nz/anaconda3/envs/hack35/lib/python3.5/site-packages":$PYTHONPATH
# Default hoomd
#export PYTHONPATH="/raid6/software/hoomd/hoomd_2.1.2/":$PYTHONPATH

# my hoomd
#export PYTHONPATH="/raid6/homes/ahy3nz/Programs/hoomd-blue/install/":$PYTHONPATH
# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=

# User specific aliases and functions
unset SSH_ASKPASS

# added by Anaconda3 4.2.0 installer
export PATH="/raid6/homes/ahy3nz/anaconda3/bin:$PATH"
