LS_COLORS=$LS_COLORS:'di=1;36:fi=1;35:ex=1;37:'; export LS_COLORS
alias ls='ls -h'
export WHITBY='yangah@whitby.vuse.vanderbilt.edu'
export FILEY='ahy3nz@filey.vuse.vanderbilt.edu'
export CORI='ahy3nz@cori.nersc.gov'
export EDISON='ahy3nz@edison.nersc.gov'
export TITAN='ahy3nz@titan.ccs.ornl.gov'
export ACCRE='yangah@login.accre.vanderbilt.edu'
export RAHMAN='ahy3nz@rahman.vuse.vanderbilt.edu'
export PYTHONPATH="/users/ahy3nz/Programs/groupy2_alex/groupy2":$PYTHONPATH
export PYTHONPATH="/users/ahy3nz/Programs/Analysis/permeability":$PYTHONPATH
alias vmd='/Applications/VMD\ 1.9.1.app/Contents/MacOS/startup.command'
export PATH="/usr/local/gromacs/bin/GMXRC":$PATH
export PATH="/usr/local/gromacs/bin/gmx":$PATH
export PATH="/usr/local/gromacs/bin/":$PATH
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
#export PATH="/raid6/homes/ahy3nz/anaconda3/bin:$PATH"
