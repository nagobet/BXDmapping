# Momix

# Date: 14 July 2021

# Goal: running momix jupyter notebook remotely

# follow instructions from https://github.com/ComputationalSystemsBiology/momix-notebook
# download notebooks from https://zenodo.org/record/4194102

# install conda from computer 3009
# create a new environment:
conda create -n momix -c conda-forge -c bioconda -c lcantini momix r-irkernel
# activate the environment
conda activate momix
# Launch the notebook with jupyter-notebook.

# then follow from https://docs.anaconda.com/anaconda/user-guide/tasks/remote-jupyter-notebook/
jupyter notebook --no-browser --port=8080

# from git bash on personal laptop
ssh -L 8080:localhost:8080 ngobet@pccig3009.unil.ch
# connect with 3009 identification

# copy the link from 3009 terminal and open in browser
http://localhost:8080/tree
# naviguate through the tree to find jupyter notebook of choice
