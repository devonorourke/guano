# Installation steps

All information for `amptk` installation is [documented online](http://amptk.readthedocs.io/en/latest/). The following information reflects changes to what was already used in previous `amptk` usage from other projects.

## amptk install
Install amptk in `~/pkgs` directory, rename directory, then and add symbolic link in `~/bin`:

```
cd $HOME/pkgs
git clone https://github.com/nextgenusfs/amptk.git
mv amptk amptk-1.0.3
cd $HOME/bin
ln -s /mnt/lustre/macmaneslab/devon/pkgs/amptk-1.0.3/amptk .
```

## other dependencies
Updated [vsearch](https://github.com/torognes/vsearch) by adding newest version; updates symbolic link in `~/bin` (kept 2.4.2 though in case of program crash)
```
cd $HOME/pkgs
wget https://github.com/torognes/vsearch/releases/download/v2.6.2/vsearch-2.6.2-linux-x86_64.tar.gz
tar xzf vsearch-2.6.2-linux-x86_64.tar.gz
cd $HOME/bin
ln -s /mnt/lustre/macmaneslab/devon/pkgs/vsearch-2.6.2-linux-x86_64/bin/vsearch .
```

> usearch 9.2.64 was previously installed and not required to be updated; no usearch-10 was installed as this is not going to be used in our work

Next, create a virtual environment using `conda`; proceed with `y`.
```
module purge
module load anaconda/colsa
conda create -n amptk_env python=2.7
```

Load the venv, then install remaining `amptk` dependencies. These include both python modules as well as R libraries.
```
source activate amptk_env
conda install biopython natsort pandas numpy matplotlib seaborn edlib biom-format psutil
conda config --add channels r
conda install --yes r-base bioconductor-dada2
conda install --yes r-base bioconductor-phyloseq
```
