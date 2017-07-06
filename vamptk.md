# ampTK install scripts
Reworking Jon Palmer's [installation instructions](https://github.com/nextgenusfs/amptk/blob/master/docs/ubuntu_install.md) so that they work for me. If you're unclear of which Python version you have, I'd suggest looking at installing either Conda or Miniconda as it will allow you to manage all of your installs as well as ensure you have an updated Python version which should cause less headaches. To do that, [click here](https://www.continuum.io/content/conda-data-science) and download either the full Conda installation or the smaller Miniconda (for our purposes, Miniconda will do just fine). As an example using Miniconda, [download your version](https://conda.io/miniconda.html) either using the browser and subsequent install, or stick to the command line following [these instructions](https://conda.io/docs/install/quick.html#linux-miniconda-install), condensed below:
```
## Downloading to $HOME directory, but you could just as well put it in '/usr/local'...
cd $HOME
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
If you're unlucky like me and don't have super user priviliges and can't manually just upgrade the python version to what you want, you can always specify to use this shiny new version of Python first, over whatever default alias is attributed to your $PATH when you type `python` or `python3` by logging into your `.bashrc` script and adding the following lines:
```
#default path for Python3 changed to that set in Miniconda install:
alias python3=/home/fosterLab/devonr/miniconda3/bin/python3
```
You'll need to exit the current session and open up a new shell, but that should ensure that every time you type `python3` it points to your updated Miniconda version. Exit the terminal for the changes to take place, and open a new shell (window).  

Why bother doing this at all? Because Conda's software includes management of packages outside of Python (whereas Pip sticks to the Python universe only). This is a huge help when keeping everything together in one place, regardless of whether it's a Python package or not. There are a few things to do once you have Conda (I'm going to call it Conda whether you have Conda or Miniconda installed - the package manger you installed was Conda either way) installed - check out using bioconda (which we'll use a ton) [here](https://bioconda.github.io/). These steps are condensed from that link:
```
(conda config --add channels r)
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```
Now we can use a `conda install ...` command to quickly install most of our programs. More on that next.  


## virtual enviornment space  and install all packages there
See [here](https://virtualenv.pypa.io/en/stable/installation/) for information about installing the software needed to create virtual environments. See [this documentation](http://python-guide-pt-br.readthedocs.io/en/latest/dev/virtualenvs/) for a brief overview of actually creating virtual environments. We're creating a sub-directory named `venv` within $HOME, and naming this environment `vamptk`. All programs will then be installed within `vamptk` directories. Note that we do not need to store data here - it just serves as the program repository so that you can update all your program information as needed (without moving your raw or filtered data around). The steps are as follows:  
- Starting in the $HOME directory, make a new folder named **venvs**:
```
cd $HOME
mkdir venvs
```
- We're going to move into that directory and create a new virtual environment, named **vamptk**. You can make a very unassuming mistake here - depending on where your `virtualenv` command points to, you might be installing a virtual environment using Python2 or Python3; furthermore you might be installing a version of Python2 or 3 that you didn't want if you happen to have multiple versions installed. Get some help using Python3 virtual environments [here](https://docs.python.org/3/library/venv.html), and check out  I'm going to demonstrate two different ways that you might get something to work. *Note that Flavor1 is what we want, while Flavor2 is NOT what you want to do*  

```
## Assumes Python2 and Python3 both installed
## Move into the new directory you just created
cd venvs

## Now create a new virtual environment. 
## Flavor1 - you specify to build a Python3 virtual environment:
python3 -m venv vamptk

## Flavor 2
virtualenv vamptk
## this creates a virtual environment with Python2.7 as your default language. This isn't what you want in our case.
```
It can be worse than that too. If you are using an older version of Python3, anything before Python 3.5 (so 3.4, 3.3, etc...) then the `venv` command won't work, as they changed the sytnax to execute that command between versions 3.4 and 3.5. If you're stuck with v3.4 or less, you're going to use a few different commands. There are two things you'll do: first, you're going to use a slightly different installation sytnax to get the virtual environment running on the older Python3 variety - see [this site](https://robinwinslow.uk/2013/12/26/python-3-4-virtual-environment/#pyvenv-3-4) for a helpful example of how you'd do that in general. However, you're potentially going to run into an error because of a bug between the Linux OS and Python version you're trying to put together - see [this message](https://askubuntu.com/questions/488529/pyvenv-3-4-error-returned-non-zero-exit-status-1). Fortunately for you, I've condensed that info below into something that works pretty easily:
```
## Which Python3 version are you using?
[my@machine]: python3 --version
Python 3.4.0

## Because I have version 3.4, I'm going to have to use the 'pyvenv' command to install my virtual environment, rather than the 'venv' command as noted above.
pyvenv-3.4 --without-pip vamptk
source vamptk/bin/activate
curl https://bootstrap.pypa.io/get-pip.py | python
deactivate
```
Regardless of what Python3 version you're using, the previous command should create a new folder within `$HOME/venvs/` called **vamtpk**. This parent directory will serve as the repository for all the program information we need moving forward. Each time you want to use this repository (and therefore programs), you need to activate the virtual environment:  
```
source activate /vamptk/bin/activate

#full path would also work...
[my@machine]: source activate $HOME/vamptk/bin/activate
```
You should notice that the name of your prompt will have changed to indicate that the virtual environment is now your working environment. So if before you looked like `[my@machine]:`, you after sourcing the virtual environment named **vamptk** it should look like `(vamptk)[my@machine]:`.  
Likewise, to get out of the virtual environment, you can always just type `deactivate` and you should see your prompt on the left side disappear from `(vamptk)[my@machine]` to `[my@machine]`.


## installation of dependencies
Once you've installed your virtual environment it's time to install a slew of dependencies. Many of these are easy to install by typing just a one line command, though some more of these require a bit more work. These packages enable Jon's Python scripts to run, but are certainly capable of completing many more operations than what are described within the **amptk** scripts. For some documentation on many of these packages, [start here](https://www.scipy.org/index.html).  

Let's do the easy stuff first, installing python modules with *pip*. Ensure that you are working within the *vamptk* directory and have that virtual environment active (look at the prompt - do you see something like ```vamptk[your@terminal]```?). By the way, these are pretty big programs so go get a cup of coffee and come back in about 20 minutes.
```
conda install biopython natsort pandas numpy matplotlib biom-format psutil
```  
We're going to install a few other programs not listed in Jon's recommendations which can be very helpful for later applications such as plotting and note taking, and we'll use the same strategy as above with *pip*:  
```
conda install ipython jupyter sympy
```  
Notably, we could have substitute `conda` with `pip` and been successful. However, keeping everything in the `conda` family keeps all the packages in one place. Up to this point we could use `conda` or `pip` because they are all Python packages. However, as we start adding more programs to the mix which fall outside of the Python universe, we'd need to start working harder to incorporate all these non-Python programs into the same directories (or we couldn't do that at all). Because I don't want to work hard, we started with `conda` deliberately to avoid these headaches.  Now on to the slightly more challenging installs.    

**USEARCH**
This is the most peculiar and yet simple of the bunch. It's the only program you can't install from the command line directly, so you can't use `conda install` (or `pip install` for that matter). It also have a bunch of possible versions you might download, and some of these work well with **amptk** and others do not (it's not **amptk** it doesn't work with, it's all the denendencies **amptk** relies on that cause the commotion). At the time of this writing (7.6.17) the USEARCH version which was playing nicely with **amptk** was v9.2.64. I'd recommend downloading that one for now.  
- First, go to the [program website](http://www.drive5.com/usearch/download.html) and enter your email address. Then go to your email account and click the link to download the program binary file, which will likely download them to your personal computer's `~/Downloads/` directory, and likely have a name along the lines of `usearch9.2.64_i86linux32`.
- And now we're back to the command line.  Rename this file `usearch9`. 
```
## assumes you downloaded file from website to your ~/Downloads directory on your local computer
[your@localmachine]: cd ~/Downloads
[your@localmachine]: mv usearch9.2.64_i86linux32 usearch9
```
Now here's where things get different for everyone - I'm installing all these programs on a remote machine. You could be doing this directly on your own local computer. If you're like me and need to move this file to your remote machine, then you need to copy this file to the remote machine, otherwise, just use a `mv` command to transfer the file from your ~/Desktop location to your virtual environment location (the path you're placing the file in is the same regardless):
```
## for a local transfer:
[your@localmachine]: mv usearch9 $HOME/venvs/vamptk/bin/.

## using rsync to move file from local to remote machine (could use scp, etc.) requires you to use the full path for where you are placing files
## my full path to $HOME is equivalent to /home/fosterlab/devonr
## you can find your equivalent path by typing: $HOME
## when you type your own command, replace the "$HOME" variable with whatever your equivalent is

[your@localmachine]: rsync usearch9 your@remotemachine:/$HOME/venvs/vamptk/bin/
```
That's it! **USEARCH** is fully executable binary file which you don't have to do anything with (and can't!). You should be able to call it by typing `usearch9` and see a very basic help menu print out on the screen.  


**Command line installs of other dependencies**
Back to our Conda world. Pretty simple for many of these programs with just a few lines - see [this link](https://www.continuum.io/blog/developer/jupyter-and-conda-r) for a brief description of the 'r-essentials' package we install in the second command:  
```
## Install non-R dependencies
conda install bedtools vsearch

## Now install R (this includes R and some helpful packages)
conda install -c r r-essentials

## Install a special R package not installed with that previous command
conda install bioconductor-dada2
```

Last thing is to actually install the **amptk** program itself, which we, alas, can't do in Conda. We're going to use **git** to install the final package. See [this page](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) if you don't have git already installed. It's another simple one-liner with git installed, but because we're not using the typical `conda install` command, we have to be specific about where we want to put this file. We're going to move to the specific directory where `conda install ...` packages have been placed, then make a linke for the specific binary files to where they should be too:
```
## change directories to where packages are installed for Miniconda:
cd $HOME/miniconda3/pkgs

## install amptk program
git clone https://github.com/nextgenusfs/amptk.git
pwd
[here you should see the full path to the directory]... for example: /home/users/yourname/miniconda3/pkgs/amptk

## create a soft link for the path to that binary in another directory
cd $HOME/miniconda3/bin
ln -s /home/fosterLab/devonr/miniconda3//pkgs/amptk/amptk
```

If you happen to get a weird message, try this:
```
conda remove --force readline
```
