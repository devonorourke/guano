# ampTK install scripts
Reworking Jon Palmer's [installation instructions](https://github.com/nextgenusfs/amptk/blob/master/docs/ubuntu_install.md) so that they work for me.  

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
It can be worse than that too. If you are using an older version of Python3, anything before Python 3.5 (so 3.4, 3.3, etc...) then the `venv` command won't work, as they changed the sytnax to execute that command between versions 3.4 and 3.5. If you're stuck with v3.4 or less, you're going to use a few different commands. This is what happens when you're stuck with old versions, and different updates of different programs get worse over time. There are two things you'll do: first, you're going to use a slightly different installation sytnax to get the virtual environment running on the older Python3 variety - see [this site](https://robinwinslow.uk/2013/12/26/python-3-4-virtual-environment/#pyvenv-3-4) for a helpful example of how you'd do that in general. However, you're potentially going to run into an error because of a bug between the Linux OS and Python version you're trying to put together - see [this message](https://askubuntu.com/questions/488529/pyvenv-3-4-error-returned-non-zero-exit-status-1). Fortunately for you, I've condensed that info below into something that works pretty easily:
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
pip install biopython natsort pandas numpy matplotlib biom-format psutil
```  
We're going to install a few other programs not listed in Jon's recommendations which can be very helpful for later applications such as plotting and note taking, and we'll use the same strategy as above with *pip*:  
```
pip install ipython jupyter sympy
```  

Now on to the slightly more challenging installs:
**USEARCH**
This is the most peculiar and yet simple of the bunch. It's the only program you can't install from the command line directly. It also have a bunch of possible versions you might download, and some of these work well with **amptk** and others do not (it's not **amptk** it doesn't work with, it's all the denendencies **amptk** relies on that cause the commotion). At the time of this writing (7.6.17) the USEARCH version which was playing nicely with **amptk** was v9.2.64. I'd recommend downloading that one for now.  
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



