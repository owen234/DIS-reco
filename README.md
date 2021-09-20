# DIS-reco




## What's in here?

The stuff in this repository is for looking at the reconstruction of Q2, x, and y, using both the scattered electron and the inclusive
hadronic final state.


## Machine Learning

There are lots of Jupyter Notebooks in here for both training DNNs and also making plots from a completed DNN training.
If you haven't used Jupyter Notebooks before, it's pretty easy to learn.  In a nutshell, you start a server on your machine
from a terminal that you then connect to with your web browser (I use Safari on a Mac).  You enter code in the web browser.
When you are ready to run it, you hit shift-enter within that block and that signals the server to run that code and then
present the results in the browser.  There are a lot of tutorials out there that will help you learn the basics of
Jupyter Notebooks.

If you get errors, it may be that you need to install some more python packages.  I use "pip install" for everything
and it works fine.  

If you need the input root files, ask me.

## Fast simulation (Delphes).


The section below will get you started with a setup for running the fast simulation (Delphes).  Once you have that, you can
make a simple reduced TTree with the root c++ code here and then analyze it with the jupyter notebook.  Note that you don't
have to do all of the below in the EIC virtual machine you use to run the full simulation.  You should be able to install
Pythia and Delphes directly on your laptop.  It worked for me on my Mac.

**0.  Setup Root environment in shell**

Go to your root installation location and do
```
cd <your root dir>
cd bin
source thisroot.sh
```

**1. Install Pythia**

Instructions are here: https://pythia.org/releases/

Here's what I did.

```
cd <where you want to put all this stuff>
mkdir pythia
cd pythia

<Download the tar file for the most recent release.>

tar -xvf pythia8306.tgz
cd pythia8306
mkdir install

<change the absolute path in the next command to point to your directory.  If you use bash, do equivalent export.>

setenv PYTHIA8 /Users/owen/work/eic/pythia/pythia8306/install

./configure --prefix=$PYTHIA8 | & tee configure.log
make install | & tee build.log

setenv LD_LIBRARY_PATH ${PYTHIA8}/lib:$LD_LIBRARY_PATH

make HAS_PYTHIA8=true

```


**2. Download and build Delphes**

Here's the link:  https://cp3.irmp.ucl.ac.be/projects/delphes  (click on Download for the tar file)

and here's what I did

```
mkdir delphes
cd delphes

<download tar file into this dir>

tar -xvf Delphes-3.5.0.tar.gz
cd Delphes-3.5.0
make HAS_PYTHIA8=true | & tee build.log

```

**3. Check out the ATHENA / EIC implementation in Delphes**

Here's the documentation :  https://delphes-eic.readthedocs.io/en/latest/

and here's the github repository :  https://github.com/eic/delphes_EIC

```
cd <your delphes directory, the one below Delphes-3.5.0 where Delphes-3.5.0.tar.gz is>

git clone https://github.com/eic/delphes_EIC.git

```



**4.  Try running Pythia + Delphes**

```
cd <your delphes directory, the one below Delphes-3.5.0 where Delphes-3.5.0.tar.gz is>

cd Delphes-3.5.0

./DelphesPythia8 ../delphes_EIC/delphes_card_allsilicon_3T.tcl ../delphes_EIC/pythia8cards/NC_DIS.cmnd out.root

```

That will produce the file out.root, which contains a big TTree with the gen-level particles and the fast simulation
of the detector reconstruction.


**5.  Make the minitree from the Delphes output**

```
<start an interactive root session>
.L minitree_maker.c
minitree_maker m("<your delphes root file here>")
m.Loop(true,10)
```
That will print out a lot of infor for the first 10 events.  To run over the whole file, do this
```
m.Loop(0,-1)
```
The output file is mini-tree.root











