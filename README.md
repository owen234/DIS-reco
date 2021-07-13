# DIS-reco

This will get you started with a setup for running the fast simulation (Delphes).

**1. Install Pythia**

Instructions are here: https://pythia.org/releases/

Here's what I did.

```
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


