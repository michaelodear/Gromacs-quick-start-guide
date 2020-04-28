# Gromacs quick start guide

1. [Setting up Gromacs](#setup) <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.1 [Install with GPU (recommended)](#gpu) <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.2 [Install with CPU](#cpu)
2. [Preprocessing (Optional)](#pre)
3. [Running your first Molecular Dynamics Simulation](#md1)
4. [Visualising your simulation](#vis)
5. [Secondary structure analysis with DSSP](#dssp)

When I was first getting started with GROMACS I was frustrated by how much of the documentation and tutorials were outdated. One notable exception to this are [Justin Lemkul's tutorials](http://www.mdtutorials.com/gmx/index.html), which are an excellent resource for running MD simulations that I highly recommend everyone check out, and which the 'Running your first MD section' borrows heavily from. 

This guide's aim is to collate a number of resources from around the web to walk a new user through the process of running their first molecular dynamics simulation, from installation to analysis. It is aimed at people without any programming background and written assuming that you are working on a Linux machine running Ubuntu 18.04. If you don't have access to a linux machine, this [very short guide from Microsoft](https://docs.microsoft.com/en-us/windows/wsl/install-win10) explains how to run Linux on a Windows machine.

### A Quick note on the Terminal/Command Line programming

Given that this guide is written for the vast majority of biologists without a programming background, I know most of you won't have any experience with how to use the command line and might find it a little intimidating. As such, with this tutorial I've tried to write it so that you don't need to know how to use the command line at all. However, just knowing the very basics of how to navigate directories and organise files will help you hugely when analysing your results. This is a [pretty good cheatsheet](https://www.linuxtrainingacademy.com/linux-commands-cheat-sheet/) that I recommend referencing when you need to.

<a name="setup"></a>
## Setting up Gromacs

<a name="gpu"></a>
### Install with GPU (recommended)

This is the default way to set-up GROMACS and will result in much faster performance than running with a CPU.

#### Check CUDA installation

First, before we can get GROMACS, we need to check whether the machine you are working on has a CUDA installed. You don't need to know what this is, just that it is necessary for doing most computing with a GPU these days. To check if you already have a CUDA, run this:

`nvcc --version`

If you see any output that is not

```
Command 'nvcc' not found, but can be installed with:

sudo apt install nvidia-cuda-toolkit
```

then you should be fine. If you do see this message, then the install instructions for CUDA on linux are as follows:

```
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/cuda-ubuntu1804.pin
sudo mv cuda-ubuntu1804.pin /etc/apt/preferences.d/cuda-repository-pin-600
sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/7fa2af80.pub
sudo add-apt-repository "deb http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/ /"
sudo apt-get update
sudo apt-get -y install cuda
```

This will take a little while to run. But once it's finished, we can move on to:

#### Downloading and installing GROMACS (GPU enabled)

To start, you'll want to download the GROMACS package for your desired version. We'll be using 2018.1 in this guide, which can be found [here](http://manual.gromacs.org/documentation/2018.1/download.html), however, other versions are available [here](http://manual.gromacs.org/documentation/). You'll want to select the http option. Then, run the following commands in the folder containing the tar file you just downloaded. This will set you up to be able to run gromacs from anywhere in the terminal simply by typing 'gmx' followed by your desired command. If you want to install a different version, simply replace 2018.1 with your desired version name. 

```
tar xfz gromacs-2018.1.tar.gz
cd gromacs-2018.1
mkdir build
cd build
cmake .. -DGMX_GPU=ON -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda
make
make check
sudo make install
source /usr/local/gromacs/bin/GMXRC
```

When that's done, to check that GROMACS has successfully installed run:

`gmx --version`

If everything has worked properly, then you should see a list of details about your GROMACS installation. Make sure to check that 'GPU = enabled' before proceeding to simulation.

<a name="cpu"></a>
### Install with CPU

If you don't have a GPU available, the quickest way to get GROMACS running on Ubuntu is to use the native package installer, sudo apt-get:

`sudo apt-get install gromacs`

This will install the 2018.1 version of GROMACS, which is the one used in this guide. 'Sudo' will prompt you for a password, which will be the login of the admin user.

<a name="pre"></a>
## Preprocessing (Optional)

##### Note: To start the MD your protein will need to be in a .pdb format. If it's already known to science, you can probably find it by searching in the [Protein Data Bank](https://www.rcsb.org/). If your specific protein is not documented but it shares significant sequence homology with other known proteins, then you can try using [SWISSMODEL](https://swissmodel.expasy.org/) to generate a pdb file from the homologous reference. If neither of those are an option, then [I-Tasser](https://zhanglab.ccmb.med.umich.edu/I-TASSER/) is a useful webserver for generating protein structures from scratch, albeit with a bit of a waiting time If you want to use the same protein that this tutorial will be working with, you can download it [here](https://www.rcsb.org/structure/1fxk)

In this tutorial we will be looking at the effect ofhosphorylation on or protein's structure. To do this, we'll need to generate a copy of our proteins with the desired post-translational modifications. A brilliant way to do this is the [Vienna-PTM webserver](http://vienna-ptm.univie.ac.at/), a tool that is simple, quick and intuitive to use. Simply upload your protein, select the 'force field' you want (We'll be using ffG54a7) and then on the next page select the residues you want to change and how to change them. To use the modified protein, you'll need to use a slightly different version of this force field, which can be found [here](http://vienna-ptm.univie.ac.at/?page_id=100) (you want the one for GROMACS 4.5x and up). Unzip that in the same directory as your .pdb file. If you're new to GROMACS, force fields are essentially a series of equations which define how the molecules in your simulation should interact. We're using ffG54a7 because of it's compatibility with Vienna-PTM, however, if your project involves different preprocessing, or doesn't require it at all I recommend researching some of the other force fields such as AMBER to find what is right for your particular problem.

While this tutorial uses phosphorylation of proteins as it's example, GROMACS is useful for simulating many other problems. For instance, if you were interested in investigating alternative splicing, this is where you could generate the different isoforms of the protein you wanted to look at.

<a name="md1"></a>
## Running your first MD

Alright, now that our proteins are prepared and in the right format we can start our first simulation! As previously mentioned, this section leans heavily on Justin Lemkul's brilliant work and I highly recommend checking out his [tutorials](http://www.mdtutorials.com/gmx/) if you want any more details on the following steps.

### Generating an inital topology of your protein

To start off, you'll want to convert your protein to a format GROMACS better understands. But first, we'll need to remove crystal water from our protein. This isn't necessary for every protein, but is in our case. To do that , just run the following:

`grep -v HOH 1fxk.pdb > 1fxk_clean.pdb`

Great! You're molecule is now water free, at least temporarily. Now you want to do the following:

`gmx pdb2gmx -f 1fxk_clean.pdb -ff gromos54a7 -o first_protein.gro -water spce`

This will produce a topology which GROMACS needs to simulate your protein. If you are using your own protein, simply replace '1fxk.pdb' with the name of your protein file. This command would normally prompt us to select our desired force field for the simulation, however, since we are not using a standard one we have provided it to gromacs with '-ff gromos54a7'. If you want to use a different force field to this you can either provide it directly in a similar way, or exclude this section of the command and choose from the ones GROMACS offers.

Next, it's time to add that water right back (sort of)

### Solvating your protein

To better simulate real life biology, we want to make sure our protein is in solution, rather than in a vacuum. This involves two steps. First, you want to build a box, like so:

`gmx editconf -f first_protein.gro -o boxed_protein.gro -c -d 1.0 -bt dodecahedron`

We'll be using a dodecahedron in this guide to reduce the number of molecules you'll have to simulate, giving quicker results. Now, we want to fill that box with solvent:

`gmx solvate -cp boxed_protein.gro -cs spc216.gro -o solv_prefoldin.gro -p topol.top`

To continue better approximating real life biology, next it's time to get rid of any net charges in our system

### Neutralising net charges

Gromacs frequently uses .mdp files to specify the parameters it cares about for a particular command, such as assembling a binary input file as we are doing here. These parameter files can be very complex and quite intimidating for a newcomer. As such, we are once again going to borrow from Justin's amazing work and use the mdp files he uses in his tutorial. These are designed for a different force field (OPSL) than we are using, but I have verified that they work for our purposes as well. However, as he stresses, the same may not be true if you are using a different force field, such as AMBER. The mdp file we need for this section can be found [here](http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp). Make sure to save it in the same directory as all the files you have generated so far. Now, run these commands.

```
gmx grompp -f ions.mdp -c solv_prefoldin.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o -solv_ions_prefoldin.gro -p topol.top -pname NA -nname CL -neutral
```

This will add NA or CL ions as appropriate until the charge in our system is neutralised. When prompted, enter 13 'SOL' so that we are only replacing the solvate with these ions and not our protein. Now, we can move onto energy minimisation.

### Energy Minimisation

Now that we've fixed the proteins molecular environment, it's time to focus on it's structure. We want to make sure that the protein is in a regular relaxed state before simulating. To do this, we are going to apply energy minimisation to our protein.
We'll need a new mdp file for this step, which can be found [here](http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp). Save it like before, then run this

```
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
```

The -v here is especially important if you are impatient like me, as it will give you an estimate of when this process is finished. This mdun shouldn't take more than a few minutes. 

### Equilibration

We've now stabilised the solvent and the protein, but we still need to make sure that they will play nice with each other. To do that, we will perform two rounds of equilibration, first stabilising the solute with regards to temperature, and then with regards to pressure and density. In practice this looks very similar to what we did in the last step. Starting with temperature, you'll need to download this [new mdp file](http://www.mdtutorials.com/gmx/lysozyme/Files/nvt.mdp) and then run the following commands.

```
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt
```

Then for pressure and density, with [yet another mdp file](http://www.mdtutorials.com/gmx/lysozyme/Files/npt.mdp):

```
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr
gmx mdrun -v -deffnm nvt
```

Both of these commands should take around half an hour to run, but times may vary depending on the power of your machine. It's important to note here that this guide isn't spending too much time talking about each step since it's aimed at getting you simulating as fast as possible, but after each of these stabilisation steps there are a few useful sanity checks that you can do to make sure that your simulation is going smoothly, which Justin details in his more indepth guide linked previously. Note that for these checks you'll want to have installed xmgrace, which is covered in the analysis section of this tutorial.

### Simulation time

Alright, now we are finally ready to produce our simulation. We'll need [one final mdp file](http://www.mdtutorials.com/gmx/lysozyme/Files/md.mdp) for this. Once that's downloaded, we're ready to go. Simply type this:

```
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -v -deffnm md
```

This is the time heavy part of MD, and will take a few hours to run. The mdp file used will create a 1ns simulation, but if you want a longer or shorter ssimulation just change the nsteps parameter as desired. If you've used the -v flag you should see an estimate of when the run will be finished after a few minutes. So, feel free to take a minute to take a break, pat yourself on the back for getting this far, and prepare to visualise your protein.

<a name="vis"></a>
## Visualising your simulation

### Installing a protein viewer

Many tools are freely available for viewing proteins, such as PyMol and Chimera. I prefer Chimera due to its versatility and ease of use. Additionally, it has quite a few analytic methods built in to start with. The linux version can be downloaded [here](https://www.cgl.ucsf.edu/chimera/download.html)(you want the 64 bit version. Once downloaded, run the following:

```
chmod +x chimera-1.14-linux_x86_64.bin
./chimera-1.14-linux_x86_64.bin
```

If the second command doesn't work you, retry with 'sudo' at the front and enter the admin password when prompted.

### Making a movie

Now that you have Chimera installed, it's time for the cool part: viewing your simulation. Simply open Chimera and select Tools > MD/Ensemble Analysis > MD Movie. This will open a new window. For trajectory format specify GROMACS and then select md.tpr and md.xtc for your run input and tranjectory files respectively. These should have been generated by your final simulation. When you click OK Chimera will analyse these files and build every frame of your simulation for you to visualise.

Speaking of visualisation, let's make your protein a little prettier. Chimera provides a ton of options for colouring your protein. Probably the two most useful ones can both be found under Tools > Depiction. Color Secondary Structure is very useful for understanding the broader structure of your protein, differentiating helixs from strands and strands from coils. However, if you'd rather your protein more closely resemble the images you see in the Protein Data Bank, you can use the Rainbow option, which can colour per residue or per chain.

There are a ton more nifty things you can do in Chimera that might be of interest but, in an effort to keep this guide somewhat brief, I'll simply link the following guide as a [handy reference for some of the basic techniques](https://www.cgl.ucsf.edu/Outreach/Tutorials/GettingStarted.html)

To save your simulation as a movie is trivial in Chimera. Simply go to the box with the playback controls for the simulation and select File > Record Movie. You'll then see a window where you can name your movie and choose your preferred format. Once you've entered your preferences, Chimera will generate your movie for you! 

Congratulations on successfully running and visualising your first Molecular Dynamics Simulation! Now, let's see what that simulation can tell you about your protein.

<a name="dssp"></a>
## Secondary structure analysis with DSSP

Now that you have successfully simulated your protein, there's all sorts of analysis that you can perform on it. In this tutorial we will focus on just one example, concerning the secondary structure of our protein. We want to see how the secondary structure of our phosphorylated protein compares to the wild-type version, especially in the region around the residues we have phosphorylated. To do this, we'll first need to download three more tools, two of which you'll find useful for most analysis, not just this one. If you've just completed your first simulation on your provided protein, this will require you to repeat those steps with the modified protein produced by Vienna-PTM.

### Installing Python 3

Python is an easy to learn programming language which plays a huge role in data science. This tutorial only needs it to run a script I wrote for comparing secondary structure differences between simulations, and so you don't need to know the language at all. However, as with the command line, learning the basics is highly recommended if you plan to spend any significant amount of time fiddling around with bioinformatics/data science. Python should be installed by default on most linux machines, but in case it isn't run:

`sudo apt install python3.7`

### Installing xmgrace

Xmgrace is a simple handy plotting software that you will find useful not only for a variety of different analysis, but also in the sanity checks mentioned during the Running your MD section. Installing it is trivial, simply run:

`sudo apt-get install grace`

Xmgrace can be used to view .xvg files. To use it, just type:

`xmgrace plot.xvg`

### Installing DSSP

Finally, DSSP is the tool we will use to assign secondary structures to our protein in order to compare how they change between simulations. To set it up, do the following:

```
sudo apt-get install dssp
sudo ln -s /usr/bin/mkdssp /usr/bin/dssp
```

We are now ready to analyse our proteins.

### DSSP pipeline

First, we'll need to assign secondary structures to our proteins with DSSP. Let's start with our wild-type protein. Run the following command with the md.xtc and md.tpr from the wild-type simulation.

`gmx do_dssp -f md.xtc -s md.tpr -o ss.xpm`

This command by itself can give you a plot of how the overall secondary structure of the simulation changes over time (just add -sc scount.xvg to the above command). However, if you want to look at specific residues, rather than the structure as a whole, you'll need to do another step. For this next step we will again be thanking Justin, the GROMACS God. He wrote a perl script which can be found [here](https://osf.io/2f3bg/?pid=bafn4) which when run will produce a plot labelled 'summary.xvg', contaning the average percentage of each secondary structure element per residue. To run it, download the script and then enter:

`perl plot_ss_xpm.pl ss.xpm`

Great! We're almost there. Let's quickly rename that plot to something else so that we can repeat these first two steps on the modified protein. How about something descriptive like 'summary_WT.xvg':

`mv summary.xvg summary_WT.xvg`

That's the same command used for moving files. It just happens to also be the best way to rename things. Alright, now repeat the first two steps on your modified protein and call it something else, perhaps 'summary_mod.xvg'. Now comes the final part. Now we are going to use a simple python script I wrote which can be downloaded at the top of this page. Once done, run:

`python3 comp_ss.py`

This script is designed to analyse the average secondary structure of a modifed residue in a polymer, and compare that to the wild-type form of the protein. It includes the 5 closest residues to either side to give a better sense of the trends in the chosen protein region. It will prompt you for five things - which residue you want to look at, how many chains are in your protein, which file you want to test, which one you want to compare it to, and what you want to call your output file. For this tutorial, enter the number of the residue you modifed, 1, summary_mod.xvg, summary_WT.xvg and comp.xvg respectively. The plot produced when viewed in xmgrace should look something like this:

![Temporary Image](https://user-images.githubusercontent.com/47711697/80465810-02db4d80-897f-11ea-9365-c0c88de0a8be.png)
