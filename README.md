# Gromacs quick start guide

1. [Setting up Gromacs](#setup) <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.1 [Install with GPU (recommended)](#gpu) <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.2 [Install with CPU](#cpu)
2. [Preprocessing (Optional)](#pre)
3. [Running your first MD](#md1)
4. [Visualising your simulation](#vis)
5. [Secondary structure analysis with DSSP](#dssp)

When I was first getting started with GROMACS I was frustrated by how much of the documentation and tutorials were outdated. One notable exception to this are [Justin Lemkul's tutorials](http://www.mdtutorials.com/gmx/index.html), which are an excellent resource for running MD simulations that I highly recommend everyone check out, and which the 'Running your first MD section' borrows heavily from. 

This guide's aim is to collate a number of resources from around the web to walk a new user through the process of running their first molecular dynamics simulation, from installation to analysis. It is aimed at people without any programming background and written assuming that you are working on a Linux machine running Ubuntu 18.04. If you don't have access to a linux machine, this [very short guide from Microsoft](https://docs.microsoft.com/en-us/windows/wsl/install-win10) explains how to run Linux on a Windows machine.

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

<a name="md1"></a>
## Running your first MD

### Converting your protein to a format Gromacs understands

### Neutralising net charges

### Energy Minimisation

### Equilibration

### Running your first MD

<a name="vis"></a>
## Visualising your simulation

<a name="dssp"></a>
## Secondary structure analysis with DSSP
