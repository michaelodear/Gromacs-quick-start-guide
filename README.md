# Gromacs quick start guide

1. [Setting up Gromacs](#setup) <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.1 [Install with GPU (recommended)](#gpu) <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.2 [Install with CPU](#cpu)
2. [Running your first MD](#md1)
3. [Secondary structure analysis with DSSP](#dssp)

When I was first getting started with GROMACS I was frustrated by how much of the documentation and tutorials were outdated. 
While Justin Lemkul's suite of Gromacs tutorials are a brilliant resource, with the running your MD section being 

Note: This guide is written assuming you are working on a linux machine running Ubuntu

<a name="setup"></a>
## Setting up Gromacs

<a name="gpu"></a>
### Install with GPU (recommended)

<a name="cpu"></a>
### Install with CPU

If you don't have a GPU available, the quickest way to get GROMACS running on Ubuntu is to use the native package installer, sudo apt-get:

`sudo apt-get install gromacs`

Note that you will need admin privileges for this and so will need to enter a password if not running as admin

<a name="md1"></a>
## Running your first MD

<a name="dssp"></a>
## Secondary structure analysis with DSSP
