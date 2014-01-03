[DRAFT]
=======

- [Install](#installation)
- [Configure](#configuring-the-pipeline)
- [Run](#running-the-pipeline)
- [Project Status](#project-status)
- [Bugs](#reporting-bugs)

Installation
============

the hcp pipeline
----------------
At this point we're targeting the development vesion of the HCP Pipeline. You'll have to make sure that the latest version is installed somewhere on the machines you'll be using to run the pipeline. If we ever lock onto one version, we'll make a note of it here.

I currently use a fork of the main HCP repository, with a couple of very small fixes for a couple of permissions issues. You're free to use that, but I should really put in a pull request soon.

hcp_preproc (this project)
--------------------------
This project was developed for python 2.7+. Has not been tested in python 3+. To install, clone this repository to your machine, and update your PATH and PYTHONPATH variables:

```bash
export PATH=$PATH:/path/to/hcp_preproc/hcp_prep
export PYTHONPATH=$PYTHONPATH:/path/to/hcp_preproc
```

To install the python requirements, use the requirements.txt file and pip:

```bash
pip install -r requirements.txt
```

Please ensure that you install Nipype > 0.9.1.

Running the Pipeline
====================
The file that contains the nipype workflow, pipe.py, can be called as a command line script. See the section on configuration, then when you're ready to run pass the -r argument, along with any others you choose to use (see below).

Configuring the Pipeline
========================

We currently use Python's ConfigParser to write and read plain text configuration files. pipe.py includes some tools to help you build and update config files pretty quickly, but since they're plain text you can always open them up with a text editor and change them by hand.

To build a new config file, call pipe.py with the -i argument. You'll be walked throught the creation of a new config file. You'll want to have already downloaded your data, and you should be sure that you have run the freesurfer and fsl setup scripts. You'll also need to know the path to the HCP Pipeline code.

TODO: a step by step of the config process. hopefully not so urgent.

While this gets you a fairly complete config file, there is also some limited capacity for passing arguements from the config file directly to the nipype workflow's nodes. To see the arguemnts that these nodes can accept, check out the /docs/interface_docs.txt file. And of the input arguments that are strings can currently be specified in the config file. Yes, I know this is kinda lame right now. I plan to move away from ConfigParser, and into something that uses JSON, so we can handle more than strings.

If you have something about your project change, and you want to use the same tool to update your config file, just call pipe.py -u.

We try to derive as much information as we can about your scans at runtime. Decisions that are explicitly surfaced in the config file are there because you should pay attention to them. The unwarpdir, for example, might need to be reversed if you find that the distortion correction leaves your images twice as distorted.

Project Status
==============

Alpha release.

This software is not a replacement for knowing what you're doing, and you should inspect all of the pipeline's output logs carefully. We make no promises whatsoever at this point that the output of this pipeline is worth using.

Reporting Bugs
==============

Please use the githup isses tab for bug reports and feature requests.
