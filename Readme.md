[DRAFT]
=======

- [Install](#installation)
- [Configure](#configuring-the-pipeline)
- [Run](#running-the-pipeline)
- [Project Status](#project-status)
- [Bugs](#reporting-bugs)

Installation
============

quick note
----------
For the sake of brevity, all instructions assume that you are using a BASH shell. If you have made an informed decision to do otherwise, I assume you know enough to translate everything into your own environment. If you have been forced to do otherwise by some powers that be, poor soul, ask your systems administrator (or local friendly nerd) for help, or take this as a learning opportunity.

the hcp pipeline
----------------
This code has been tested against HCP Pipeline v3.0RC3 (commit 058c132fc, Tue, Jan 14 2014). You'll have to make sure that the this or a later compatible version of the HCP code is installed on all of the machines on which you want to run the workflow. Installation is fairly easy - as long as you already have all of the HCP Pipeline's dependencies installed :). Check the HCP Pipeline readme.txt for more information on how to get that done, paying special attention to FSL and FreeSurfer versions, and installing all of the dependencies of gradient_unwarp.py. If you have multiple versions of gradient_unwarp.py on your machine - be careful! Make sure that the version you call from the command line is the first version found on your python path, otherwise you might see some crashes.

One final note: because the HCP Pipelines include some pretty large files, your systems admin would probably appriciate it if there weren't a new installation for every user. Check around with anyone else who uses the systems you're planning to use who might also use the HCP Pipelines. If they're already installed, it'll save you some pain.

hcp_preproc (this project)
--------------------------
This project was developed for python versions > 2.7 and < 3.0, Nipype > 0.9.1.

To install, clone this repository to your machine, and update your PATH and PYTHONPATH variables:

```bash
export PATH=$PATH:/path/to/hcp_preproc/hcp_prep
export PYTHONPATH=$PYTHONPATH:/path/to/hcp_preproc
```

To install the python requirements, use the requirements.txt file and pip:

```bash
pip install -r requirements.txt
```

This may take more or less time depending on how many of the dependencies are already installed, the type of environment you're running, and your own level of expertise. If you're working on a community machine, talk to your systems administrator about the contents of this file, whether or not it is already installed, and any modifications that you may need to make to your environment to make sure that it runs correctly. 

environment variables
---------------------
The HCP Pipelines make heavy use of environment variables - most of that is taken care of by the nipype workflow. But there are still a couple of variables that it's important to set correctly: $FREESURFER_HOME and $FSLDIR. What's more, it's important that you call the FSL and FreeSurfer setup scripts from your .bashrc or .bash_profile file. Check the HCP Pipeline readme for information regarding which particular versions of FreeSurfer and FSL their code currently targets. 

running the pipeline
====================
The file that contains the nipype workflow, pipe.py, can be called as a command line script. See the section on configuration, then when you're ready to run pass the -r argument, along with any others you choose to use (see below).

For help, call:

```bash
pipe.py --help
```

configuring the pipeline
========================

We currently use [configobj](https://pypi.python.org/pypi/configobj/) to write and read files. pipe.py includes some tools to help you build and update config files pretty quickly, but since they're plain text you can always open them up with a text editor and change them by hand (more on this below).

To build a new config file, call pipe.py with the -i or --init argument. You'll be walked throught the creation of a new config file. You'll want to have already downloaded your data, and you should be sure that you have run the freesurfer and fsl setup scripts. You'll also need to know the path to the HCP Pipeline code. Let's quickly walk through the config steps as they are now. I'll take some time to discuss each of the questions, and how to figure out the appropriate answers.

step-by-step
------------

We start by initializing a new config file.

```
pipe.py --init
New config file name [hcp.conf]: 
```

Decisions decisions - what shall we call the new config file? The square brackets mean that whatever they contain is the default value. So if we just press return, we'll choose the name "hcp.conf". Sounds good to me, so let's just press return. If the file already exists, the script will exit.

```
The subjects directory should contain all raw data for all subjects.
Subjects Directory [./]: /data/nil-external/ccp/MOOD_RISK/DICOMs
```

The subjects directory is the lowest directory below which we can find all of your experiment's dicoms. So, if you had dicoms sorted into /some/dir/sub_a, /some/dir/sub_b, etc., then the subjects directory would be /some/dir. Here, I've chosen something appropriate for my current experiment.

```
The DICOM template should be a format string for a glob which, when combined
with an individual subject ID, will get us all of the subject's DICOM files.
DICOM template [data/raw_dicom/%s/*.dcm]: DR%s/SCANS/*/DICOM/*.dcm
```

This is perhaps the trickiest one to answer, so I'm going to walk you through it step by step. The goal here is to tell the workflow how to find all of the dicoms for a given subject. We do that by giving it what is called a format string (google FMI), which allows us to substitute in each subject number in place of what is called a string format specifier, in this case '%s'. If that all seems like jargon, just follow along and hopefully things will start to make sense.

The first step towards finding my format string is understanding how my data is organized. Let's take a look. We know my subject folder is /data/nil-external/ccp/MOOD_RISK/DICOMs, so let's make a quick exploration of that directory's organization:

```bash
$> ls /data/nil-external/ccp/MOOD_RISK/DICOMs
DR060  DR061  DR064
$> ls /data/nil-external/ccp/MOOD_RISK/DICOMs/DR060/
SCANS
$> ls /data/nil-external/ccp/MOOD_RISK/DICOMs/DR060/SCANS/
1  10  11  12  13  14  15  16  17  18  19  2  20  21  22  23  24  25  26  27
28  29  3  30  31  32  33  34  35  36  37  38  39  4  40  5  6  7  8  9
$> ls /data/nil-external/ccp/MOOD_RISK/DICOMs/DR060/SCANS/1/
DICOM/     SNAPSHOTS/ 
$> ls /data/nil-external/ccp/MOOD_RISK/DICOMs/DR060/SCANS/1/DICOM/
DR060.MR.Barch_MoodRisk.1.1.20131205.173445.s02tx.dcm    DR060.MR.Barch_MoodRisk.1.3.20131205.173445.19qm3q2.dcm
DR060.MR.Barch_MoodRisk.1.2.20131205.173445.10iky0u.dcm  scan_1_catalog.xml
```

I might do the same looking into other subjects to verify that the organization is consistent. If that's not the case, you'll need to do some cleanup, then come back to this point. Assuming that we're satisfied on this point for now, we can see that a valid path to a specific dicom might be:

```bash
/data/nil-external/ccp/MOOD_RISK/DICOMs/DR060/SCANS/1/DICOM/DR060.MR.Barch_MoodRisk.1.1.20131205.173445.s02tx.dcm
```

So how to we get from here to a list of *all* of the experiments dicoms? First, we make use of the wildcard, *. Because we use this string as what is known as a 'glob,' the character * will be expanded to match any number of characters, with a few exceptions (like '/'). So using this, we can begin to shrink our string:

```bash
/data/nil-external/ccp/MOOD_RISK/DICOMs/DR060/SCANS/*/DICOM/*.dcm
```

So here we've used two globs, one to replace the specific file name (we want everything that ends in .dcm), and another to replace the scan number. So now what about the subject number? This case is a little special. Since we want the script to be able to substitute in specific subject numbers, we use a string format specifier here instead of a wildcard. In the case of this data, the pattern seems to be .../DICOMs/DR<SUBJECT_NUMBER>/SCANS..., so let's put a string format specifier in place of the subject number:

```bash
/data/nil-external/ccp/MOOD_RISK/DICOMs/DR%s/SCANS/*/DICOM/*.dcm
```

One last change. We don't need to include the subject directory prepended to the DICOM template (in fact, it is important that we do not). So let's get that out of there, leaving us with:

```
DR%s/SCANS/*/DICOM/*.dcm
```

Which is what we hand to the script! Moving right along:

```
Subjects should be a comma separated list of subject ids.
Subject list ['']: 060, 061, 064
```

Feel free to provide nothing here. If you want to store a particular list of users to whom this config script should be applied, you can supply them here or later by hand. You can also specify them on the command line when you call pipe.py -r.

After pressing enter, the script will look through all of the DICOM files that it can find. If you need to speed this up, I'll leave it as an exercise for the reader to figure out how to make it so the script only finds one or two subjects worth of data.

```
Checking series names (this may take some time) (41 chunks remaining)...
```

After it gets to 0, we can start in on the fun stuff! The script will print a numbered list of all the Series Descriptions it was able to find. Our job now is to tell it how to use that information to feed our data through the HCP Pipelines. Let's take a look at what it found in my case:

```
Found 25 unique series descriptions.
-------
Series:
-------
0:  AAHScout
1:  AAHScout_MPR_cor
2:  AAHScout_MPR_sag
3:  AAHScout_MPR_tra
4:  BIAS_32CH
5:  BIAS_BC
6:  BOLD_FACE1
7:  BOLD_FACE1_SBRef
8:  BOLD_FACE2
9:  BOLD_FACE2_SBRef
10: BOLD_REWARD1
11: BOLD_REWARD1_SBRef
12: BOLD_REWARD2
13: BOLD_REWARD2_SBRef
14: BOLD_REWARD3
15: BOLD_REWARD3_SBRef
16: BOLD_TEST
17: BOLD_TEST_SBRef
18: FieldMap
19: Localizer
20: Localizer_aligned
21: SpinEchoFieldMap_AP
22: SpinEchoFieldMap_PA
23: T1w_MPR_08mm
24: T2w_SPC_08mm
```

Alright - the first couple of questions are pretty easy. Your SBRef images might be called Scout or something else instead:

```
Which series do you use for 'bold'?
[None] or comma separated values 0-24: 6,8,10,12,14

Which series do you use for 'bold_sbref'?
[None] or comma separated values 0-24: 7,9,11,13,15
```

For the next two, you'll often provide the same answer twice:

```
Which series do you use for 'fieldmap_magnitude'?
[None] or comma separated values 0-24: 18

Which series do you use for 'fieldmap_phase'?
[None] or comma separated values 0-24: 18
```

If you collected Spin Echo Fieldmaps, they're probably either PA/AP, or RL/LR. Just leave blank responses for those you didn't collect:

```
Which series do you use for 'fieldmap_ap'?
[None] or comma separated values 0-24: 21

Which series do you use for 'fieldmap_lr'?
[None] or comma separated values 0-24: 

Which series do you use for 'fieldmap_pa'?
[None] or comma separated values 0-24: 22

Which series do you use for 'fieldmap_rl'?
[None] or comma separated values 0-24:
```

T1/T2 are pretty easy:

```
Which series do you use for 't1'?
[None] or comma separated values 0-24: 23

Which series do you use for 't2'?
[None] or comma separated values 0-24: 24
```

The next one, 'polarity_swapped,' requires some explanation. In some experiments, you might acquire both RL and LR (or AP and PA) BOLD images. We're always trying to improve the list of values that the workflow derives at runtime, but for various reasons detecting this switch is difficult to do reliably. So we need your help, dear user. If you acquire images with opposing polarities, choose one of them, say "LR" if you're RL/LR, or "PA" if you're AP/PA, to call "swapped." We don't see that in this experiment, so we'll leave this blank. But if I did have two of each bold image, one with suffix _AP and one with suffix _PA, I would list here the numbers of all of the "swapped" series, ie those that ended with _PA.

```
Which series do you use for 'polarity_swapped'?
[None] or comma separated values 0-24:
```

If the FreeSurfer and FSL environment variables are set correctly, the next two questions should provide correct default options. This works for me, so I don't supply an alternate value:

```
Path for FREESURFER_HOME [/usr/local/freesurfer]: 

Path for FSLDIR [/usr/share/fsl/5.0]:
```

Cool beans. For the next one, I need to know the location of the HCP Pipeline code:

```
Path to HCP Pipelines dir [ ]: /scratch_cl1/hcp_pipe/Pipelines
```

Pop quiz: do you know the resolution of your structural data?

```
What is your structural image resolution (mm)?
[Skip] or one of (.7, .8, 1): .8
```

Now it asks if I want to use the default template files for my t1 resolution, and the default config files. Yes and yes.

```
Use default template files for resolution 0.8 [y]/n? 

Use default config files [y]/n?
```

Getting close! The next step was added to handle cases in which you may acquire more than one Spin Echo Fieldmap. In these cases you have a choice between two policies. Either we'll just choose the first set we find in each session (AP/PA or RL/LR pair), or we'll do something a little subtler. If you want the simple option, just choose 'first.' If you choose 'most_recent,' then for each bold we'll either use the SE Fieldmap pair most recently acquired prior to that particular BOLD, or the first one acquired thereafter. In our case, we want the more complex option:

```
If you have more than one ep fieldmap set, you may either
want to use the first, or always use the most recent.
Which policy would you like to use:

0: first
1: most_recent

Select 0-1: 1
```

This next one is pretty self explanatory:

```
If you collect multiple t1 or t2 images, and averaging them yields warped
results, try blocking structural image averaging.

Block averaging of structural images [y]/n? y
```

The last thing that the config setup script does is to make a guess at the unwarp direction.

```
Very weak guess that your primary unwarp direction is y.
Did I mention this is a GUESS?

When finished, please open your config file check the value for ep_unwarp_dir.
```

At this point, unfortunately, we do not have fully automated derivation of unwarpdir in place. So this really is just a guess. You should look over the final results carefully, taking care to check for any untoward distortions (like swirls, or unlikely overall shapes). If you see any of these things, try opening the config file and changing the unwarp direcion from x to y, from -x to x, or from y to -y, or any other combination. You might even try z in a pinch - but I wouldn't try it first!

That's it for the config script at this point. To re-run the later part (optionally skipping the series mapping), call pipe.py -u or pipe.py --update.

by hand
-------

If you choose to edit the config file by hand, please note that we use configobj in unrepr mode, which means that it expects the values to be in python format. The means that strings "need to be quoted," lists ["must","be","in","brackets"], and other primitives like 2.3 (floats), 2 (integets), True and False (boolean values) can be used just as you would use them in a python script. Because you know and love to write python, riiiight? Right!

customizing node configuration
------------------------------

If you choose to get fancy, it might help to know a little more about how configuration works. In particular, the nifti_wrangler and each of the hcp nodes have a special trick: just before running, they check to see if any of their values have been set in the config file. To do this, I first check to see if the config file contains sections with names matching any of these nodes. This includes: nifti_wrangler, pre_freesurfer, freesurfer, post_freesurfer, volume_processing, and surface_processing. Then I check to see if any of the settings in that section have names that match any of the respective node's input values. If so, I try to set the value. We can see this at work in the default config script, where we supply several arguments to nifti_wrangler:

```
[nifti_wrangler]
ep_fieldmap_selection = 'most_recent'
block_struct_averaging = True
ep_unwarp_dir = 'y'
```

For a full list of the inputs of each node, check out interface_docs.txt in the docs directory. This file is generated by looking at the current interface specifications, so it should be accurate.

Take note that not all attributes are config-file-settable. Those that are not include those that are derived at runtime. This list is undocumented at the moment, so you can either open the code or bug me about it.

validation
----------

Some time in the future, we'll put in some nice config file validation stuff. Maybe :).

Project Status
==============

Alpha release.

This software is not a replacement for knowing what you're doing, and you should inspect all of the pipeline's output logs carefully. We make no promises whatsoever at this point that the output of this pipeline is worth using.

Reporting Bugs
==============

Please use the githup isses tab for bug reports and feature requests.
