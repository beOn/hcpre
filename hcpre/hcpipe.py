#!/usr/bin/env python
# encoding: utf-8

"""
hcpipe.py

Copyright (c) 2013 Ben Acland
"""

import os
import sys
import getopt

from hcpre.workflows import *

help_message = """
Provides command line tools to build and work on workflow confiuration files,
and for launching the pipeline using a specific config. You can either import
HCPrepWorkflow, then configure and run in from your own code, or use the
built-in command line tools to launch the workflow on your data.

Commands
--------
-h, --help
    Prints out this message.
-i, --init
    Creates a new config file in your current directory.
-u, --update
    Reruns part of the config file setup script on any config file in your current directory.
-g, --graph
    Draws a graph of the nipype workflow to the current directory.
-r, --run
    Runs the workflow using any config file in the current directory.

Parameters
----------
-c, --config (path)
    The config file to use.
-s (comma separated subject numbers)
    The list of subjects who you'd like to put through the workflow. Overrides
    setting in the config file.
-n (integer)
    The number of threads you would like to use. Higher numbers are faster to
    a point, but depending on how large your data is (chances are, it is quite
    large), you may well want to limit yourself to something below 8 for
    starters if you're working on a large server. Default is 1. Ignored if
    you use -p.
-p, --pbs
    Causes nipype to try to use the PBS plugin. For use on the cluster only.
    Experimental.
-o (path)
    The directory to put preprocessed data in. Default is current directory.
-v 
    Verbose. At this point, just makes -g print a graph that expands iterables.
"""

class Usage(Exception):
    def __init__(self, msg=help_message):
        self.msg = msg

def main(argv=None):
    import os
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "pvrhiugs:n:o:c:", ["run", "help", "init", "update", "graph", "config", "pbs"])
        except getopt.error, msg:
            raise Usage(msg="\n"+str(msg))
        # option processing
        update = False
        graph = False
        subs = None
        N_PROCS = None
        run = False
        verbose = False
        c_file = None
        use_pbs = False
        out_dir = os.getcwd()
        for option, value in opts:
            if option in ("-h", "--help"):
                raise Usage()
            if option in ("-i", "--init"):
                setup_conf()
                return
            if option in ("-u", "--update"):
                update = True
            if option in ("-g", "--graph"):
                graph = True
            if option in ("-r", "--run"):
                run = True
            if option in ("-c", "--config"):
                c_file = value
            if option in ("-p", "--pbs"):
                use_pbs = True
            if option in ("-s"):
                subs = [sub.strip() for sub in value.split(",")]
            if option in ("-v"):
                verbose = True
            if option in ('-n'):
                N_PROCS = int_or_none(value)
            if option in ('-o'):
                out_dir = value
        # select config file
        if not c_file:
            try:
                c_file = select_conf()
            except Exception, e:
                raise Usage(msg="Could not find a config file.")
        # update if necessary
        if update:
            update_conf(c_file)
            return
        # make sure we're going to do something
        if not run and not graph:
            raise Usage(msg="Nothing to do...")
        # validate the config for running
        conf = get_config_dict(c_file)
        if not conf:
            raise Usage(msg="Could not parse config file.")
        if not validate_config(conf):
            raise Usage(msg="Invalid config file.")
        # build the workflow, pass it subjects if they were given in the command line
        wk = HCPrepWorkflow(name="hcp_prep_workflow", config=conf)
        if subs:
            wk.subjects_node.iterables = ("subject", subs)
        # set the output dir
        out_dir = os.path.abspath(out_dir)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        wk.data_sink.inputs.base_directory = out_dir
        # graph if you like
        if graph:
            g2u = "exec" if verbose else "orig"
            wk.write_graph(dotfilename="hcp_pipe_graph", graph2use=g2u)
            return
        if not run:
            return
        if use_pbs:
            print "running using PBS"
            # note - the default qsub args are very modest. use this structure
            # to scale up if necessary.
            ques = [
                [[wk.dicom_convert, wk.dicom_info, wk.nii_wrangler],
                    {"qsub_args":"-l nodes=1:ppn=1,mem=1gb,walltime=1:00:00"}],
                [[wk.hc_pre_fs],
                    {"qsub_args":"-l nodes=1:ppn=2,mem=10gb,vmem=10gb,walltime=6:00:00"}],
                [[wk.hc_fs],
                    {"qsub_args":"-l nodes=1:ppn=2,mem=5gb,walltime=24:00:00"}],
                [[wk.hc_post_fs],
                    {"qsub_args":"-l nodes=1:ppn=4,mem=10gb,walltime=4:00:00"}],
                [[wk.hc_volume, wk.hc_surface],
                    {"qsub_args":"-l nodes=1:ppn=4,mem=10gb,walltime=12:00:00"}],
                ]
            for q in ques:
                p_args = dict(q[1],**{"overwrite":True})
                for n in q[0]:
                    n.plugin_args = p_args
            wk.run(plugin="PBS",
                   plugin_args={"qsub_args":"-l nodes=1:ppn=1,mem=1gb,walltime=1:00:00"})
        elif N_PROCS > 0:
            print "running with %d processes" % N_PROCS
            wk.run(plugin="MultiProc", plugin_args={"n_procs" : N_PROCS, "non_daemon" : True})
        else:
            print "running single process"
            wk.run()
    except Usage, err:
        f_str = sys.argv[0].split("/")[-1] + ":"
        lfs = len(f_str)
        f_str = "%s\n%s\n%s\n" % ("-"*lfs, f_str, "-"*lfs)
        print >> sys.stderr, f_str + str(err.msg)
        print >> sys.stderr, "-------------------\nfor help use --help\n-------------------"
        return 2

if __name__ == "__main__":
    sys.exit(main())

