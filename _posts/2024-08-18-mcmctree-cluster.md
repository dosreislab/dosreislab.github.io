---
layout: post
title:   "Using MCMCtree in QMUL's Apocrita"
author: "Mario dos Reis"
---

In this tutorial we will run a simple MCMCtree analysis at [QMUL's Apocrita high-performance computer cluster](https://docs.hpc.qmul.ac.uk/). QMUL students and staff can request an Apocrita account [here](https://docs.hpc.qmul.ac.uk/intro/hpc-account/). Apocrita uses job scheduling software for fair allocation of computational resources to users of the cluster. To submit a computational analysis to Apocrita, you must first write a script with instructions for you analyses. This script (known as a job) is submitted to the scheduler's job queue. When a computing node with enough resources becomes available, your job is released from the queue and allocated to the node, where it is allowed to run. There is extensive documentation at Apocrita's home site on how to use the cluster and the scheduler.

This tutorial assumes you are familiar with the command line, SSH, MCMCtree (i.e, you have already done the MCMCtree tutorials), and basic shell scripting with BASH.

We will follow these steps:

1. Download and compile MCMCtree in Apocrita.
2. Prepare data and control files for our MCMCtree analysis.
3. Submit the MCMCtree job to Apocrita's job queueing system.
4. Check the MCMC results once the job has finished.

# 1. Download and compile MCMCtree in the cluster

First, let's login onto the cluster. Please replace `xyz123` below with your QMUL user name.

```
ssh xyz123@login.hpc.qmul.ac.uk
```

Downloading and compiling software are resource intensive tasks. Thus, we will first log onto the queue, so that a suitable node is allocated to us. In that way, we will not interfere with the work of other users of Apocrita. In the command line type:

```
qlogin
```

The command will wait until a node becomes available, and then it will log you onto an interactive job.

When working in the cluster it is very important to keep our files well organised. We will create three directories, `bin`, `software` and `projects`. The `bin` directory will contain scripts and executable files we need to use often; `software` will contain programs for major analyses (such as MCMCtree, RAxML, or MrBayes), and `projects` will contain our analyses, i.e., it is the place where we will run our _HPC jobs_.

Make the three directories with the `mkdir` command:

```
mkdir bin software projects
```
You can verify the directories have been created using `ls`:

```
ls
```

We will now download and install MCMCtree. First we go into the `software` directory and then we use `git` to clone the MCMCtree software:

```
cd software
git clone https://github.com/dosreislab/mcmctree.git
```

You now have downloaded dos Reis' lab modified version of MCMCtree. The modifications have been made to simplify certain HPC tasks (such as running unlimited checkpoint cycles, for computational heavy analyses, such as trees with too many species). Let's go into the `src` directory of mcmctree and compile the executable.

```
cd mcmctree/src
```

The `src` directory contains the _source_ files with the code that makes up MCMCtree. MCMCtree is written in C, an efficient computer language that needs to be compiled into machine code. We will use the Intel compiler to generate fast, effcient executables that are fine tuned to QMUL's cluster architecture. The instructions for the compiler are inside the `Makefile` file. Using your favorite text editor, open `Makefile` (for example, you can try `module load nano; nano Makefile`). Edit `Makefile` so that it looks like this:

```
PRGS =  mcmctree
CC = icx # cc, gcc, cl # if running on an intel cluster, check if you can compile with the intel compiler, icc

#CFLAGS = -O3
CFLAGS = -fast # this flag is for the intel compiler

LIBS = -lm # -lM

all : $(PRGS)

mcmctree : mcmctree.c  tools.c treesub.c treespace.c paml.h
        $(CC) $(CFLAGS) -o $@ mcmctree.c tools.c $(LIBS)
        $(CC) $(CFLAGS) -o infinitesites -D INFINITESITES mcmctree.c tools.c $(LIBS)
```

Three lines have changed. The second line is changed to `CC = icx`, which tells the system to use the Intel compiler. Line five has a comment character `#` added at the beginning. Line six is changed to `CFLAGS = -fast`, which tells the Intel compiler to generate fast code. 

[**Note:** `icx` is the new name for the intel compiler executable, with the old `icc` now deprecated.]

Go back to the command line and ask the system to list the Intel compilers available:

```
module avail intel
```

As of April 2025, the command above lists the following compilers:

```
intel-classic/2021.10.0  intel-mkl/2024.1.0  intel-mpi/2021.12.1  intel-tbb/2021.9.0-gcc-12.2.0  intel/2023.2.4  intel/2024.1.0  
```

The mpi compilers are used to generate MPI parallelisable code. MCMCtree does not use MPI so we don't need to use these. The latest compiler available in my list is `intel/2024.1.0`. Load the compiler with:

```
module load intel/2024.1.0
```

You can use a later version of the compiler if available. Now we can compile MCMCtree:

```
make
```

The `make` program reads the instructions in `Makefile`, and calls the `icx` compiler as appropriate. Let's print the contents of the `src` directory to our screen, which now contains the `mcmctree` executable:

```
ls
```

This should produce the following output:

```
Makefile  infinitesites  mcmctree  mcmctree.c  paml.h  tools.c  treespace.c  treesub.c
```

That is, two new executable files have been created, `mcmctree` and `infinitesites`. The latter is a special version of MCMCtree used for asymptotic analysis assuming an infinite amount of data.

**Tip:** To get out of the interactive queue started with the `qlogin` command, simply type `exit` in the command line.

**Tip:** You can find the manual for a command using `man`. For example, try `man qlogin` or `man ls` to learn more about the `qlogin` and `ls` commands.

**Tip:** Sometimes an ssh session to the HPC cluster becomes unresponsive. For example, you went away from your computer and when you came back you found you could not type anything into the terminal. If that happens, simply press the follwing keys in sequence: `Enter` `~` `.`. That will terminate the shell process. You will then need to login again. You can try it now.

# 2. Prepare data and control files for our MCMCtree analysis

We will estimate divergence times on a four-species hominid dataset including human, Neanderthal, chimp and gorilla. First, we create the appropriate project directory.

```
cd ~/projects/
mkdir hominid
```

The alignment, tree file, and control file we need for our analyses are in the mcmctree git directory, so we copy them (it is not good pratice to carry out analyses within git directories, as these should be used for software development and testing only).

```
cd hominid
cp ~/software/mcmctree/data/ape4s/* .
ls
```
You should see the following files:

```
ape4s.phy ape4s.tree mcmctree.ctl
```

We will run four MCMC chains in parallel. So we need to create four directories:

```
mkdir `seq 1 4`
```

The `seq` command simply generates sequences of numbers. In this case, 1, 2, 3, and 4, which are used by the `mkdir` command as the names of the directories to be generated.

We now need to link the alignment, tree and control files inside the directories that will run each independent MCMC chain:

```
for d in `seq 1 4`
  do 
  cd $d
  ln -s ../mcmctree.ctl .
  ln -s ../ape4s.phy .
  ln -s ../ape4s.tree .
  cd ..
done
```

The sequence of commands above loops among the four directories, goes inside each one, creates soft links to the control, alignment and tree files, and then leaves each directory before moving to the next one. Soft links are useful: instead of copying each file inside each directory, we simply put a link that tells the system where the true file is located. Soft links are efficient because they use very little disk space, and they are synchronised to the main file. If you edit the tree file (for example, you change a fossil calibration), all the soft links are updated automatically.

You can check the commands above worked correctly by listing the contents of one the directories. For example:

```
ls -l 1
```

Which should generate an output similar to this one:

```
lrwxrwxrwx 1 xyz123 qmul 12 Aug 18 19:31 ape4s.phy -> ../ape4s.phy
lrwxrwxrwx 1 xyz123 qmul 13 Aug 18 19:31 ape4s.tree -> ../ape4s.tree
lrwxrwxrwx 1 xyz123 qmul 15 Aug 18 19:31 mcmctree.ctl -> ../mcmctree.ctl
```

Note the arrows indicating the files are soft links.

# 3. Submit the MCMCtree job to Apocrita's job queueing system

Now we prepare our job submission script. First we create the job file inside the `bin` directory:

```
cd ~/bin
touch mcmctree-job.sh
```

Now using your favorite text editor, open the file, which will be empty, and add the following lines:

```
#! /bin/bash          # use the BASH shell
#$ -cwd               # change to current working directory
#$ -V                 # export environment to job
#$ -j y               # merge standard out and standard error
#$ -l h_rt=01:00:00   # set time limit
#$ -l h_vmem=100M     # amount of memory needed 

FMCMC=mcmc.txt

cd $SGE_TASK_ID

# Check if mcmc.txt file exists, if it does, abort
# mission, otherwise run mcmctree
if [ -f "$FMCMC" ]; then
    echo "ERROR: $FMCMC exists, nothing done."
else
    # This is MCMCtree 4.9j compiled with icc -fast and with NS = 1500
    date >START 
    ~/software/mcmctree/src/mcmctree mcmctree.ctl   
    date >DONE
fi
```

The first six lines have general instructions for the scheduler, such as the shell type (BASH), and the amount of memory and time needed to run the jobs. Here we request 1 hour of time and 100 Mb of memory. It is important to get these numbers right, as the scheduler uses them to decide on the priority of your job in the queue. A job that needs 10 days to run may spend longer in the queue than a job that requires just one hour.

Line `FMCMC=mcmc.txt` creates a variable with the name of the output of our MCMC chain. Next line `cd $SGE_TASK_ID` tells the script to change to the directory defined in `SGE_TASK_ID`. This variable will contain the directories 1 to 4 (see below). 

The next block, starting with `if [ -f "$FMCMC" ]; then` checks whether the `mcmc.txt` file exists. If it does, the job is aborted. This is done for safety. A common mistake is to submit a job in the wrong directory. If that happened, then you would have two MCMCtree runs interfering with each other. When designing submission scripts, it is important to think about safety features to make your analyses robust to errors. 

The next block, starting with `else` calls the `mcmctree` executable that we just compiled. The `date` command prints the current date and time to the screen. Here we catch the output of `date` and we redirect it to files called `START` and `DONE`. In this way, we record when our `mcmctree` analysis started, and when it finished. 

Now we submit our script. Because we need to carry out four analyses with identical setting of parameters, we will submit an array job. In the command line type:

```
cd ~/projects/hominid
qsub -t 1-4 ~/bin/mcmctree-job.sh
```

The `qsub` command submits our script to the queue. The `-t` flag activates the array job option, and tells the scheduler we are running four jobs with numbers 1 to 4. Internally, the scheduler will set the `SGE_TASK_ID` variable sequentially to 1, 2, 3, and 4. Thus, line `cd $SGE_TASK_ID` in our submission script will have the effect of changing to the appropriate directory as set by the scheduler.

Our array job is now in the queue. To check our job status type

`qstat`

This will list the array job. Initally, the job is queued, indicated with `qw`, in the output of `qstat`:

```
3803592 0.00000 mcmctree-j xyz123       qw    08/18/2024 19:39:33                                                                   1 1-4:1
```

The job will eventually be split into the array, and you will see four jobs with running status `r`. These jobs will run very quickly, so you may not have chance to see them split in the `qstat` output.

# 4. Check the MCMC results once the job has finished

Congratulations, you have submitted your firt job to the cluster! Once the jobs have finished, you can check the results. You will see four files have been created:

```
mcmctree-job.sh.o3803592.1 mcmctree-job.sh.o3803592.2 mcmctree-job.sh.o3803592.3 mcmctree-job.sh.o3803592.4
```

The long number, `3803592`, is the ID number the scheduler assigned to the job, which will look different in your case. The last number is the array number. These files contain the output that MCMCtree would usually print to the screen. If anything went wrong with your analysis, examining these files is a good place to start, as MCMCtree may print errors and issues there.

The next step is to check your MCMC chains converged. If you want, you can transfer the result files to your local computer to check them with Tracer, R, or any other useful software. The command `rsync` (for remote synchronisation) is useful for this. It is available in Linux and MacOS. In Windows, you can use `rsync` within Windows' Linux Subsystem. Alternatively, you can use many of the GUI SSH options available.

Within the cluster, you can do a quick check by inspecting the output files (`mcmc.txt`, `out` and `FigTree.tre`) generated by MCMCtree. For example,

```
cd 1/
less out
```

will let you examine the `out` file. The `less` command is very useful to examine text files as it allows you to search for specific strings of text and move quickly to specific parts of the document.

Now that you have finished your HPC tutorial, take the time to go through [Apocrita's documentation on the web](https://docs.hpc.qmul.ac.uk/), where many job example scripts are available. Also, learning BASH (shell scripting) as well as the commands `awk`, `sed` and `grep` will improve your skills in preparing scripts and processing text output. The `vim` text editor is also extremely powerful. There are many tutorials on the web available for all these commands. It is worth spending a day learning each one of them. Other simple, yet useful, commands include `rev`, `head`, `tail`, `cat` and `paste`. You can learn more about them with their manual pages (e.g., `man cat`). Enjoy!

_Last Updated 19th Aug 2024._

