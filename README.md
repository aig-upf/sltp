


# The SLTP Generalized Planning Framework: Sample, Learn, Transform & Plan


## Installation

The whole pipeline runs in *Python3* and relies on some dependencies; most notably:

* The [FS planner](https://github.com/aig-upf/fs-private/) (actually, a simplified version of it).
* Some MaxSat solver such as [OpenWBO](http://sat.inesc-id.pt/open-wbo/)
* The [Tarski](https://github.com/aig-upf/tarski/) planning problem definition module.
* CMake

### Installing required packages

If running on Ubuntu, you will need to install CMake: 
    
    sudo apt-get -y install cmake


### Installing FS

We use the FS planner to fully expand state space samples. It is a bit of overkill, but might prove
useful down the road if we need to use compact FSTRIPS encodings. To install the planner, you need 
the following packages (assuming Ubuntu):


    sudo apt-get install -y --no-install-recommends \
         build-essential g++ python3 git scons libboost-all-dev


Now checkout the `sltp-lite` branch (which contains a simplified version of the planner for our purposes)
and run the build script:

    git clone git@github.com:aig-upf/fs-private.git -b sltp-lite fs-planner
    cd fs-planner
    git submodule update --init
    ./build.py -p

Once the build has finished, you'll need to create an environment variable named `$FS_PATH` to point
to the installation location, e.g. adding this to your `.bashrc` file (change the path to your actual
location):

    export FS_PATH="${HOME}/projects/code/fs"


### Installing openWBO
Here you should simply download and follow 
[OpenWBO](http://sat.inesc-id.pt/open-wbo/)'s installation instructions, and make sure that the resulting
`open-wbo_static` binary is on your `$PATH`.

### Installing the project code

We recommend installing the current project *in development mode*, and *on a fresh Python3 virtual 
environment*. This will help debugging and troubleshooting. Assuming that `virtualenv` is installed in you
machine, create a new environment for `sltp` by running (usually in some user directory outside the project tree):

    mkdir -p ~/virtualenvs/sltp
    virtualenv -p python3 ~/virtualenvs/sltp
    
You should then be able to "enter" this virtualenv by running `source ~/virtualenvs/sltp/bin/activate`.
Our setup script requires `pip >= 18.1`. If `pip --version` shows that you have an older version, 
then you should be able to upgrade it with `python -m pip install -U pip` 

Finally, to perform the installation (from inside the virtual environment) 
issue the following command on the root directory of the project:

    pip install -e .

This should build and install all the required dependencies, including some the C++ feature generator
module, which should be compiled and installed transparently.


## Usage

Individual experiments are on the `experiments` folder, e.g. see `experiments/blocks.py` for an example.
A file such as `blocks.py` contains different experiment configurations for learning in the blocks domain.
We invoke the pipeline with <script-name> <experiment-name> <pipeline-steps-to-be-executed>,
where the last parameter is an optional  list of experiment step IDs (e.g.: 1 2 3). All steps are run if no IDs are specified.
Example invocations:

    # AAAI Experiment: BW with goal clear(x)
    ./run.py blocks:aaai_clear_x_simple_hybrid  # Theory T_G
    ./run.py blocks:aaai_clear_x_no_marking     # Theory T

    # AAAI Experiment: BW with goal on(x,y)
    ./run.py blocks:aaai_bw_on_x_y_completeness_opt  # Theory T_G
    ./run.py blocks:aaai_bw_on_x_y_completeness_opt_no_marking     # Theory T

    # AAAI Experiment: Gripper
    ./run.py gripper:aaai_prob01  # Theory T_G
    ./run.py gripper:aaai_prob01_no_marking     # Theory T

    # ...

The experiment names are not as informative as could be.
The configuration of each experiment can be inspected by looking at the experiment file.

## Using the SLTP Docker image 
We can also run experiments from within the Docker image. Assuming you want to run experiment `p1` from the `Visitall`
domain, leaving all intermediate files and results in a `workspace` directory in the host machine: 
    
    mkdir workspace
    docker pull gfrancesm/sltp
    docker run -it --mount src=`pwd`/workspace,target=/root/projects/workspace,type=bind \
        gfrancesm/sltp sltp gripper:aaai_prob01 --workspace /root/projects/workspace    

You can also run experiments that are not integrated within the standard SLTP command runner, just as you would run
any normal script within a Docker image. Suppose you have _in the host machine_ a SLTP experiment script such as the one 
found in the `examples` folder:
    
    /tmp/workspace$ ls
    gripper.py

You can run that script _within the Docker container_ as follows:

    /tmp/workspace$ docker pull gfrancesm/sltp
    /tmp/workspace$ docker run -it --mount src=`pwd`,target=/root/projects/workspace,type=bind \
                       gfrancesm/sltp /root/projects/workspace/gripper.py aaai_prob01 --workspace /root/projects/workspace
    
where `aaai_prob01` is the experiment name configured in your `gripper.py` script.

### Sampling
The SLTP pipeline uses a modified version of the FS planner (see installation instructions above)
in order to sample the (reachable part of the) state spaces of given PDDL instances. The sampling
steps of the pipeline write a few files with a textual representation of the state space. To see
an example, run e.g. from the `experiments` directory:

    $ ./run.py gripper:aaai_prob01 1 2 
    ================================================================================
    (pid: 17590) STARTING STEP #1: Sampling of the state space
    ================================================================================
    ...
    ================================================================================
    (pid: 17606) STARTING STEP #2: Generation of the training sample
    ================================================================================
    ...
    2020-02-23 20:21:18 INFO     Sample after resampling: roots: 2, states: 381, transitions: 812 (28 optimal), goals: 4, unsolvable: 0
    2020-02-23 20:21:18 INFO     Resampled states logged at "/home/frances/projects/code/sltp/workspace/2020022320.gripper.prob01_sample02.2000.cs-8/resampled.txt"
    2020-02-23 20:21:18 INFO     Printing transition matrix with 381 states and 812 transitions to '/home/frances/projects/code/sltp/workspace/2020022320.gripper.prob01_sample02.2000.cs-8/transition-matrix.dat'

The above command runs two different pipeline steps: one for calling `n` times the FS planner on `n`
given problem instances; the second step, reads all outputs, and consolidates them into a single set
of transition samples, assigning all states a unique ID, etc. The result of this consolidation is 
output, in the example above, to file
`/home/frances/projects/code/sltp/workspace/2020022320.gripper.prob01_sample02.2000.cs-8/resampled.txt`
The format of the file is relatively human-readable. Each state has an id such as #7, plus a number
of symbols that denote different (possibly complementary) properties of the state. These can be:

    "*": the state is a goal
    "^": the state has been fully expanded
    "º": the state is a dead-end
    "=": the state is a root of one of the sampled instances
    "+": the state lies on an optimal path to the goal of its instance

This sampling pipeline step also generates a second output file, in a somewhat more machine-readable
file, but also with less informatoin. In the example above, this file is 
`/home/frances/projects/code/sltp/workspace/2020022320.gripper.prob01_sample02.2000.cs-8/transition-matrix.dat` 
This file contains one line per every expanded state, such as e.g.:

    11 2 5 31

This line corresponds to state with id 11, and denotes that state 11 has children states 2, 5 and
31. 

### Generating Features
In order to be able to build the C++ feature generator, run `cmake . && make` on directory
`src/features`. As it stands now, this generator generates features up to a certain syntactic 
complexity bound and *based on a given set of transition samples* (i.e. sampled state space).
Ideally, one could also want to generate features independently of any transitions or state spaces,
but at the moment the generator creates only features that are non-redundant based on the denotation
over the whole set of sampled transitions, hence the requirement of having those transitions.
Two features are considered redundant if they have exactly the same denotation over all sampled
states.

As an example, the feature generation step for some set of Gripper can be run with  
 
    $ ./run.py gripper:aaai_prob01  1 2 3
    ================================================================================
    (pid: 18828) STARTING STEP #1: Sampling of the state space
    ================================================================================
    ...
    ...
    ================================================================================
    (pid: 18848) STARTING STEP #3: C++ feature generation module
    ================================================================================
    2020-02-23 20:54:08 INFO     Generating non-redundant concepts from sample set: roots: 2, states: 381, transitions: 812 (28 optimal), goals: 4, unsolvable: 0
    2020-02-23 20:54:08 INFO     Invoking C++ feature generation module
    ...
    DL::Factory: #concepts-final=60
    A total of 0 features were marked as goal-identifying
    FEATURES: #features=7, #nullary=0, #boolean=2, #numerical=5, #distance=0, #conditional=0
    ...
    2020-02-23 20:54:08 INFO     Reading feature information from /home/frances/projects/code/sltp/workspace/2020022320.gripper.prob01_sample02.2000.cs-8/feature-info.io
    ...
    2020-02-23 20:54:08 INFO     Reading denotation matrix from /home/frances/projects/code/sltp/workspace/2020022320.gripper.prob01_sample02.2000.cs-8/feature-matrix.io

The above command runs 3 different steps; the first two generate the sample (see above section);
the third is the one that actually generates all non-redundant concepts and features.
The relevant output of this step is, in the example above, file 
`/home/frances/projects/code/sltp/workspace/2020022320.gripper.prob01_sample02.2000.cs-8/feature-matrix.io`
This file is the denotation matrix of the generated features. The value at row i, column j
is the denotation of the j-th generated feature in state #i.
File `/home/frances/projects/code/sltp/workspace/2020022320.gripper.prob01_sample02.2000.cs-8/feature-info.io`
also contains the textual representation of each feature. 



### Software Requirements
SLTP has been tested on Python 3.6+ / Ubuntu

