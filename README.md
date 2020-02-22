


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
	python ./blocks.py aaai_clear_x_simple_hybrid  # Theory T_G
	python ./blocks.py aaai_clear_x_no_marking     # Theory T

    # AAAI Experiment: BW with goal on(x,y)
	python ./blocks.py aaai_bw_on_x_y_completeness_opt             # Theory T_G
	python ./blocks.py aaai_bw_on_x_y_completeness_opt_no_marking  # Theory T

    # AAAI Experiment: Gripper with goal on(x,y)
	python ./gripper.py aaai_prob01             # Theory T_G
	python ./gripper.py aaai_prob01_no_marking  # Theory T

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



### Software Requirements

* Python 3.6+ with the following dependencies:
  - pip
  - setuptools
  - python-dev
  - numpy

