


# The SLTP Framework: Sample, Learn, Transform & Plan




## Installation

To install in development mode, issue the following command on the root directory of the project:

    pip install -e . --process-dependency-links

This should install all the required dependencies.
If you have a modern version of setuptools, you may see some deprecation warnings about the 
`--process-dependency-links` options, you can safely ignore them for the moment being.

## Usage

Individual experiments are on the `experiments` folder, e.g. see `experiments/blocks.py` for an example.
A file such as `blocks.py` contains different experiment configurations for learning in the blocks domain.
We invoke the pipeline with <script-name> <experiment-name> <pipeline-steps-to-be-executed>,
where the last parameter is usually `--all`, to run all steps in the pipeline.
Example invokations would include:


        # AAAI Experiment: BW with goal clear(x)
	./blocks.py aaai_clear_x_simple_hybrid --all  # Theory T_G
	./blocks.py aaai_clear_x_no_marking --all     # Theory T

        # AAAI Experiment: BW with goal on(x,y)
	./blocks.py aaai_bw_on_x_y_completeness_opt --all             # Theory T_G
	./blocks.py aaai_bw_on_x_y_completeness_opt_no_marking --all  # Theory T

        # AAAI Experiment: Gripper with goal on(x,y)
	./gripper.py aaai_prob01 --all             # Theory T_G
	./gripper.py aaai_prob01_no_marking --all  # Theory T

        # ...

The experiment names are not as informative as could be.
The configuration of each experiment can be inspected by looking at the experiment file.


### Software Requirements

* Python 3.6+ with the following dependencies:
  - pip
  - setuptools
  - python-dev
  - numpy

