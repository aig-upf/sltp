import os


def compute_instance_tag(instance, num_states, **kwargs):
    inst = os.path.splitext(os.path.basename(instance))[0]
    dom = os.path.basename(os.path.dirname(instance))
    tag = "{}.{}.{}".format(dom, inst, num_states)
    return tag


def compute_experiment_tag(instance_tag, max_concept_size, **kwargs):
    return "{}.cs-{}".format(instance_tag, max_concept_size)


def compute_serialization_name(basedir, name):
    return os.path.join(basedir, '{}.pickle'.format(name))


def compute_maxsat_filename(config):
    return compute_info_filename(config, "maxsat_encoding.cnf")


def compute_maxsat_variables_filename(config):
    return compute_info_filename(config, "maxsat_variables.txt")


def compute_info_filename(config, name):
    return os.path.join(config["experiment_dir"], name)
