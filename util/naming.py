import os


def compute_instance_tag(instance, num_states, **kwargs):
    inst = os.path.splitext(os.path.basename(instance))[0]
    dom = os.path.basename(os.path.dirname(instance))
    tag = "{}.{}.{}".format(dom, inst, num_states)
    return tag


def compute_experiment_tag(instance_tag, concept_depth, **kwargs):
    return "{}.cd-{}".format(instance_tag, concept_depth)


def compute_serialization_name(basedir, name):
    return os.path.join(basedir, '{}.json'.format(name))


def compute_maxsat_filename(config):
    return os.path.join(config.experiment_dir, 'encoding.cnf'.format())
