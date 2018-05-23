import logging
import os
import argparse


def parse_arguments(args):
    parser = argparse.ArgumentParser(description="Learn generalized features and concepts")
    parser.add_argument('-k', help='Number of iterations to derive concepts and roles', action='store', default=0,
                        type=int)
    parser.add_argument('transitions', help='Name of file containing transitions (output from planner)')
    parser.add_argument('-d', '--domain', required=True, help='The PDDL domain filename')
    parser.add_argument('--debug', action='store_true', help='Print additional debugging information')
    return parser.parse_args(args)


def configure_logging(args):
    level = logging.DEBUG if args.debug else logging.INFO
    filename = os.path.basename(args.transitions)
    args.result_filename = '.'.join(filename.split('.')[:-1]) + ".{}it".format(args.k)
    filename = os.path.join('logs', args.result_filename + '.log')
    logging.basicConfig(filename=filename, filemode='w', level=level)


def bootstrap(arguments):
    args = parse_arguments(arguments)
    configure_logging(args)
    return args


