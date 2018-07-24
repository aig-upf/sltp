#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import os

from tarski.fstrips import create_fstrips_problem, language
from tarski.io import FstripsWriter
from tarski.theories import Theory
from tarski.syntax import land

_CURRENT_DIR_ = os.path.dirname(os.path.realpath(__file__))


def main():
    lang = language(theories=[Theory.EQUALITY, Theory.ARITHMETIC])
    problem = create_fstrips_problem(domain_name='gridworld',
                                     problem_name='10x10',
                                     language=lang)

    coord_t = lang.interval('coordinate', lang.Integer, 1, 10)
    xpos = lang.function('X', coord_t)
    ypos = lang.function('Y', coord_t)

    # Create the actions
    problem.action(name='move-up', parameters=[],
                   precondition=ypos() < coord_t.upper_bound,
                   effects=[ypos() << ypos() + 1])

    problem.action(name='move-down', parameters=[],
                   precondition=ypos() > coord_t.lower_bound,
                   effects=[ypos() << ypos() - 1])

    problem.action(name='move-left', parameters=[],
                   precondition=xpos() > coord_t.lower_bound,
                   effects=[xpos() << xpos() - 1])

    problem.action(name='move-right', parameters=[],
                   precondition=xpos() < coord_t.upper_bound,
                   effects=[xpos() << xpos() + 1])

    problem.init.set(xpos, (), 1)
    problem.init.set(ypos, (), 10)

    problem.goal = land(xpos() == 2, ypos() == 3)

    writer = FstripsWriter(problem)
    writer.write(domain_filename=os.path.join(_CURRENT_DIR_, "domain.pddl"),
                 instance_filename=os.path.join(_CURRENT_DIR_, "instance.pddl"))


if __name__ == "__main__":
    main()
