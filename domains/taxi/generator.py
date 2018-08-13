#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import itertools
import os
import random

from tarski.fstrips import create_fstrips_problem, language, DelEffect, AddEffect
from tarski.io import FstripsWriter
from tarski.theories import Theory
from tarski.syntax import forall, implies, exists, land, Tautology

_CURRENT_DIR_ = os.path.dirname(os.path.realpath(__file__))


def create_noop(problem):
    # A hackish no-op, to prevent the planner from detecting that the action is useless and pruning it
    lang = problem.language
    cell_t, loct = lang.get("cell", "loct")
    x = lang.variable("x", cell_t)
    problem.action(name='noop', parameters=[x], precondition=(loct() == x), effects=[loct() << x])


def create_actions(problem, add_fuel):
    lang = problem.language
    cell_t = lang.get("cell")
    loc_taxi, loc_passenger, inside, adjacent = lang.get("loct", "locp", "inside_taxi", "adjacent")

    problem.action(name='pick-passenger', parameters=[],
                   precondition=loc_passenger() == loc_taxi(),
                   effects=[loc_passenger() << inside])

    problem.action(name='drop-passenger', parameters=[],
                   precondition=loc_passenger() == inside,
                   effects=[loc_passenger() << loc_taxi()])

    to = lang.variable("to", cell_t)

    if add_fuel:  # Add the refill action plus a fuel-aware move action

        fuel_level_t, loc_fuel, current_fuel, max_fuel_level, succ = \
            lang.get("fuel_level", "loc_fuel", "current_fuel", "max_fuel_level", "succ")

        problem.action(name='refill', parameters=[],
                       precondition=loc_taxi() == loc_fuel(),
                       effects=[current_fuel() << max_fuel_level()])

        x0 = lang.variable("x0", fuel_level_t)
        x1 = lang.variable("x1", fuel_level_t)
        problem.action(name='move', parameters=[to, x0, x1],
                       precondition=(adjacent(loc_taxi(), to)) & (succ(x0, x1)) & (current_fuel() == x1),
                       effects=[loc_taxi() << to,
                                current_fuel() << x0])

    else:  # Add a fuel-independent move action
        problem.action(name='move', parameters=[to],
                       precondition=(adjacent(loc_taxi(), to)),
                       effects=[loc_taxi() << to])


def generate_domain(gridsize, add_noop=False, add_fuel=True):
    lang = language(theories=[Theory.EQUALITY, Theory.ARITHMETIC])
    problem = create_fstrips_problem(domain_name='taxi-fs',
                                     problem_name="taxi-fs-{}x{}".format(gridsize, gridsize),
                                     language=lang)

    cell_t = lang.sort('cell')

    loc_taxi = lang.function('loct', cell_t)
    loc_passenger = lang.function('locp', cell_t)

    adjacent = lang.predicate('adjacent', cell_t, cell_t)

    # Create the actions
    it = lang.constant("inside_taxi", cell_t)
    create_actions(problem, add_fuel)
    if add_noop:
        create_noop(problem)

    rng = range(0, gridsize)
    coordinates = list(itertools.product(rng, rng))

    def cell_name(x, y):
        return "c_{}_{}".format(x, y)

    # Declare the constants:
    coord_objects = [lang.constant(cell_name(x, y), cell_t) for x, y in coordinates]

    # Declare the adjacencies:
    adjacent_coords = [(a, b, c, d) for (a, b), (c, d) in itertools.combinations(coordinates, 2)
                       if abs(a-c) + abs(b-d) == 1]

    for a, b, c, d in adjacent_coords:
        problem.init.add(adjacent, cell_name(a, b), cell_name(c, d))
        problem.init.add(adjacent, cell_name(c, d), cell_name(a, b))

    cd = coord_objects[:]
    random.shuffle(cd)

    # Initial positions
    problem.init.set(loc_taxi, cd.pop())
    problem.init.set(loc_passenger, cd.pop())

    # Set the problem goal
    problem.goal = loc_passenger() == cd.pop()

    if add_fuel:
        # Our approach is not yet too int-friendly :-(
        # fuel_level_t = lang.interval('fuel_level', lang.Integer, lower_bound=0, upper_bound=10)
        fuel_level_t = lang.sort('fuel_level')

        current_fuel = lang.function('current_fuel', fuel_level_t)
        loc_fuel = lang.function('loc_fuel', cell_t)
        max_fuel_level = lang.function('max_fuel_level', fuel_level_t)
        min_fuel_level = lang.function('min_fuel_level', fuel_level_t)

        # The whole succ-predicate stuff
        succ = lang.predicate("succ", fuel_level_t, fuel_level_t)
        levels = ["f{}".format(i) for i in range(0, 11)]
        _ = [lang.constant(c, fuel_level_t) for c in levels]  # Create the "integer" objects
        _ = [problem.init.add(succ, x, y) for x, y in zip(levels, levels[1:])]

        problem.init.set(current_fuel, random.choice(levels))
        problem.init.set(min_fuel_level, levels[0])
        problem.init.set(max_fuel_level, levels[-1])
        problem.init.set(loc_fuel, cd.pop())


    return problem, [it]


def main():

    for gridsize in [3, 5]:
        problem, constants = generate_domain(gridsize, add_noop=True, add_fuel=False)
        writer = FstripsWriter(problem)
        writer.write(domain_filename=os.path.join(_CURRENT_DIR_, "domain.pddl"),  # We can overwrite the domain
                     instance_filename=os.path.join(_CURRENT_DIR_, "instance_{}.pddl".format(gridsize)),
                     domain_constants=constants)


if __name__ == "__main__":
    main()
