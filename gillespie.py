"""Implementation of the Gillespie Algorithm"""

import math
import numpy as np
from matplotlib import pyplot as plt
import random
import pdb

from parser import gillespie_parse


def calculate_a(propensity_functions, species_concentration, species_usage):
    # Calculate a_i
    a = species_usage * species_concentration
    # Multiply the elements of each row that are bigger than 0
    am = np.ma.MaskedArray(a, mask=(species_usage <= 0))
    a = am.prod(axis=1)
    a = a.filled(0)

    a = a * propensity_functions

    # Calculate a_0 (sum of all a_i)
    a0 = np.sum(a)

    return a, a0


def next_reaction(a, a0):
    # The probability of choosing a reaction is proportional to the concentration
    # of its reactants in the mixture
    r2 = random.random()
    mu = 1
    while (mu < len(a)) and (np.sum(a[:mu]) < r2*a0):
        mu += 1
    return mu - 1


def next_time(a0):
    # The next time leap is proportional to the number of reactants in the mixture
    r1 = random.random()
    tau = (1 / a0) * math.log(1 / r1)
    return tau


def update_concentrations(reaction, species_concentration, update_matrix):
    return np.add(species_concentration,
                  update_matrix[reaction])


def isBlocked(species_concentration, species_usage):
    # If there not enough reactants in the mixture for any reaction to occur,
    # the system is blocked

    # Detect with reactans are not enough for earch reaction
    insufficient = np.logical_and(species_concentration < species_usage, species_usage > 0)

    # A reaction is blocked if any of its reactants is not sufficient
    blocked_reaction = np.any(insufficient, axis=1)

    # The system is blocked if all its reactions are blocked
    return np.all(blocked_reaction)


def gillespie(initial_time, final_time,
              update_matrix, species_names, propensity_functions, species_concentration):
    """
    Simulate the evolution of the system using the Gillespie Algorithm and plot the evolution of the
    species concentration over the simulation time.

    Args:
        initial_time: initial time of the simulation.

        final_time: final time of the simulation.

        update_matrix: num_reactions x num_species matrix with the species consumption
                       and generation of each reaction. Bidirectional reaction must be considered
                       as two separated unidirectional ones.

        species_names: Array of legth num_species with the names of the species
                       (in the same order as the update_matrix)

        propensity_functions: Array of length num_reactions with the propensity function value
                              of each reaction (in the same order as the update_matrix).

        species_concentration: Array of length num_species with the initial concentrations of
                               each specie in the system (in the same order as the update_matrix).

    Example:
        System:
            Initial time = 0
            Final time = 100

            A -> B; k1
            B -> C; k2

            k1 = 0.11
            k2 = 0.1

            A = 100
            B = 0
            C = 0

        Args:
            initial_time = 0
            final_time = 100
            update_matrix: [[-1,  1,  0],
                            [ 0, -1,  1]]
            species_names: ["A", "B", "C"]
            propensity_functions: [0.11, 0.1]
            species_concentration: [100, 0, 0]
    """
    species_usage = update_matrix.copy()
    species_usage[species_usage > 0] = 0
    species_usage = -species_usage

    t = initial_time

    times = [t]
    concentrations = [species_concentration]

    # Iterate while the simulation time has not been consumed
    while t < final_time:
        if isBlocked(species_concentration, species_usage):
            print("The system is blocked!")
            break

        a, a0 = calculate_a(propensity_functions, species_concentration, species_usage)
        reaction = next_reaction(a, a0)
        t += next_time(a0)

        species_concentration = update_concentrations(reaction, species_concentration, update_matrix)

        times.append(t)
        concentrations.append(species_concentration)


    print("t:", t)
    print("species_concentration:", species_concentration)

    # Plot species concentration
    concentrations = np.array(concentrations)
    concentrations = concentrations.transpose()

    for c in concentrations:
        # plt.scatter(times, c, s=0.2)
        plt.plot(times, c)
    plt.xlabel('Tiempo (ms)')
    plt.ylabel('Número de moléculas')
    plt.legend(species_names, markerscale=10)
    plt.show()


def run_simulation(data):
    """
    Parse the description of a chemical system and simulate it using the gillespie() function.

    Args:
        data: string containing the description of the system
    """
    (initial_time, final_time, update_matrix, species_names, propensity_functions, species_concentration) = gillespie_parse(data)

    gillespie(initial_time, final_time,
              update_matrix, species_names, propensity_functions, species_concentration)
