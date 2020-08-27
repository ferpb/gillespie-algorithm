# Implementación del algoritmo de Gillespie

import math
import numpy as np
from matplotlib import pyplot as plt
import random


# Entrada del algoritmo
# =====================

# Concentración inicial de cada una de las especies del sistema

# Tiempo inicial y final de la simulación

# Matriz de actualización de las especies en el sistema:
#    Las reacciones deben ser unidireccionales, así que si una reacción es
#    bidirección hay que considerarla como dos reacciones unidireccionales
#    separadas.

#    Ejemplo:
#
#      Reacciones:
#      A + B   ->  AB
#      AB      ->  A + B
#      AB      ->  A + C
#
#      Matriz:
#               A   B   AB  C
#      reacc0 [-1, -1,  1,  0]
#      reacc1 [ 1,  1, -1,  0]
#      reacc2 [ 1,  0, -1,  1]

# Funciones de velocidad para cada reacción en un array en el mismo orden
# que la matriz de actualización

update_matrix = np.array([[-1, -1,  1,  0],
                          [ 1,  1, -1,  0],
                          [ 1,  0, -1,  1]])

propensity_functions = np.array([1, 2, 3])
species_concentration = np.array([25, 583, 5, 13])

# Tiempo en milisegundos
initial_time = 0
final_time = 100


# Algoritmo
# =========

# Eliminar elementos positivos de la matriz de actualización
# para obtener una matriz de uso de reactivos
species_usage = update_matrix
species_usage[species_usage > 0] = 0
species_usage = -species_usage

# La probabilidad de elección de cada reacción es proporcional a la concentración
# de sus reactivos en la mezcla
def next_reaction(propensity_functions, species_concentration, species_usage):
    prop = species_usage * species_concentration
    prop = np.sum(prop, axis=1)
    prop = prop * propensity_functions
    return np.argmax(prop)


def update_concentrations(reaction, species_concentration, propensity_functions):
    return np.add(species_concentration,
                  propensity_functions[reaction])


def gillespie(initial_time, final_time,
              update_matrix, species_usage, propensity_functions, species_concentration):

    t = initial_time

    times = [t]
    concentrations = [species_concentration]

    # Iterar mientras que hay reactivos y no se ha pasado el tiempo de la simuación
    while t <= final_time and np.sum(species_concentration) > 0:
        reaction = next_reaction(propensity_functions, species_concentration, species_usage)
        species_concentration = update_concentrations(reaction, species_concentration, propensity_functions)

        # incrementar el tiempo de forma aleatoria y siguiendo una
        # distribución uniforme
        t += random.random()

        times.append(t)
        concentrations.append(species_concentration)


    # Plot species concentration
    concentrations = np.array(concentrations)
    concentrations = concentrations.transpose()

    for c in concentrations:
        plt.scatter(times, c)
    plt.show()


gillespie(initial_time, final_time,
          update_matrix, species_usage, propensity_functions, species_concentration)
