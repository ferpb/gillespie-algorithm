# Implementación del algoritmo de Gillespie

import math
import numpy as np
from matplotlib import pyplot as plt
import random
import pdb


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


def calculate_a(propensity_functions, species_concentration, species_usage):
    # Calculamos las a_i
    a = species_usage * species_concentration
    # No se si hay que sumar multiplicar
    # a = np.sum(a, axis=1)
    a[a == 0] = 1
    a = np.prod(a, axis=1)
    a = a * propensity_functions

    a0 = np.sum(a)

    return a, a0


# La probabilidad de elección de cada reacción es proporcional a la concentración
# de sus reactivos en la mezcla
def next_reaction(propensity_functions, species_concentration, species_usage):
    # # Calculamos las a_i
    # a = species_usage * species_concentration
    # # No se si hay que sumar multiplicar
    # # a = np.sum(a, axis=1)
    # a[a == 0] = 1
    # a = np.prod(a, axis=1)
    # a = a * propensity_functions

    # a0 = np.sum(a)

    a, a0 = calculate_a(propensity_functions, species_concentration, species_usage)
    r2 = random.random()

    # pdb.set_trace()

    # Calculamos mu
    mu = 1
    while (np.sum(a[:mu]) < r2*a0) and (np.sum(a[:mu + 1]) >= r2*a0):
        mu += 1

    print("Encontrado mu", mu - 1)
    # print(species_concentration)

    return mu - 1


def next_time(propensity_functions, species_concentration, species_usage):
    # a = species_usage * species_concentration
    # a[a == 0] = 1
    # a = np.prod(a, axis=1)
    # a = a * propensity_functions

    _, a0 = calculate_a(propensity_functions, species_concentration, species_usage)
    r1 = random.random()

    # Calculamos tau
    tau = (1 / a0) * math.log(1 / r1)

    return tau


def update_concentrations(reaction, species_concentration, update_matrix):
    return np.add(species_concentration,
                  update_matrix[reaction])


def gillespie(initial_time, final_time,
              update_matrix, species_usage, propensity_functions, species_concentration):

    t = initial_time

    times = [t]
    concentrations = [species_concentration]

    # Iterar mientras que hay reactivos y no se ha pasado el tiempo de la simuación
    while t <= final_time and np.sum(species_concentration) > 0:
        reaction = next_reaction(propensity_functions, species_concentration, species_usage)
        species_concentration = update_concentrations(reaction, species_concentration, update_matrix)

        # incrementar el tiempo de forma aleatoria y siguiendo una
        # distribución uniforme
        t += next_time(propensity_functions, species_concentration, species_usage)

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
