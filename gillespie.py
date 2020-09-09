# Implementación del algoritmo de Gillespie

import math
import numpy as np
from matplotlib import pyplot as plt
import random
import pdb

from parser import gillespie_parse


# Entrada del algoritmo
# =====================

# data = """
# initial_time = 0
# final_time = 100

# r1: A + B -> AB; k1
# r2: 2AB -> C; k2

# k1 = 1.0

# A = 100
# C = 1
# B = 10
# AB = 0
# """

data = """
# Ejemplo de sistema

initial_time = 0
final_time = 100

r1: A + B -> AB; k1
r2: AB -> A + B; k2
r3: AB -> A + C; k3

k1 = 1
k2 = 2
k3 = 3

A = 25
B = 583
AB = 5
C = 13
"""


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

# 1
# update_matrix = np.array([[-1, -1,  1,  0],
#                           [ 1,  1, -1,  0],
#                           [ 1,  0, -1,  1]])

# species_names = ["A", "B", "AB", "C"]

# propensity_functions = np.array([5, 40, 3])
# species_concentration = np.array([25, 583, 5, 13])

# 2
# update_matrix = np.array([[-1, -1,  1],
#                           [1, 1,  -1]])
# species_names = ["A", "B", "AB"]
# propensity_functions = np.array([2, 1])
# species_concentration = np.array([10, 10, 0])

# 3
# update_matrix = np.array([[-1, 1, 0],
#                          [0, -1, 1]])
# species_names = ["A", "B", "C"]
# propensity_functions = np.array([0.11, 0.1])
# species_concentration = np.array([100, 0, 0])

# Tiempo en milisegundos
# initial_time = 0
# final_time = 100



# Algoritmo
# =========

def calculate_a(propensity_functions, species_concentration, species_usage):
    # Calculamos las a_i
    a = species_usage * species_concentration

    # Multiplicar los elementos de cada fila mayores que 0
    am = np.ma.MaskedArray(a, mask=(species_usage <= 0))
    a = am.prod(axis=1)
    a = a.filled(0)

    a = a * propensity_functions

    # Calculamos a_0 (suma de todas las a_i)
    a0 = np.sum(a)

    return a, a0


# La probabilidad de elección de cada reacción es proporcional a la concentración
# de sus reactivos en la mezcla
def next_reaction(a, a0):
    r2 = random.random()
    mu = 1
    while (mu < len(a)) and (np.sum(a[:mu]) < r2*a0):
        mu += 1
    return mu - 1


def next_time(a0):
    r1 = random.random()
    tau = (1 / a0) * math.log(1 / r1)
    return tau


def update_concentrations(reaction, species_concentration, update_matrix):
    return np.add(species_concentration,
                  update_matrix[reaction])


# Si no hay reactivos suficientes para que se produzca ninguna reacción,
# el sistema está bloqueado
def isBlocked(species_concentration, species_usage):
    # Detectar qué reactivos no son suficientes para llevar a cabo la reacción
    insufficient = np.logical_and(species_concentration < species_usage, species_usage > 0)

    # Una reacción está bloqueada si alguno de sus reactivos es insuficiente
    blocked_reaction = np.any(insufficient, axis=1)

    # El sistema está bloqueado si todas las reacciones están bloqueadas
    return np.all(blocked_reaction)


def gillespie(initial_time, final_time,
              update_matrix, species_names, species_usage, propensity_functions, species_concentration):

    t = initial_time

    times = [t]
    concentrations = [species_concentration]

    # Iterar mientras que hay reactivos y no se ha pasado el tiempo de la simuación
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


# pdb.set_trace()
(initial_time, final_time, update_matrix, species_names, propensity_functions, species_concentration) = gillespie_parse(data)
species_usage = update_matrix.copy()
species_usage[species_usage > 0] = 0
species_usage = -species_usage

print()
print(initial_time)
print(final_time)
print(update_matrix)
print(species_names)
print(propensity_functions)
print(species_concentration)
print(species_usage)

gillespie(initial_time, final_time,
          update_matrix, species_names, species_usage, propensity_functions, species_concentration)
