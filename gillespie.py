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

species_names = ["A", "B", "AB", "C"]

propensity_functions = np.array([1, 2, 3])
species_concentration = np.array([25, 583, 5, 13])

# Tiempo en milisegundos
initial_time = 0
final_time = 100



# Algoritmo
# =========

species_usage = update_matrix.copy()
species_usage[species_usage > 0] = 0
species_usage = -species_usage


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
    # Calculamos mu
    mu = 1
    while (mu < len(a)) and (np.sum(a[:mu]) < r2*a0):
        mu += 1
    # print("Encontrado mu", mu - 1)
    # print(species_concentration)
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
              update_matrix, species_usage, propensity_functions, species_concentration):

    t = initial_time

    times = [t]
    concentrations = [species_concentration]

    # Iterar mientras que hay reactivos y no se ha pasado el tiempo de la simuación
    while t <= final_time:
        if isBlocked(species_concentration, species_usage):
            print("The system is blocked!")
            print("t:", t)
            print("species_concentration:", species_concentration)
            break

        a, a0 = calculate_a(propensity_functions, species_concentration, species_usage)
        reaction = next_reaction(a, a0)
        t += next_time(a0)

        species_concentration = update_concentrations(reaction, species_concentration, update_matrix)

        times.append(t)
        concentrations.append(species_concentration)


    # Plot species concentration
    concentrations = np.array(concentrations)
    concentrations = concentrations.transpose()

    for c in concentrations:
        plt.scatter(times, c, s=0.2)
    plt.xlabel('Tiempo (ms)')
    plt.ylabel('Concentración de especies')
    plt.legend(species_names, markerscale=10)
    plt.show()


# pdb.set_trace()
gillespie(initial_time, final_time,
          update_matrix, species_usage, propensity_functions, species_concentration)
