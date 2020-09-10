# bioinformatica

## Input

A string describing the system, with the next format:

```
# Example system

# Initial and final simulation times
initial_time = 0
final_time = 100

# System reactions
r1: A -> B; k1
r2: B -> C; k2

# Values of the reactions constants
k1 = 0.11
k2 = 0.1

# Initial species concentrations
A = 100
B = 0
C = 0
```

It's posible to write comments using `#`.

## Usage
```
from gillespie import run_simulation
run_simulation(sistema)
```

## Examples
The file `examples.py` contains three example systems. They can be run with the following command:

```
python3 examples.py
```

## References
* Daniel T. Gillespie, *A General Method for Numerically Simulating the Stochastic Time Evolution of Coupled Chemical Reactions* (1976): <http://web.mit.edu/endy/www/scraps/signal/JCompPhys(22)403.pdf>
* Daniel T. Gillespie, *Exact Stochastic Simulation of Coupled Chemical Reactions* (1977): <http://www.stat.yale.edu/~jtc5/GeneticNetworksGroup/Gillespie1977.pdf>
