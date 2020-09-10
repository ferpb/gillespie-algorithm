from gillespie import run_simulation

# Example 1
run_simulation(
"""
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
)

# Example 2
run_simulation(
"""
initial_time = 0
final_time = 100

r1: A + B -> AB; kb
r2: AB -> A + B; kd

kb = 1
kd = 2

A = 20
B = 10
AB = 0
"""
)

# Example 3
run_simulation(
"""
initial_time = 0
final_time = 100

r1: A -> B; k1
r2: B -> C; k2

k1 = 0.11
k2 = 0.1

A = 100
B = 0
C = 0
"""
)
