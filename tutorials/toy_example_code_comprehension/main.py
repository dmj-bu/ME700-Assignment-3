import numpy as np
from solver import newton_solver
from equations import residual_fcn, jacobian_fcn

# Initial guess (explicitly clear)
u_initial = np.array([1.0, 1.0])

# Solve nonlinear system using Newton solver
u_solution = newton_solver(residual_fcn, jacobian_fcn, u_initial)

# Clearly print final solution
print(f"Solution found: u = {u_solution}")
