import numpy as np


def newton_solver(residual_func, jacobian_func, u_init, tol=1e-8, max_iter=50):
    """
    Implements Newton's method explicitly for nonlinear systems.
    """
    u = u_init.copy()
    for iter in range(max_iter):
        res = residual_func(u)
        norm_res = np.linalg.norm(res)
        print(f"Iteration {iter}: Residual norm = {norm_res:.3e}")

        if norm_res < tol:
            print("Converged successfully!")
            return u

        J = jacobian_func(u)
        try:
            delta_u = np.linalg.solve(J, -res)
        except np.linalg.LinAlgError:
            raise RuntimeError("Jacobian is singular!")

        u += delta_u

    raise RuntimeError("Newton solver did not converge within max iterations")
