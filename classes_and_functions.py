import numpy as np
import scipy.sparse as ssp
from scipy.special import factorial

def envelope(t, Tend):
    """
    Can be done with np.piecewise as well.
    Input: 
        Time (real number): t
        time at which the puls ends (real number >0): T
    Returns:
        envelope function f(t) (real number)
    """
    assert(Tend > 0)
    if 0 <= t <= Tend: return np.sin(np.pi * t / Tend)**2
    else: return 0.0
venvelope = np.vectorize(envelope) # vectorize this function so it accepts numpy arrays

def potential(x, V_0 = 1.82, alpha = 0.44, beta = 1.0):
    """
    Implementation of a soft-core potential. 1 / x doesn't work in 1D. Therefore, one has to adapt.
    """
    return -alpha * V_0 / np.sqrt((beta * x)**2 + alpha**2)
vpotential = np.vectorize(potential) # vectorize this function so it accepts numpy arrays

def coefficients_equi(derivative, order = 2):
    """
    Calculation of the finite difference coefficients for an uniform grid
    """
    mat = np.full((2 * order + 1, 2 * order + 1), 1, dtype = np.float64)
    vec = np.zeros((2 * order + 1,), dtype = np.float64)
    dummy_list = np.arange(-order, order + 1, 1)
    for i in range(1, 2 * order + 1):
        mat[i] = dummy_list**i
    vec[derivative] = factorial(derivative)
    return np.linalg.solve(mat, vec)

def coeffs_from_difference_map(dxs, size, order):
    """
    Get the coefficients for the first and second derivative from a difference array
    """
    first_deriv_coefficients = np.empty((size, 2 * order + 1), dtype = np.float64)
    second_deriv_coefficients = np.empty((size, 2 * order + 1), dtype = np.float64)

    workspace = np.empty((2 * order + 1, 2 * order + 1), dtype = np.float64)
    workspace[0] = 1.0

    rhs_first_deriv     = np.zeros((2 * order + 1,), dtype = np.float64)
    rhs_first_deriv[1]  = 1.0
    rhs_second_deriv    = np.zeros((2 * order + 1,), dtype = np.float64)
    rhs_second_deriv[2] = 1.0
    for i_x in range(size):
        dx_pts = np.concatenate((-np.flip(np.cumsum(np.flip(dxs[i_x: i_x + order]))), [0], np.cumsum(dxs[i_x + order: i_x + 2 * order])))
        for i in range(1, 2 * order + 1): workspace[i] = np.power(dx_pts, i) / factorial(i)
        first_deriv_coefficients[i_x] = np.linalg.solve(workspace, rhs_first_deriv)
        second_deriv_coefficients[i_x] = np.linalg.solve(workspace, rhs_second_deriv)
    return first_deriv_coefficients, second_deriv_coefficients

def calculate_FDM_coefficients_DBC(xs, order = 2, return_dx = False):
    """
    Implementation of the finite difference method on a nonuniform grid for Dirichlet boundary condition.
    """
    assert(hasattr(xs, '__len__'))
    assert(isinstance(return_dx, bool))
    size = len(xs)
    dxs  = np.diff(xs)
    dx_outer_left  = np.flip(dxs[:order])
    dx_outer_right = np.flip(dxs[-order:])
    total_dxs = np.concatenate((dx_outer_left, dxs, dx_outer_right))

    first_deriv_coefficients, second_deriv_coefficients = coeffs_from_difference_map(total_dxs, size, order)
    
    if return_dx:
        return first_deriv_coefficients, second_deriv_coefficients, total_dxs
    return first_deriv_coefficients, second_deriv_coefficients

def calculate_FDM_coefficients_radial(rs, return_dx = False):
    """
    Implementation of the finite difference method on a nonuniform grid for radial grid of first order (results in tridiagonal matrix).
    """
    assert(hasattr(rs, '__len__'))
    assert(isinstance(return_dx, bool))
    size = len(rs)
    order = 1
    drs = np.diff(rs)
    total_drs = np.concatenate(([2 * rs[0]], drs, [drs[-1]]))
    
    first_deriv_coefficients, second_deriv_coefficients = coeffs_from_difference_map(total_drs, size, order)
    
    if return_dx:
        return first_deriv_coefficients, second_deriv_coefficients, total_drs
    return first_deriv_coefficients, second_deriv_coefficients

def build_jacobi_det_grid(dxs, size, order = 2):
    """
    Function to calculate the Jacobi determinant of the grid for the integration
    """
    xs_shifted = np.cumsum(np.concatenate(([0], dxs)))
    offsets = range(2 * order + 1)
    jac = ssp.diags(coefficients_equi(1, order = order), offsets = offsets, shape = (size, size + 2 * order), dtype = np.float64).dot(xs_shifted)
    return jac

def vector_to_sparse_mat(vec):
    return ssp.diags(vec)

def build_derivative_matrices(coeffs_table, size, order, neumann = False):
    offsets = range(-order, order + 1)
    if neumann:
        assert(order == 1)
        coeffs_table[0, 1] += coeffs_table[0, 0]
    return ssp.diags([coeffs_table[max(order - i, 0):, i] for i in range(2 * order + 1)], offsets = offsets, shape = (size, size), dtype = np.float64).tocsr()


class Hamiltonian:
    Tend = 1.0
    omega = None
    def __init__(self, xs, pot_func = vpotential, order = 2):
        """
        Input:
            Grid points (numpy array of real numbers (sorted)): xs
            Angular frequency of the pulse (real number): omega
            Function for the envelope of the pulse (optional python function accepting time and end of pulse): envelope_func(t, Tend)
            Order of the FDM algorithm (Integer number >0, 2*order+1 is the number of points used): order
        """
        # check if xs has the len attribute and can thus be transformed to an array (doesn't work with dictionaries)
        assert(hasattr(xs, '__len__'))
        self.xs = np.array(xs)
        # Check if xs is 1D
        assert(self.xs.ndim == 1)
        # Check if xs is sorted
        assert(np.all(self.xs[:-1] <= self.xs[1:]))
        self.size = len(self.xs)

        self.pot_func = pot_func
        coefficients1, coefficients2, dxs = calculate_FDM_coefficients_DBC(self.xs, order = order, return_dx = True)
        self.jac = build_jacobi_det_grid(dxs, self.size, order = order)

        self._build_static_parts([coefficients1, coefficients2], order = order)

    def _build_static_parts(self, FDM_coefficients, order = 2):
        """
        Implementation of the kinetic energy and the time independent potential term
        """
        first_deriv = FDM_coefficients[0]
        second_deriv = FDM_coefficients[1]
        offsets = range(-order, order + 1)
        self.momentum = -1.j * ssp.diags([first_deriv[max(order - i, 0):, i] for i in range(2 * order + 1)], offsets = offsets, shape = (self.size, self.size), dtype = np.complex128).tocsr()
        self.kin_energy = -0.5 * ssp.diags([second_deriv[max(order - i, 0):, i] for i in range(2 * order + 1)], offsets = offsets, shape = (self.size, self.size), dtype = np.float64).tocsr()
        pot_energy = ssp.diags(self.pot_func(self.xs), dtype = np.float64)
        self.H_0 = (self.kin_energy + pot_energy).tocsr()
        return

    def set_pulse_length(self, Tend):
        """
        Set the pulse length of the system
        """
        self.Tend = Tend
        return

    def minimal_coupling_phi(self, E_0, omega, envelope_func = venvelope):
        """
        Couple an external electric field with a time envelope function to the static Hamiltonian
        """
        self.minimal_coup_phi = ssp.diags(self.xs * E_0, dtype = np.float64).tocsr()
        self.omega = omega
        self.envelope_func = envelope_func
        return
    
    def build_at_time(self, t):
        """
        Build the Hamiltonian at a certain time
        """
        if self.omega is not None:
            H = self.H_0 + self.minimal_coup_phi * np.cos(self.omega * t) * self.envelope_func(t, self.Tend)
        else:
            H = self.H_0
        return H

