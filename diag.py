from typing import Callable
from scipy import sparse
def get_preconditioner(name: str) -> Callable:
    """Return pre-conditioner based on name."""
    if name == 'jacobi':
        return jacobi
    raise KeyError('No pre-conditioner for provided name = {}'.format(name))


def jacobi(a_matrix):
    """Return csc_matrix obtained by inverse(diag(a_matrix))"""
    _to_return = {}  # type: dict
    if 'inverted' not in _to_return:
        # we want to calculate matrix inversion only once...
        diag = 1.0 / (a_matrix.diagonal() + 1e-6)
        pre_1 = sparse.csc_matrix(sparse.diags(diag))
        _to_return['inverted'] = pre_1
    return _to_return['inverted']
