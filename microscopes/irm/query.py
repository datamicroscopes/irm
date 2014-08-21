"""The query interface for IRM.

Note that the methods of this interface all take a list of latent state objects
(as opposed to a single latent).

"""

from microscopes.common import query


def zmatrix(domain, latents):
    """Compute a z-matrix (cluster co-assignment matrix). The ij-th entry of a
    z-matrix is a real value scalar between [0, 1] indicating the frequency of
    how often entities i and j appear in the same cluster.

    Parameters
    ----------
    domain : int
        The domain ID to compute the z-matrix for
    latents : list of mixturemodel latent objects
        The latents should all be points in the state space of the same
        structural model. The implementation currently does not check for this.

    Returns
    -------
    zmat : (N, N) ndarray

    Notes
    -----
    Currently does not support a sparse zmatrix representation, so only use
    this for small N.

    """
    return query.zmatrix([latent.assignments(domain) for latent in latents])
