import numpy as np
import itertools

def generate_basis(L, N):
    """
    Generate all Fock states |n0, n1, ..., n_{L-1}> with sum n_i = N.
    Returns:
        basis: list of tuples (n0, n1, ..., n_{L-1})
        index: dict mapping tuple -> integer index
    """
    basis = []
    # compositions of N into L parts
    for occ in itertools.product(range(N+1), repeat=L):
        if sum(occ) == N:
            basis.append(occ)
    index = {occ: i for i, occ in enumerate(basis)}
    return basis, index


def build_hamiltonian(L, N, t, U, mu):
    """
    Build the full Bose-Hubbard Hamiltonian in the fixed-N sector.
    """
    basis, index = generate_basis(L, N)
    dim = len(basis)
    H = np.zeros((dim, dim), dtype=float)

    # nearest-neighbor bonds on a chain
    bonds = [(i, i+1) for i in range(L-1)]

    for i, occ in enumerate(basis):
        occ = list(occ)

        # diagonal part
        E_diag = 0.0
        for s in range(L):
            n = occ[s]
            E_diag += 0.5 * U * n * (n - 1) - mu * n
        H[i, i] += E_diag

        # kinetic part
        for (a, b) in bonds:
            na, nb = occ[a], occ[b]

            # b_a^\dagger b_b
            if nb > 0:
                new = occ.copy()
                new[a] += 1
                new[b] -= 1
                j = index[tuple(new)]
                amp = -t * np.sqrt((na + 1) * nb)
                H[i, j] += amp
                H[j, i] += amp

            # b_b^\dagger b_a
            if na > 0:
                new = occ.copy()
                new[a] -= 1
                new[b] += 1
                j = index[tuple(new)]
                amp = -t * np.sqrt((nb + 1) * na)
                H[i, j] += amp
                H[j, i] += amp

    return H, basis


def bh_ed(L, N, t, U, mu, beta):
    """
    Compute diagonal, kinetic, and total energies at finite temperature.
    """
    H, basis = build_hamiltonian(L, N, t, U, mu)
    dim = H.shape[0]

    # diagonalize
    eigvals, eigvecs = np.linalg.eigh(H)

    # thermal weights
    w = np.exp(-beta * eigvals)
    Z = np.sum(w)

    # operators
    D_op = np.zeros(dim)
    K_op = H.copy()  # total H = diag + kinetic; we'll subtract diag later

    # build diagonal operator
    for i, occ in enumerate(basis):
        E_diag = 0.0
        for n in occ:
            E_diag += 0.5 * U * n * (n - 1) - mu * n
        D_op[i] = E_diag

    # kinetic operator = H - diag(D_op)
    for i in range(dim):
        K_op[i, i] -= D_op[i]

    # thermal averages
    E_diag = 0.0
    E_kin = 0.0

    for k in range(dim):
        psi = eigvecs[:, k]
        prob = psi * psi

        E_diag += w[k] * np.sum(prob * D_op)
        E_kin  += w[k] * (psi @ (K_op @ psi))

    E_diag /= Z
    E_kin  /= Z
    E_tot  = E_diag + E_kin

    return E_diag, E_kin, E_tot


# ------------------------------------------------------------
# Example usage
# ------------------------------------------------------------
if __name__ == "__main__":
    L = 2
    N = 2
    t = 1.0
    U = 4.0
    mu = 0.0
    beta = 4.0

    E_diag, E_kin, E_tot = bh_ed(L, N, t, U, mu, beta)

    print("=== Exact Diagonalization (L-site BH, fixed N) ===")
    print(f"L = {L}, N = {N}, t = {t}, U = {U}, mu = {mu}, beta = {beta}")
    print()
    print(f"E_diag = {E_diag}")
    print(f"E_kin  = {E_kin}")
    print(f"E_tot  = {E_tot}")
