import numpy as np


# ------------------------------------------------------------
# Generate all Fock states with sum n_i = N
# ------------------------------------------------------------
def generate_fock_basis(L, N):
    def compositions(n, k):
        if k == 1:
            yield (n,)
        else:
            for i in range(n + 1):
                for rest in compositions(n - i, k - 1):
                    yield (i,) + rest
    return list(compositions(N, L))


# ------------------------------------------------------------
# Build Bose–Hubbard Hamiltonian
# ------------------------------------------------------------
def build_bose_hubbard_hamiltonian(L, N, t, U, mu=0.0, periodic=False):
    basis = generate_fock_basis(L, N)
    dim = len(basis)
    index = {state: i for i, state in enumerate(basis)}

    H = np.zeros((dim, dim), dtype=np.complex128)

    # Diagonal part
    for i, state in enumerate(basis):
        n = np.array(state, dtype=int)
        diag = 0.5 * U * np.sum(n * (n - 1)) - mu * np.sum(n)
        H[i, i] = diag

    # Hopping part
    def neighbors(j):
        if periodic:
            return [(j - 1) % L, (j + 1) % L]
        else:
            out = []
            if j - 1 >= 0:
                out.append(j - 1)
            if j + 1 < L:
                out.append(j + 1)
            return out

    for i, state in enumerate(basis):
        n = list(state)
        for j in range(L):
            for k in neighbors(j):
                if n[k] > 0:
                    new_n = n.copy()
                    new_n[k] -= 1
                    new_n[j] += 1
                    new_state = tuple(new_n)
                    col = index[new_state]
                    amp = -t * np.sqrt((n[j] + 1) * n[k])
                    H[i, col] += amp

    H = 0.5 * (H + H.conj().T)
    return H, basis


# ------------------------------------------------------------
# Diagonalize
# ------------------------------------------------------------
def diagonalize_bose_hubbard(L, N, t, U, mu=0.0, periodic=False):
    H, basis = build_bose_hubbard_hamiltonian(L, N, t, U, mu, periodic)
    evals, evecs = np.linalg.eigh(H)
    return evals, evecs, basis, H


# ------------------------------------------------------------
# Entanglement entropies
# ------------------------------------------------------------
def entanglement_entropies(psi, basis, L, L_A):
    configs_A = {}
    configs_B = {}
    list_A = []
    list_B = []

    for state in basis:
        a = state[:L_A]
        b = state[L_A:]
        if a not in configs_A:
            configs_A[a] = len(list_A)
            list_A.append(a)
        if b not in configs_B:
            configs_B[b] = len(list_B)
            list_B.append(b)

    dim_A = len(list_A)
    dim_B = len(list_B)

    Psi = np.zeros((dim_A, dim_B), dtype=np.complex128)

    for idx, state in enumerate(basis):
        a = configs_A[state[:L_A]]
        b = configs_B[state[L_A:]]
        Psi[a, b] += psi[idx]

    rho_A = Psi @ Psi.conj().T
    rho_A = 0.5 * (rho_A + rho_A.conj().T)
    rho_A /= np.trace(rho_A).real

    evals = np.linalg.eigvalsh(rho_A)
    evals = np.clip(evals.real, 0, 1)

    mask = evals > 0
    S_vN = -np.sum(evals[mask] * np.log(evals[mask]))
    S_2 = -np.log(np.sum(evals**2))

    return S_vN, S_2


# ------------------------------------------------------------
# Expectation values: diagonal and kinetic energies
# ------------------------------------------------------------
def energy_components(psi, H, basis, L, N, t, U, mu=0.0, periodic=False):
    dim = len(basis)
    psi = psi.reshape(dim)

    # Diagonal operator
    E_diag = 0.0
    for i, state in enumerate(basis):
        n = np.array(state)
        diag = 0.5 * U * np.sum(n * (n - 1)) - mu * np.sum(n)
        E_diag += np.conj(psi[i]) * diag * psi[i]
    E_diag = E_diag.real

    # Total energy
    E_tot = np.vdot(psi, H @ psi).real

    # Kinetic = total - diagonal
    E_kin = E_tot - E_diag

    return E_tot, E_diag, E_kin


# ------------------------------------------------------------
# Example usage
# ------------------------------------------------------------
if __name__ == "__main__":
    L = 2
    N = 2
    t = 0.5
    U = 1.0
    mu = 0.0

    evals, evecs, basis, H = diagonalize_bose_hubbard(L, N, t, U, mu)
    psi0 = evecs[:, 0]

    print("Ground state energy:", evals[0])

    E_tot, E_diag, E_kin = energy_components(psi0, H, basis, L, N, t, U, mu)
    print("Diagonal energy:", E_diag)
    print("Kinetic energy:", E_kin)

    S_vN, S_2 = entanglement_entropies(psi0, basis, L, L_A=1)
    print("S_vN =", S_vN)
    print("S_2  =", S_2)
