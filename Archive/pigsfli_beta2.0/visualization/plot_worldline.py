import matplotlib.pyplot as plt

def plot_worldline(worldline, M, beta):
    """
    worldline: list of lists
        worldline[site] = list of (tau, n) sorted by tau
    M: number of sites
    beta: imaginary time extent
    """

    fig, ax = plt.subplots(figsize=(10, 6))

    for site in range(M):
        kinks = worldline[site]

        # If no kinks, skip
        if len(kinks) == 0:
            continue

        # Sort by tau
        kinks = sorted(kinks, key=lambda x: x[0])

        # Draw segments between kinks
        for i in range(len(kinks) - 1):
            tau1, n1 = kinks[i]
            tau2, n2 = kinks[i + 1]

            # Horizontal segment: occupation constant
            ax.plot([site, site], [tau1, tau2], color='black')

            # Vertical jump if occupation changes
            if n1 != n2:
                ax.plot([site - 0.2, site + 0.2], [tau2, tau2], color='red')

        # Last segment to beta
        tau_last, n_last = kinks[-1]
        ax.plot([site, site], [tau_last, beta], color='black')

    ax.set_xlabel("Site index")
    ax.set_ylabel("Imaginary time τ")
    ax.set_title("Continuous-Time Worldline Configuration")
    ax.set_ylim(0, beta)
    ax.invert_yaxis()  # τ=0 at top, τ=β at bottom
    plt.show()

    # --------------------------------------------------------------------------
    
import numpy as np

data = np.loadtxt("worldline.txt")

# data columns: site, tau, n
sites = data[:,0].astype(int)
taus  = data[:,1]
ns    = data[:,2].astype(int)

M = sites.max() + 1
beta = taus.max()

# Build worldline structure
worldline = [[] for _ in range(M)]
for site, tau, n in zip(sites, taus, ns):
    worldline[site].append((tau, n))

plot_worldline(worldline, M, beta)

