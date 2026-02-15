# Path-Integral Monte Carlo Engine  
## Architecture & API Documentation  
### Version: Foundation Layer (Initialization + Worldline Structure)

---

# 1. Overview

This document describes the architecture and functionality of the current PIMC (Path-Integral Monte Carlo) engine. It covers:

- System and lattice representation  
- Worldline data structures  
- Replica and configuration management  
- Initialization logic  
- Monte Carlo framework scaffolding  
- Testing and validation  

This documentation reflects the state of the code **after the foundational layer is complete and validated**.

---

# 2. Utilities

## 2.1 RNG

**Class:** `pimc::RNG`  
A thin wrapper around `std::mt19937_64` providing:

- `double uniform()` → random number in \([0,1)\)  
- `int randint(a,b)` → integer in \([a,b]\) inclusive  

Used throughout the simulation for stochastic decisions.

---

## 2.2 random_fock_state(M, N, rng)

Generates a random Fock state with:

- `M` lattice sites  
- `N` total bosons  
- multiple occupancy allowed  

Returns a `std::vector<int>` of size `M`.

Used to initialize worldlines before any Monte Carlo updates.

---

# 3. Static Model

## 3.1 System

Represents the **static physical parameters** of the Bose–Hubbard model:

- `L` — linear lattice size  
- `D` — dimension  
- `M = L^D` — total sites  
- `t` — hopping amplitude  
- `U` — onsite interaction  
- `mu` — chemical potential  

This class contains **no dynamic state**. It is read-only after construction.

---

## 3.2 Lattice

Builds a **periodic hypercubic lattice** of size \(L^D\).

Provides:

- adjacency lists for each site  
- dimension and size information  

Used by Monte Carlo moves to determine allowed hops.

---

# 4. Worldline Data Structures

## 4.1 Kink

A kink represents a **continuous-time event** on a worldline.

Fields:

| Field | Meaning |
|-------|---------|
| `tau` | imaginary time of event |
| `n` | occupation after event |
| `src`, `dest` | source/destination sites (for hops) |
| `prev`, `next` | linked-list pointers on same site |
| `src_replica`, `dest_replica` | replica indices |
| `partner` | index of paired kink (for hops) |

This is the fundamental building block of the worldline representation.

---

## 4.2 Worldline

Manages all kinks for a single replica.

### Responsibilities:

- Maintain **one linked list per site**  
- Store all kinks in a dynamic vector  
- Insert/remove kinks safely  
- Maintain partner links for hops  
- Provide consistency checks  

### Important methods:

- `addKink(k)`  
- `insertHop(k1, k2)`  
- `removeKink(idx)`  
- `deleteHop(idx)`  
- `checkConsistency()`  

This class contains **all dynamic state** for a replica.

---

## 4.3 initialize_worldline_from_fock

Creates a “flat” worldline:

- one kink per site  
- tau = 0  
- occupation = fock_state[site]  
- no hops, no partners  

This is the starting point before any Monte Carlo updates.

---

# 5. Replica and Configuration

## 5.1 Replica

A thin wrapper around a single `Worldline`.

Why it exists:

- future algorithms (replica coupling, SWAP moves, entanglement estimators) operate on multiple replicas  
- keeps `Configuration` clean and extensible  

---

## 5.2 Configuration

The **full dynamic state** of the simulation.

Contains:

- `System` (static physics)  
- `Lattice` (static geometry)  
- `std::vector<Replica>` (dynamic worldlines)  

Responsibilities:

- initialize replicas  
- provide access to system, lattice, and worldlines  
- serve as the object passed to Monte Carlo moves  

---

# 6. Monte Carlo Framework

## 6.1 Move (abstract)

Base class for Monte Carlo updates.

Defines:

```cpp
virtual bool attempt(Configuration&, RNG&) = 0;
```

Derived classes will implement:

- kink–antikink insertion  
- deletion  
- worm moves  
- swap moves  
- etc.

---

## 6.2 Estimator (abstract)

Base class for observables.

Defines:

```cpp
virtual double measure(const Configuration&) = 0;
```

Derived classes will implement:

- energy estimator  
- density estimator  
- winding number  
- entanglement estimators  

---

## 6.3 Simulation

High-level driver that:

- holds a `Configuration`  
- applies moves  
- evaluates estimators  

This is the outer loop of the Monte Carlo algorithm.

---

# 7. Testing Infrastructure

A dedicated test suite validates:

### ✔ Initialization correctness  
- exactly one kink per site  
- tau = 0  
- no partners  
- no prev/next  

### ✔ Fock state correctness  
- sum of occupations = N  

### ✔ Structural integrity  
- linked lists valid  
- partner links valid  
- no cycles  
- no invalid indices  

### ✔ Extreme cases  
- all bosons on one site  
- one boson per site  
- zero bosons  

### ✔ Stress test  
- 10,000 random initializations  
- consistency check after each  

All tests passed successfully.

---

# 8. Current Status Summary

You now have:

- a clean, modular architecture  
- a robust worldline data structure  
- correct initialization  
- validated consistency  
- a Monte Carlo framework scaffold  
- a test suite ensuring correctness  

This is the ideal foundation for implementing the first real Monte Carlo move.

---

# 9. Next Steps

You are now ready to implement:

### **Move 1: Kink–Antikink Insertion**

This will test:

- partner links  
- adjacency  
- time sampling  
- acceptance ratios  
- structural updates  

