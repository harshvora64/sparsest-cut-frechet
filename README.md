# Sparsest Cut Problem

This repository implements exact and approximation algorithms for solving the Sparsest Cut Problem, which involves finding a cut in a graph that minimizes the ratio of cut costs to separated demands.

## Problem Statement
Given:
- A graph \( G = (V, E) \)
- Costs \( c_e \) on edges
- Demands \( d_i \) on vertex pairs \((s_i, t_i)\)
## Features
- **Exact Solution**: Solves the problem using a Mixed Integer Non-Linear Program (MINLP) Solver SCIP.
- **Approximation Algorithm**:
  - LP Relaxation with polynomial constraints.
  - Frechet Embeddings for generating candidate sets.
  - \( O(\log k) \)-approximation with provable bounds.
- **Visualization**: Plots comparing sparsity ratio and runtime for various graph sizes.

## Usage
1. Install the required packages:
```bash
pip install -r requirements.txt
```

2. Run graphgen.py to generate a random graph:
```bash
python graphgen.py
```

3. Run main.py to solve the Sparsest Cut Problem using the approximation algorithm:
```bash
python frechet_code.py
```

4. Run exact.py to solve the Sparsest Cut Problem using the exact algorithm:
```bash
python scip_demands.py
```
