# TODO List

- Affinity model - use ODE solver

- More Constant
    - V = baseline
    - B = 0
    - M = previous step

- Differential 
    - dV = 2V - model responses (inverse distance times B-cells)
    - dB = V/D + XMV/D - decay (where X is memory-to-Bcell transition)
    - dM = aB - XMV/D (where a is Bcell-to-memory transition)


ultimately, at the end, we show the memory-cell levels of individuals under different treatment histories
- Total M-cell response? Sum equal?