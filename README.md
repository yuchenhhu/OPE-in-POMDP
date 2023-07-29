# Off-Policy Evaluation in Partially Observed Markov Decision Processes under Sequential Ignorability

Code for reproducing the numerical experiments in Hu and Wager (2023).

## Reproducing the experiments

Folder Structure:
- `Simple-Example`: code used to reproduce the results in Section 5.1 of Hu and Wager (2023a).
  - `main.jl`
  - `results.jl`
  - `analysis.R`
 
- `mHealth-Study`: code used to reproduce the results in Section 5.2 of Hu and Wager (2023a).
  - `main.jl`
  - `results.jl`
  - `analysis.R`
 
- `Simple-Example-Random`: code used to reproduce the results in Section C of Hu and Wager (2023b).
  - `main.jl`
  - `results.jl`
  - `analysis.R`

Start by running `main.jl` in each folder 100 times with different seeds. The file `results.jl` in each folder collects the results produced by `main.jl`, and the `analysis.R` file produces the plots reported in Hu and Wager (2023a) and Hu and Wager (2023b) based on the pulled results.

## References

HU, Y. and WAGER, S. (2023a). <b>Off-Policy Evaluation in Partially Observed Markov Decision Processes under Sequential Ignorability.</b> [[arXiv](https://arxiv.org/abs/2110.12343)]

HU, Y. and WAGER, S. (2023b). <b>Supplement to ``Off-Policy Evaluation in Partially Observed Markov Decision Processes under Sequential Ignorability.</b> [[arXiv](https://arxiv.org/abs/2110.12343)]
