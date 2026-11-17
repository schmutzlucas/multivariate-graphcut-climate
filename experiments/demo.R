# Level 2: a tiny deterministic example of the actual native graph-cut method.
source("R/load.R")
load_mvgc()
require_production_solver()
set.seed(11)
pmfs <- array(runif(4 * 2 * 8 * 2), c(4, 2, 8, 2))
pmfs <- sweep(pmfs, c(1, 2, 4), apply(pmfs, c(1, 2, 4), sum), "/")
unary <- array(c(.05,.08,.35,.30,.05,.08,.35,.30,
                 .35,.30,.05,.08,.35,.30,.05,.08), c(4, 2, 2))
lat <- c(-30, 30)
solution <- GraphCutHellinger_nD_lat(pmfs, unary, 1, .125, 8L, lat, 1L, FALSE, FALSE)
energy <- graph_cut_energy(solution$label_attribution, unary, pmfs, lat)
out <- mvgc_output_dir("demo")
saveRDS(list(solution = solution, reference_energy = energy, pmfs = pmfs, unary = unary),
        file.path(out, "demo.rds"))
print(solution$label_attribution)
print(energy)
cat("Saved compact native-solver demonstration to", out, "\n")
