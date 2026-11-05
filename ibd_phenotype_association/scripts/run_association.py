import os
# one thread per process (set before numpy is imported), since the windows are run in parallel processes
for var in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ[var] = "1"
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from common import (COVARIATES, FDR, FOCAL_SPECIES, MIN_PER_SEX, MIN_TEST, RESULTS, TRAITS, WINDOW_COVERAGE, eigen_grm, lmm_lrt, load_bim, load_birds,
                    load_chromosomes, load_loco_grm, load_trimmed_segments, loco_ibd_fractions, make_windows, window_coverage)

N_WORKERS = 30  # species x chromosome jobs run at once

########################### LOAD DATA ###########################
# birds tested (all six measurements and known sex), traits fit within each species
birds = {spp: load_birds(spp) for spp in FOCAL_SPECIES}
chr_names = load_chromosomes()

# fraction of each 100 kb window covered by each tested bird's IBD segments, and each bird's IBD outside every chromosome (LOCO)
windows = make_windows(load_bim())
tested = pd.concat(birds.values())["ind"]
seg = load_trimmed_segments()
coverage = window_coverage(seg[seg["ind"].isin(tested)], windows)
loco_ibd = loco_ibd_fractions()


########################### ASSOCIATION TEST ###########################
def run_chromosome(job):
    """LMM of every trait on IBD in each window of one chromosome, for one species."""
    spp, chrom = job
    start = time.time()
    print(f"[{datetime.now():%H:%M:%S}] START {spp} chr{chr_names[chrom]}", flush=True)

    # birds in the order of this chromosome's LOCO GRM, with their IBD outside the chromosome
    K, grm_ids = load_loco_grm(spp, chrom)
    df = birds[spp].set_index("ind")
    assert set(grm_ids) == set(df.index), f"{spp} birds differ from the GRM birds (rerun 01_make_grm.sh)"
    df = df.loc[grm_ids]
    df["loco_ibd"] = loco_ibd[chrom].reindex(df.index).fillna(0).to_numpy()
    cov = coverage[coverage["ind"].isin(df.index)]
    cov_by_window = dict(list(cov.groupby("window")))

    results = []
    last_toss, eig = None, None
    for w, win in windows[windows["chrom"] == chrom].iterrows():
        # test group: segments cover >= 80% of the window; toss group (left out): some IBD, but < 80%; control group: everyone else
        win_cov = cov_by_window.get(w)
        if win_cov is None:
            continue
        test = set(win_cov.loc[win_cov["coverage"] >= WINDOW_COVERAGE, "ind"])
        toss = frozenset(win_cov.loc[win_cov["coverage"] < WINDOW_COVERAGE, "ind"])

        # skip windows with too few birds in the test group, or too few of either sex
        sex_counts = df.loc[list(test), "sex"].value_counts()
        if len(test) < MIN_TEST[spp] or sex_counts.reindex([0, 1], fill_value=0).min() < MIN_PER_SEX:
            continue

        # relatedness of the birds kept (eigendecomposition reused while the toss group is unchanged, e.g. across neighboring windows)
        keep = ~df.index.isin(toss)
        if toss != last_toss:
            eig = eigen_grm(K[np.ix_(keep, keep)])
            last_toss = toss
        sub = df[keep]
        in_ibd = sub.index.isin(test).astype(float)

        for trait in TRAITS:
            X = np.column_stack([np.ones(len(sub)), in_ibd] + [sub[c].to_numpy(dtype=float) for c in COVARIATES[trait]])
            beta, p = lmm_lrt(sub[trait].to_numpy(dtype=float), X, eig)
            results.append({"species": spp, "trait": trait, "chrom": chr_names[chrom], "window_start": win["start"], "window_end": win["end"],
                            "beta": beta, "p_value": p, "n_test": len(test), "n_control": len(sub) - len(test), "n_toss": len(toss)})

    print(f"[{datetime.now():%H:%M:%S}] DONE  {spp} chr{chr_names[chrom]}: {len(results) // len(TRAITS)} windows tested ({time.time() - start:.0f} s)", flush=True)
    return pd.DataFrame(results)


# species x chromosome jobs, largest chromosomes first so the long jobs start early
chrom_sizes = windows.groupby("chrom").size().sort_values(ascending=False)
jobs = [(spp, chrom) for chrom in chrom_sizes.index for spp in FOCAL_SPECIES]
with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
    results = pd.concat(executor.map(run_chromosome, jobs), ignore_index=True)

# Benjamini-Hochberg FDR within each species x trait (q value: the lowest FDR at which the window would be called significant)
results["q_value"] = results.groupby(["species", "trait"])["p_value"].transform(lambda p: multipletests(p, method="fdr_bh")[1])

# windows in genome order (chromosome order as in the IBIS chromosome map: 1, 1A, 2, 3, ...)
chrom_order = {name: i for i, name in enumerate(chr_names.values())}
results = results.sort_values(["species", "trait", "chrom", "window_start"], key=lambda col: col.map(chrom_order) if col.name == "chrom" else col)
results.to_csv(RESULTS, index=False)
print(f"saved {len(results)} tests to {RESULTS}")
print(results.groupby(["species", "trait"]).agg(windows=("p_value", "size"), significant=("q_value", lambda q: (q < FDR).sum())))
