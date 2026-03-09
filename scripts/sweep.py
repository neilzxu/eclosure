"""Exploratory parameter sweeps and plotting helpers for simulation studies."""

import argparse
from itertools import combinations
from multiprocessing import Pool, cpu_count
import os
from typing import Optional

import numpy as np

try:
    import matplotlib.pyplot as plt
    import pandas as pd
    import scipy
    import seaborn as sns
    from tqdm import tqdm
except ImportError as exc:  # pragma: no cover - exercised by optional installs
    raise ImportError(
        "scripts/sweep.py requires plotting and experiment dependencies. "
        "Install them with `pip install -e .[experiments,plotting]`."
    ) from exc


def make_p_values(K, rho, pi, mu):
    n0 = int(K * pi)
    cov_matrix = np.eye(K)
    for i in range(K):
        for j in range(i + 1, K):
            cov_matrix[i, j] = rho ** np.abs(i - j)

    Z = np.random.multivariate_normal(
        mean=np.concatenate([np.zeros(n0), np.full(K - n0, mu)]),
        cov=cov_matrix,
    )
    P = 1 - scipy.stats.norm.cdf(Z)
    return P


def simple_bh(P, alpha):
    sorted_idx = np.argsort(P)
    sorted_P = P[sorted_idx]
    for i in range(len(P), 0, -1):
        threshold = alpha * i / len(P)
        if sorted_P[i - 1] <= threshold:
            return set(sorted_idx[:i])
    return set()


def closed_bh_brute(P, alpha):
    K = len(P)
    rej_sizes = [K] + [len(simple_bh(P, alpha * K / A_size)) for A_size in range(1, K + 1)]
    sorted_indices = np.argsort(P)

    for k in range(K, -1, -1):
        fail = False
        top_k_rejset = set(sorted_indices[:k])
        for A_size in range(1, K + 1):
            for A in combinations(range(K), A_size):
                e_rejset = set(sorted_indices[:rej_sizes[A_size]])
                e_fdp = len(e_rejset & set(A)) / rej_sizes[A_size] if rej_sizes[A_size] > 0 else 0
                evalue = e_fdp / alpha
                fdp = len(top_k_rejset & set(A)) / k if k > 0 else 0
                if fdp / alpha > evalue:
                    fail = True
                break
            if fail:
                break
        if not fail:
            return set(np.argsort(P)[:k])
    return set()


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse CLI arguments for the exploratory sweep script."""
    parser = argparse.ArgumentParser(description="Run a small closed-BH brute-force sweep example.")
    parser.add_argument("--K", type=int, default=6, help="Number of hypotheses")
    parser.add_argument("--rho", type=float, default=0.2, help="Correlation decay parameter")
    parser.add_argument("--pi", type=float, default=0.5, help="Null proportion")
    parser.add_argument("--mu", type=float, default=2.0, help="Alternative mean shift")
    parser.add_argument("--alpha", type=float, default=0.05, help="Target FDR level")
    parser.add_argument("--seed", type=int, default=0, help="Random seed")
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    """Run a small exploratory brute-force closed-BH example and print results."""
    args = parse_args(argv)
    np.random.seed(args.seed)
    pvalues = make_p_values(args.K, args.rho, args.pi, args.mu)
    bh = simple_bh(pvalues, args.alpha)
    closed = closed_bh_brute(pvalues, args.alpha)
    print("p-values:", np.array2string(np.sort(pvalues), precision=4))
    print("BH discoveries:", sorted(map(int, bh)))
    print("closed-BH discoveries:", sorted(map(int, closed)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
