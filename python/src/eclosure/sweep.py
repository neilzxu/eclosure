from typing import Optional
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
import os
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
import pandas as pd
import scipy


def make_p_values(K, rho, pi, mu):
    n0 = int(K * pi)
    cov_matrix = np.eye(K)  # Start with identity matrix
    for i in range(K):
        for j in range(i + 1, K):
            cov_matrix[i, j] = rho ** (np.abs(i - j))

    Z = np.random.multivariate_normal(mean=np.concatenate(
        [np.zeros(n0),
         np.full(K - n0, mu)]), cov=cov_matrix)
    P = 1 - scipy.stats.norm.cdf(Z)
    return P

def closed_bh_brute(P, alpha):
    K = len(P)
    rej_sizes = [K] + [
        len(compute_ebh(E, alpha * K / A_size))
        for A_size in range(1, K + 1)
    ] # for each A_size, we hav

    # for each A, we compute E_A = FDP_A(R_{alpha * K / |A|}) --- note that R_{alpha * K / |A|} is constant for the same size |A|
    #
    for k in range(K, -1, -1):
        fail = False
        for A_size in range(1, K + 1):
            for start in range(0, K - A_size):
                e_fdp = min(rej_size[A_size] - start, A_size) / rej_size[A_size]
                evalue = e_fdp / alpha
                fdp = min(k - start, A_size) / k
                if fdp / alpha > e_fdp:
                    fail = True
                break
            if fail:
                break
        if not fail:
            return set(np.argsort(P)[:k])
    return set()


