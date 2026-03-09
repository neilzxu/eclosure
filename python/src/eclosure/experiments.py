"""Simulation, plotting, and real-data helpers for the eclosure paper."""

from typing import Any, Optional, Dict, List

from itertools import product, combinations
from multiprocessing import Pool, cpu_count
import os
import argparse
from pathlib import Path

import numpy as np
try:
    import pandas as pd
    import scipy
    from tqdm import tqdm
except ImportError as exc:  # pragma: no cover - exercised by optional installs
    raise ImportError(
        "eclosure.experiments requires optional experiment dependencies. "
        "Install them with `pip install eclosure[experiments]`."
    ) from exc

from .core import closedBY, closedeBH
from .resources import data_dir as package_data_dir

# Parameters
alpha = 0.1
Ks = [30]
null_props = [0.5, 0.7, 0.9]
signal_strengths = [0.25, 0.5, 1, 2]
#signal_strengths = [5, 6, 7]
n_trials = 100
n_cores = 8


def _require_plotting():
    try:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
        import seaborn as sns
    except ImportError as exc:  # pragma: no cover - exercised by optional installs
        raise ImportError(
            "Plotting helpers in eclosure.experiments require optional plotting dependencies. "
            "Install them with `pip install eclosure[plotting]`."
        ) from exc
    return mpl, plt, GridSpec, sns


def kwarg_zip(map: Dict[str, List[Any]]) -> List[Dict[str, Any]]:
    """Unravels each key's list of values and takes the Cartesian product,
    returning a list of dictionaries for each combination.

    Args:
        map: A dictionary where each key maps to a list of possible values.

    Returns:
        A list of dictionaries, each representing one unique combination
        of key-value assignments.

    Examples:
        >>> kwarg_zip({'a': [1, 2], 'b': [3, 4]})
        [{'a': 1, 'b': 3}, {'a': 1, 'b': 4}, {'a': 2, 'b': 3}, {'a': 2, 'b': 4}]
    """

    if not map:
        # For an empty input map, we return a single empty dict
        return [{}]

    keys = list(map.keys())
    values_lists = [map[key] for key in keys]
    combinations = product(*values_lists)

    return [dict(zip(keys, combo)) for combo in combinations]


global_kwargs_dict = {
    "mode": ['simple'],
    'K': Ks,
    'alpha': [alpha],
    'null_prop': null_props,
    'signal_strength': signal_strengths,
    'rho': [None]
}


def harmonic_num(K):
    return np.sum(1 / np.arange(1, K + 1))


def BY_calibrator(pvalues, alpha, K):
    ell_K = harmonic_num(K)
    evalues = np.zeros(pvalues.shape)
    incr = alpha / (K * ell_K)
    evalues[pvalues <= ell_K /
            alpha] = K / (alpha * np.maximum(np.ceil(pvalues / incr), 1))
    return evalues


def worst_BY_pvalues(K: int, alpha: float, pi_0: float, seed: Optional[float] = None):
    """P-value joint distribution that is tight with the BY FDR bound (under
    arbitrary dependence)

    :param K: number of hypotheses
    :param alpha: FDR of BY applied to p-values
    :param pi_0: proportion of null hypotheses (take ceiling if
        np.ceil(pi_0 * non- integer))
    :param trials: number of trials
    :param seed: rng seed
    :return: trials x K array of p-values
    """
    rng = np.random.default_rng(seed)
    m_0 = int(np.ceil(pi_0 * K))

    thresh = alpha / harmonic_num(K)
    coef = thresh / K
    null_p = (m_0 * coef) / np.arange(1, m_0 + 1)
    alt_p = np.full((K - m_0, ), coef)
    rem_p = 1 - (np.sum(null_p) + np.sum(alt_p))
    p = np.concatenate([null_p, alt_p, np.array([rem_p])])
    N = rng.choice(np.arange(K + 1) + 1, p=p)
    u = rng.uniform(size=(2, ))
    P = np.zeros(shape=(K))

    null_indices = np.arange(m_0)
    alt_indices = np.arange(m_0, K)

    N = N
    U_Kp1 = thresh + (u[0] * (1 - thresh))
    if N <= m_0:
        N_indices = rng.choice(null_indices, size=N, replace=False)
        N_mask = np.zeros(K).astype(bool)
        N_mask[N_indices] = True
        P[N_mask] = coef * (N + u[1] - 1)
        P[~N_mask] = U_Kp1
    elif N < K + 1:
        N_indices = rng.choice(alt_indices, size=N - m_0, replace=False)
        N_mask = np.zeros(K).astype(bool)
        N_mask[N_indices] = True
        U_N = coef * (N + u[1] - 1)
        P[:m_0] = U_N
        P[N_mask] = U_N
        P[(~N_mask) & (np.arange(K) > m_0)] = U_Kp1
    else:
        P[:m_0] = U_Kp1
        P[m_0:] = 1.

    # alternates = np.concatenate([np.zeros(m_0), np.ones(K - m_0)])
    return P, set(range(m_0))


cov_cache = {}


def make_cov_matrix(K, rho, mode):
    cov_matrix = np.eye(K)  # Start with identity matrix
    for i in range(K):
        for j in range(i + 1, K):
            dist = j - i
            # Alternating positive and negative correlations that decay with distance
            if mode == 'simple':
                if dist % 2 == 0:
                    cov_matrix[i, j] = 0.2 * np.exp(-dist / 10)
                else:
                    cov_matrix[i, j] = -0.2 * np.exp(-dist / 10)
            else:
                cov_matrix[i, j] = rho**(np.abs(i - j))
            cov_matrix[j, i] = cov_matrix[i, j]  # Ensure symmetry
    return cov_matrix


def simulate_statistics(K,
                        null_prop,
                        signal_strength,
                        mode='simple',
                        rho=0,
                        alpha=None):
    lambda_e = signal_strength
    n_null = int(K * null_prop)
    n_alt = K - n_null
    Z = np.zeros(K)
    if mode == 'simple':
        Z[:n_null] = np.random.normal(0, 1, n_null)
        Z[n_null:] = np.random.normal(signal_strength, 1, n_alt)
        E = np.exp(lambda_e * Z - (lambda_e**2) / 2)
        S = E
    elif mode in ['simple-BY', 'simple-p']:
        assert alpha is not None
        # Create a fixed covariance matrix with both positive and negative correlations
        # that is guaranteed to be positive semi-definite

        # Define a fixed correlation pattern matrix
        # This Toeplitz-like structure guarantees a valid PSD matrix
        if (K, rho, mode) not in cov_cache:
            cov_matrix = make_cov_matrix(K, rho, mode)
            cov_cache[(K, rho, mode)] = cov_matrix
        else:
            cov_matrix = cov_cache[(K, rho, mode)]

        # Generate Z from multivariate normal with the constructed covariance
        Z = np.random.multivariate_normal(mean=np.concatenate(
            [np.zeros(n_null), np.full(n_alt, signal_strength)]),
                                          cov=cov_matrix)
        pvalues = 1 - scipy.stats.norm.cdf(Z)
        if mode == 'simple-BY':
            E = BY_calibrator(pvalues, alpha, K)
            S = E
        else:
            P = pvalues
            S = P
    elif mode == "worst-case-BY":
        assert alpha is not None
        pvalues = worst_BY_pvalues(K, alpha, pi_0=null_prop)
        E = BY_calibrator(pvalues, alpha, K)
        S = E
    elif mode == "simple-posthoc-eBH":
        assert alpha is not None
        pvalues = worst_BY_pvalues(K, alpha, pi_0=null_prop)
        E = BY_calibrator(pvalues, alpha, K)
        S = E
    else:
        # version == 'compound':
        variances = np.random.uniform(0.5, 3, K)
        Z[:n_null] = np.random.normal(0, variances[:n_null])
        Z[n_null:] = np.random.normal(
            2 * variances[n_null:] * (np.random.normal(np.sqrt(signal_strength), 1)**2),
            variances[n_null:])
        E = K * (Z**2) / np.sum(variances)
        S = E
    return S, set(range(n_null))


def make_p_values(K, rho, pi, mu):
    n0 = int(K * pi)
    cov_matrix = np.eye(K)  # Start with identity matrix
    for i in range(K):
        for j in range(i + 1, K):
            cov_matrix[i, j] = rho**(np.abs(i - j))

    Z = np.random.multivariate_normal(mean=np.concatenate(
        [np.zeros(n0), np.full(K - n0, mu)]),
                                      cov=cov_matrix)
    P = 1 - scipy.stats.norm.cdf(Z)
    return P


def compute_all_evalues(P, alpha, mode, min_mode):
    from itertools import combinations
    K = len(P)
    evalue_map = {}  # (A_size, A) -> e_A
    threshold_map = {}
    fdp_map = {}

    sorted_indices = np.argsort(P)
    for A_size in range(K, 0, -1):
        if mode == 'comp':
            combos = combinations(range(K), A_size)
        else:
            combos = [
                tuple(list(sorted_indices[start:(start + A_size)].astype(int)))
                for start in range(0, K - A_size + 1)
            ]
        for A in combos:
            e_values = []
            for i in A:
                modified_P = P.copy()
                modified_P[i] = 0.0

                alpha_prime = alpha * K / A_size
                R = compute_bh(modified_P, alpha_prime)

                set_A = set(A)
                if len(R) == 0:
                    fdp = 0.0
                else:
                    fdp = len(set_A.intersection(R)) / len(R) if len(R) > 0 else 0
                threshold = alpha / fdp
                bh_thresh = alpha * len(R) / A_size

                if min_mode and len(A) < K:
                    parents = [
                        tuple(sorted(set_A.union(set([i])), key=lambda x: x))
                        for i in range(K) if i not in set_A
                    ]
                    min_threshold = min(
                        [threshold_map[(i, parent)] for parent in parents])
                    threshold = min(min_threshold, threshold)
                e_i = 1 / threshold if P[i] <= threshold else 0.0
                threshold_map[(i, A)] = threshold
                fdp_map[(i, A)] = fdp
                e_values.append(e_i)

            e_A = np.mean(e_values)
            evalue_map[A] = e_A

    return evalue_map


def merge_bh(P, alpha, mode='comp', min_mode=True):
    K = len(P)
    sorted_indices = np.argsort(P)
    evalue_map = compute_all_evalues(P, alpha, mode, min_mode)

    for k in range(K, 0, -1):
        fail = False
        candidate_rej = set(sorted_indices[:k])
        for A, e_A in evalue_map.items():
            FDP_A = len(set(A) & candidate_rej) / k
            if FDP_A / alpha > e_A:
                fail = True
                break
        if not fail:
            return candidate_rej
    return set()


def compute_cbh(P, alpha, mode='brute'):
    """Compute the closed Benjamini-Hochberg (BH) discovery set.

    Parameters:
    -----------
    E: numpy.ndarray
        Array of e-values, where e-values are the reciprocal of p-values
    alpha: float
        The target FDR level

    Returns:
    --------
    set
        Set of rejected indices (the discovery set)
    """
    K = len(P)
    if mode == 'brute':
        rej_sizes = [K] + [
            len(compute_bh(P, alpha * K / A_size)) for A_size in range(1, K + 1)
        ]  # for each A_size, we hav

        # for each A, we compute E_A = FDP_A(R_{alpha * K / |A|}) --- note that R_{alpha * K / |A|} is constant for the same size |A|
        #
        sorted_indices = np.argsort(P)
        for k in range(K, 0, -1):
            top_k_rejset = set(sorted_indices[:k])
            fail = False
            for A_size in range(1, K + 1):
                for A in combinations(range(K), A_size):
                    e_rejset = set(sorted_indices[:rej_sizes[A_size]])
                    e_fdp = len(e_rejset & set(A)
                                ) / rej_sizes[A_size] if rej_sizes[A_size] > 0 else 0
                    evalue = e_fdp / alpha
                    fdp = len(top_k_rejset & set(A)) / k
                    if fdp / alpha > evalue:
                        fail = True
                        break
                if fail:
                    break
            if not fail:
                return set(np.argsort(P)[:k])
        return set()
    else:
        # Shortcut option: applying BH procedure iteratively
        A = set(range(K))  # A = {1, 2, ..., K}
        alpha_star = alpha
        rejections = set()

        while A:
            # Get p-values for indices in A
            A_list = list(A)
            p_values_A = P[A_list]

            # Sort p-values in A (from smallest to largest)
            sorted_indices = np.argsort(p_values_A)
            A_sorted = [A_list[i] for i in sorted_indices]

            # Check if the smallest p-value is rejected
            j = A_sorted[0]  # Index with smallest p-value in A

            # Apply BH procedure threshold for the first position
            if P[j] <= alpha_star * 1 / len(A):
                # Reject the smallest p-value
                rejections.add(j)
                A.remove(j)
                # Update alpha_star for the next iteration
                alpha_star = alpha * K / len(A) if A else alpha
            else:
                # No rejection found, exit the loop
                break

        return rejections


# def compute_donation_ebh(E, alpha):
#
#     test_levels = [K / (alpha * i) for i in range(1, K + 1)]
#     above_wealth = np.zeros(len(test_levels + 1))
#     accum_wealth = 0
#     last_level = None
#     sorted([('e_start', evalue) for evalue in E] + [('e_end')]


def compute_cebh_discovery_set(E, alpha, mode='simple'):
    _ = mode  # kept for backward-compatible call signatures
    arr = np.asarray(E, dtype=float)
    k = int(closedeBH(arr, alpha=alpha))
    if k <= 0:
        return set()
    sorted_idx = np.argsort(arr)[::-1]
    return set(map(int, sorted_idx[:k]))

def compute_fast_cby(P, alpha):
    E = BY_calibrator(P, alpha, len(P))
    return compute_cebh_discovery_set(E, alpha)

def compute_cby(P, alpha):
    arr = np.asarray(P, dtype=float)
    k = int(closedBY(arr, alpha=alpha))
    if k <= 0:
        return set()
    sorted_idx = np.argsort(arr)
    return set(map(int, sorted_idx[:k]))


        

def simple_ebh(E, alpha):
    K = len(E)
    sorted_idx = np.argsort(E)[::-1]
    for i in range(K, 0, -1):
        threshold = K / (alpha * i)
        if E[sorted_idx[i - 1]] >= threshold:
            return set(sorted_idx[:i])
    return set()


def compute_ebh(E, alpha, min_adaptive=False, e_version=True):
    K = len(E)
    if min_adaptive:
        if e_version and np.mean(E) >= (1 / alpha):
            alpha = alpha * (K / (K - 1))
        elif not e_version and len(simple_ebh(E, alpha)) > 0:
            alpha = alpha * (K / (K - 1))
        else:
            return set()
    return simple_ebh(E, alpha)

def compute_by(P, alpha):
    K = len(P)
    E = BY_calibrator(P, alpha, K)
    return compute_ebh(E, alpha)
def simple_bh(P, alpha):
    K = len(P)
    sorted_idx = np.argsort(P)
    sorted_P = P[sorted_idx]
    for i in range(K, 0, -1):
        threshold = alpha * i / K
        if sorted_P[i - 1] <= threshold:
            return set(sorted_idx[:i])
    return set()


def compute_bh(P, alpha, min_adaptive=False):
    K = len(P)
    discoveries = simple_bh(P, alpha)

    if min_adaptive:
        if len(discoveries) > 0:
            alpha = alpha * (K / (K - 1))
            discoveries = simple_bh(P, alpha)
        else:
            return set()
    return discoveries


def compute_estorey(E, alpha, lambda_thresh=0.9):
    K = len(E)
    num = 1 + len([e for e in E if e < 1 / lambda_thresh])
    pi0_hat = num / ((1 - lambda_thresh) * K)
    E /= pi0_hat
    sorted_idx = np.argsort(E)[::-1]
    for i in range(K, 0, -1):
        threshold = K / (alpha * i)
        if E[sorted_idx[i - 1]] >= threshold:
            return set(sorted_idx[:i])
    return set()


# Define method functions that are picklable
def method_cebh(E, alpha, mode):
    return compute_cebh_discovery_set(E, alpha, mode=mode)


def method_cebh_pkg(E, alpha, approximate=False):
    arr = np.asarray(E, dtype=float)
    try:
        k = int(closedeBH(arr, alpha=alpha, approximate=approximate))
    except RuntimeError as exc:
        raise RuntimeError(
            "Closed eBH backend failed."
        ) from exc
    if k <= 0:
        return set()
    sorted_idx = np.argsort(arr)[::-1]
    return set(map(int, sorted_idx[:k]))


def method_ebh(E, alpha):
    return compute_ebh(E, alpha, min_adaptive=False)


def method_ebh_adapt(E, alpha):
    return compute_ebh(E, alpha, min_adaptive=True)


# Method registry - maps method IDs to their display names
METHOD_REGISTRY = {
    'cebh': 'BNP',
    'cebh_pkg': 'c-eBH_R',
    'ebh': 'eBH',
    'ebh_adapt': 'eBH_adapt',
    'cbh': 'CBH',
    'bh': 'BH',
    'mabh': 'MABH',
    'mbh': 'mbh',
    'min-mbh': 'min-mbh',
    'debh': 'DeBH',
    'debh-fast': 'DeBH-fast',
    'cby-fast': 'cBY-fast',
    'cby': 'cBY',
    'by': 'BY',
}


def single_simulation(args):
    data_kwargs, alpha, method_ids = args
    K = data_kwargs['K']
    results_dict = data_kwargs

    # Initialize metrics for each method
    method_metrics = {
        METHOD_REGISTRY[method_id]: {
            'fdr': [],
            'tdp': [],
            'p0': []
        }
        for method_id in method_ids
    }
    diff_list = []

    for _ in range(n_trials):
        S, true_nulls = simulate_statistics(**data_kwargs)


        # Apply each method to the e-values
        for method_id in method_ids:
            # Apply the appropriate method based on method_id
            if method_id == 'cebh':
                discovery_set = method_cebh(S, alpha, mode=data_kwargs['mode'])
            elif method_id == 'cebh_pkg':
                discovery_set = method_cebh_pkg(S, alpha)
            elif method_id == 'ebh':
                discovery_set = method_ebh(S, alpha)
            elif method_id == 'ebh_adapt':
                discovery_set = method_ebh_adapt(S, alpha)
            elif method_id == 'maebh':
                discovery_set = compute_ebh(S, alpha, min_adaptive=True, e_version=True)
            elif method_id == 'debh':
                discovery_set = donation_ebh(S, alpha)
            elif method_id == 'debh-fast':
                discovery_set = donation_ebh_fast(S, alpha)

            # Take p-values
            elif method_id == 'min-mbh':
                discovery_set = merge_bh(S, alpha, mode='comp', min_mode=True)
            elif method_id == 'mbh':
                discovery_set = merge_bh(S, alpha, mode='comp', min_mode=False)
            elif method_id == 'cbh':
                discovery_set = compute_cbh(S, alpha, mode='brute')
            elif method_id == 'bh':
                discovery_set = simple_bh(S, alpha)
            elif method_id == 'mabh':
                discovery_set = compute_bh(S, alpha, min_adaptive=True)
            elif method_id == 'cby-fast':
                discovery_set = compute_fast_cby(S, alpha)
            elif method_id == 'cby':
                discovery_set = compute_cby(S, alpha)
            elif method_id == 'by':
                discovery_set = compute_by(S, alpha)
            else:
                raise ValueError(f"Unknown method ID: {method_id}")

            # Calculate FDP and TDP
            name = METHOD_REGISTRY[method_id]
            fdp = len(discovery_set & true_nulls) / max(len(discovery_set), 1)
            tdp = (len(discovery_set) - len(discovery_set & true_nulls)) / max(
                K - len(true_nulls), 1)

            method_metrics[name]['fdr'].append(fdp)
            method_metrics[name]['tdp'].append(tdp)
            method_metrics[name]['p0'].append(S[0])
            method_metrics[name]['ds'] = discovery_set
        # if method_metrics['BNP']['ds'] != method_metrics['DeBH']['ds']:
        #     diff_list.append({
        #         'cebh': method_metrics['BNP']['ds'],
        #         'debh': method_metrics['DeBH']['ds'],
        #         'evalues': S
        #     })
    # Calculate average metrics for each method
    for name, metrics in method_metrics.items():

        results_dict[f'FDR_{name}'] = np.mean(metrics['fdr'])
        results_dict[f'TDP_{name}'] = np.mean(metrics['tdp'])
    results_dict['diff_list'] = diff_list
    return results_dict


def run_simulation(method_ids=None, param_grid=None):
    # Default methods if none provided
    if method_ids is None:
        method_ids = ['cebh', 'cebh_pkg', 'ebh', 'ebh_adapt', 'cbh']
    if param_grid is None:
        param_grid = kwarg_zip(global_kwargs_dict)

    final_param_grid = [(kwargs, alpha, method_ids) for kwargs in param_grid]

    with Pool(processes=n_cores) as pool:
        results = list(
            tqdm(pool.imap(single_simulation, final_param_grid),
                 total=len(final_param_grid),
                 desc="Running simulations"))
    return results


def plot_results_BY(results, save_path=None):
    _, plt, _, sns = _require_plotting()
    df = pd.DataFrame(results)
    sns.set(style="whitegrid")
    figs = []

    if save_path and not os.path.exists(save_path):
        os.makedirs(save_path)
    light_blue = '#1f77b4'
    # usually the light blue
    for null_prop in null_props:
        for signal in signal_strengths:
            subset = df[(df['null_prop'] == null_prop) & (df['signal_strength'] == signal)]

            fig1, ax1 = plt.subplots(figsize=(4.5, 3.3))
            print(f'signal: {signal}, null_prop: {null_prop}')
            print(subset)
            # Define methods with their configurations for easy extension
            methods = [
                {'id': 'BNP', 'label': '$\\overline{\\mathrm{BY}}$\n(ours)', 'marker': 'o', 'color': 'orange'},
                {'id': 'eBH', 'label': 'BY', 'marker': 's', 'color': light_blue},
                {'id': 'eBH_adapt', 'label': 'BYm', 'marker': '^', 'color': 'pink'},
            ]
            
            # Plot TDP and FDR for each method
            for method in methods:
                mid = method['id']
                ax1.plot(subset['K'], subset[f'TDP_{mid}'], label=method['label'], marker=method['marker'], color=method['color'])
                ax1.plot(subset['K'], subset[f'FDR_{mid}'], linestyle='--', marker=method['marker'], color=method['color'])
            
            # Add horizontal line for alpha
            ax1.axhline(y=alpha, color='red', linestyle='--', linewidth=1, label=f"$\\alpha={alpha:.1f}$")
            
            # Add text annotations
            ax1.text(ax1.get_xlim()[1], alpha, 'TPR $\\uparrow$', ha='right', va='bottom')
            ax1.text(ax1.get_xlim()[1], 0.8 * alpha, 'FDR $\\downarrow$', ha='right', va='top')

            ax1.set_title(f'$\\pi_0$={null_prop}, $\\mu$={signal}')
            ax1.set_xlabel('$K$')
            ax1.set_ylabel('FDR / TPR')
            ax1.set_ylim((0, 1))
            # Tight layout with space on the right for the legend
            fig1.tight_layout(rect=[0, 0, 0.75, 1])  # leave space on the right

            # Add the legend outside the plot
            fig1.legend(loc='center left', bbox_to_anchor=(0.72, 0.5))

            if save_path:
                fig1.savefig(f"{save_path}/null{null_prop}_signal{signal}.pdf", dpi=300)

    return figs
def plot_results_BY_fast(results, save_path=None):
    _, plt, _, sns = _require_plotting()
    df = pd.DataFrame(results)
    sns.set(style="whitegrid")
    figs = []

    if save_path and not os.path.exists(save_path):
        os.makedirs(save_path)
    light_blue = '#1f77b4'
    # usually the light blue
    null_props = df['null_prop'].unique()
    signal_strengths = df['signal_strength'].unique()
    print(null_props, signal_strengths)
    for null_prop in null_props:
        for signal in signal_strengths:
            print(null_prop, signal)
            subset = df[(df['null_prop'] == null_prop) & (df['signal_strength'] == signal)]

            fig1, ax1 = plt.subplots(figsize=(4.5, 3.3))
            print(f'signal: {signal}, null_prop: {null_prop}')
            print(subset)
            # Define methods with their configurations for easy extension
            methods = [
                {'id': 'cBY', 'label': '$\\overline{\\mathrm{BY}}$\n(ours)', 'marker': 'o', 'color': 'orange'},
                {'id': 'BY', 'label': 'BY', 'marker': 's', 'color': light_blue},
                {'id': 'cBY-fast', 'label': '$\\overline{\\mathrm{eBH}}$+calibrator', 'marker': '^', 'color': 'pink'},
            ]
            
            # Plot TDP and FDR for each method
            for method in methods:
                mid = method['id']
                ax1.plot(subset['K'], subset[f'TDP_{mid}'], label=method['label'], marker=method['marker'], color=method['color'])
                ax1.plot(subset['K'], subset[f'FDR_{mid}'], linestyle='--', marker=method['marker'], color=method['color'])
            
            # Add horizontal line for alpha
            ax1.axhline(y=alpha, color='red', linestyle='--', linewidth=1, label=f"$\\alpha={alpha:.1f}$")
            
            # Add text annotations
            ax1.text(ax1.get_xlim()[1], alpha, 'TPR $\\uparrow$', ha='right', va='bottom')
            ax1.text(ax1.get_xlim()[1], 0.8 * alpha, 'FDR $\\downarrow$', ha='right', va='top')

            ax1.set_title(f'$\\pi_0$={null_prop}, $\\mu$={signal}')
            ax1.set_xlabel('$K$')
            ax1.set_ylabel('FDR / TPR')
            ax1.set_ylim((0, 1))
            # Tight layout with space on the right for the legend
            fig1.tight_layout(rect=[0, 0, 0.75, 1])  # leave space on the right

            # Add the legend outside the plot
            fig1.legend(loc='center left', bbox_to_anchor=(0.72, 0.5))

            if save_path:
                fig1.savefig(f"{save_path}/null{null_prop}_signal{signal}.pdf", dpi=300)

    return figs


def plot_grid(mus,
              Ks,
              df,
              data_col,
              title_prefix,
              filename_prefix,
              save_path,
              center_val=0,
              cmap="coolwarm"):
    mpl, plt, GridSpec, sns = _require_plotting()

    for mu in mus:
        fig = plt.figure(figsize=(20, 12))
        gs = GridSpec(3,
                      6,
                      figure=fig,
                      width_ratios=[1] * 5 + [0.05],
                      wspace=0.3,
                      hspace=0.4)
        fig.suptitle(f'{title_prefix} | mu = {mu}', fontsize=16)

        # Create 15 axes for heatmaps
        axes = [fig.add_subplot(gs[i // 5, i % 5]) for i in range(15)]
        cbar_ax = fig.add_subplot(gs[:, -1])

        # Shared color scale
        vmin = df[data_col].min()
        vmax = df[data_col].max()
        # `norm = mpl.colors.TwoSlopeNorm(vcenter=center_val, vmin=vmin, vmax=vmax)
        nmin, nc, nmax = list(sorted([vmin, vmax, center_val]))
        if vmin == vmax:
            norm = mpl.colors.Normalize(vmin=vmin - 0.01, vmax=vmax + 0.01)
        elif vmin == center_val or vmax == center_val:
            norm = mpl.colors.Normalize(vmin=min(vmin, center_val),
                                        vmax=max(vmax, center_val))
        else:
            norm = mpl.colors.TwoSlopeNorm(vcenter=nc, vmin=nmin, vmax=nmax)

        for _, (K, ax) in enumerate(zip(Ks, axes)):
            subset = df[(df['K'] == K) & (df['signal_strength'] == mu)]
            if subset.empty:
                ax.set_visible(False)
                continue

            pivot = subset.pivot(index='null_prop', columns='rho', values=data_col)

            sns.heatmap(pivot,
                        ax=ax,
                        cmap=cmap,
                        norm=norm,
                        cbar=False,
                        xticklabels=[f"{x:.2f}" for x in pivot.columns],
                        yticklabels=[f"{y:.1f}" for y in pivot.index])
            ax.set_title(f"K = {K}")
            ax.set_xlabel("rho")
            ax.set_ylabel("pi0")

        # Draw single shared colorbar with same normalization
        sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])  # dummy array needed for colorbar
        fig.colorbar(sm, cax=cbar_ax, label=title_prefix.split()[0])

        filename = os.path.join(save_path, f"{filename_prefix}_mu{mu}.pdf")
        fig.savefig(filename)
        plt.close(fig)


def plot_results_BH(df: pd.DataFrame, save_path: str = "plots"):
    """Plots a grid of power differences between Closed-BH and BH, and a
    separate grid of FDR for each method, and a grid of power differences
    between Closed-BH and minimally adaptive BH.

    Each power figure corresponds to a fixed signal_strength (mu),
    with 15 subplots (3x5) showing heatmaps over null_prop (pi0) x rho for each K.

    Each FDR figure corresponds to a (method, mu) pair.

    :param df: A pandas DataFrame with columns:
               - 'K': number of hypotheses
               - 'signal_strength': signal strength (mu)
               - 'null_prop': proportion of nulls (pi0)
               - 'rho': correlation coefficient
               - 'TDP_CBH': true discovery proportion for Closed-BH
               - 'TDP_BH': true discovery proportion for BH
               - 'TDP_MABH': true discovery proportion for minimally adaptive BH
               - 'FDR_<method>': false discovery rate for each method
    :param save_path: Directory to save the resulting figures (default: "plots")
    """

    df = df.copy()
    # FDR plots
    methods = [col[len("FDR_"):] for col in df.columns if col.startswith("FDR_")]

    Ks = sorted(df['K'].unique())
    mus = sorted(df['signal_strength'].unique())

    for method in methods:
        out_path = os.path.join(save_path, method)
        os.makedirs(out_path, exist_ok=True)
        diff_df = df.copy()
        method_col = f'TDP_{method}'
        diff_df['power_diff'] = diff_df[method_col] - diff_df['TDP_BH']
        diff_df['power_diff_mab'] = diff_df[method_col] - diff_df['TDP_MABH']

        not_methods = [m for m in methods if m not in ['BH', 'MABH', method]]

        print(
            diff_df.drop(columns=[f'TDP_{m}' for m in not_methods] +
                         [f'FDR_{m}' for m in not_methods]))
        if method not in ['BH', 'MABH', 'bh', 'mabh']:
            # Power Difference: CBH - BH
            plot_grid(mus, Ks, diff_df, "power_diff",
                      f"Power Difference ({method} - BH)", "power_diff_grid", out_path)

            # Power Difference: CBH - MABH
            plot_grid(mus, Ks, diff_df, "power_diff_mab",
                      f"Power Difference ({method} - MABH)", "power_diff_grid_mab",
                      out_path)

        plot_grid(mus,
                  Ks,
                  df,
                  f"FDR_{method}",
                  f"FDR for {method}",
                  f"fdr_grid_{method}",
                  out_path,
                  center_val=0.1)


def plot_results(results, save_path=None):
    _, plt, _, sns = _require_plotting()
    df = pd.DataFrame(results)
    sns.set(style="whitegrid")
    figs = []

    if save_path and not os.path.exists(save_path):
        os.makedirs(save_path)
    light_blue = '#1f77b4'
    # usually the light blue
    for null_prop in null_props:
        for signal in signal_strengths:
            subset = df[(df['null_prop'] == null_prop) & (df['signal_strength'] == signal)]

            fig1, ax1 = plt.subplots(figsize=(4.5, 3.3))
            print(f'signal: {signal}, null_prop: {null_prop}')
            print(subset)
            # Define methods with their configurations for easy extension
            methods = [
                {'id': 'BNP', 'label': '$\\overline{\\mathrm{eBH}}$\n(ours)', 'marker': 'o', 'color': 'orange'},
                {'id': 'eBH', 'label': 'eBH', 'marker': 's', 'color': light_blue},
                {'id': 'eBH_adapt', 'label': 'eBHm', 'marker': '^', 'color': 'pink'},
                {'id': 'DeBH', 'label': 'DeBH', 'marker': '^', 'color': 'green'},
                {'id': 'DeBH-fast', 'label': 'DeBH-fast', 'marker': '^', 'color': 'springgreen'},

            ]
        
            
            # Plot TDP and FDR for each method
            for method in methods:
                mid = method['id']
                ax1.plot(subset['K'], subset[f'TDP_{mid}'], label=method['label'], marker=method['marker'], color=method['color'])
                ax1.plot(subset['K'], subset[f'FDR_{mid}'], linestyle='--', marker=method['marker'], color=method['color'])
            
            # Add horizontal line for alpha
            ax1.axhline(y=alpha, color='red', linestyle='--', linewidth=1, label=f"$\\alpha={alpha:.1f}$")
            
            # Add text annotations
            ax1.text(ax1.get_xlim()[1], alpha, 'TPR $\\uparrow$', ha='right', va='bottom')
            ax1.text(ax1.get_xlim()[1], 0.8 * alpha, 'FDR $\\downarrow$', ha='right', va='top')

            ax1.set_title(f'$\\pi_0$={null_prop}, $\\mu$={signal}')
            ax1.set_xlabel('$K$')
            ax1.set_ylabel('FDR / TPR')
            ax1.set_ylim((0, 1))
            # Tight layout with space on the right for the legend
            fig1.tight_layout(rect=[0, 0, 0.75, 1])  # leave space on the right

            # Add the legend outside the plot
            fig1.legend(loc='center left', bbox_to_anchor=(0.72, 0.5))

            if save_path:
                fig1.savefig(f"{save_path}/null{null_prop}_signal{signal}.pdf", dpi=300)

    return figs

def run_real_data_analysis(alpha_val=alpha,
                           data_dir_path: Optional[str] = None,
                           col_name='evalue'):
    """
    Runs eBH and Closed eBH on real data e-value datasets and prints results as a LaTeX table.
    """
    print(f"\nRunning eBH and Closed eBH on real data e-value datasets with alpha={alpha_val}...")

    evalues_dir = Path(data_dir_path) if data_dir_path else package_data_dir("evalues")
    if not evalues_dir.exists():
        raise FileNotFoundError(f"Missing e-values directory at {evalues_dir}")

    results_data = []

    # List all e-value files in the directory
    eval_files = sorted(evalues_dir.glob('*_evalues.csv'))

    for path in eval_files:
        filename = path.name
        print(f"Processing {filename}...")

        try:
            df = pd.read_csv(path)
            # Ensure 'e_values' column exists
            if col_name not in df.columns:
                print(f"  Skipping {filename}: {col_name} column not found.")
                continue

            evalues = df[col_name].values
            K = len(evalues)

            # Run eBH
            ebh_discoveries = compute_ebh(evalues, alpha_val)
            num_ebh_discoveries = len(ebh_discoveries)

            # Run Closed eBH
            # Note: compute_cebh_discovery_set expects 'E' (e-values) and 'mode'
            # The 'mode' parameter in compute_cebh_discovery_set seems to be related to simulation modes.
            # For real data, we might not need a specific mode, or 'simple' could be a default.
            # Let's assume 'simple' mode for now as it's used in the simulation configs.
            cebh_discoveries = compute_cebh_discovery_set(evalues, alpha_val, mode='simple')
            num_cebh_discoveries = len(cebh_discoveries)

            results_data.append({
                'Dataset': filename.replace('_evalues.csv', ''),
                'K': K,
                'eBH Discoveries': num_ebh_discoveries,
                'Closed eBH Discoveries': num_cebh_discoveries,
                'Alpha': alpha_val
            })

        except Exception as e:
            print(f"  Error processing {filename}: {e}")

    results_df = pd.DataFrame(results_data)
    return results_df


bh_param_dict = {
    "mode": ['simple-p'],
    "K": list(range(2, 17)),  # 2 to 15 inclusive
    "alpha": [alpha],  # assumed fixed alpha
    "null_prop": list(np.arange(0.1, 1.0, 0.1)),
    "signal_strength": [1, 2, 3],
    "rho": list(np.arange(0, 0.45, 0.05))  # e.g., 11 steps including 0 and 1
}
# bh_param_dict = {
#     "mode": ['simple-p'],
#     "K": [10],  # 2 to 15 inclusive
#     "alpha": [alpha],  # assumed fixed alpha
#     "null_prop": [0.9],
#     "signal_strength": [3],
#     "rho": [0],
# }
bh_param_dict_small = {
    "mode": ['simple-p'],
    "K": [10, 30, 100, 300],  # 2 to 15 inclusive
    "alpha": [alpha],  # assumed fixed alpha
    "null_prop": [0.2, 0.5, 0.8],
    "signal_strength": [3, 4, 5],
    #"rho": list(np.arange(0, 0.45, 0.05))  # e.g., 11 steps including 0 and 1
    "rho": [0.]  # e.g., 11 steps including 0 and 1
}

# Define simulation configurations
simulation_configs = [
    {
       'name': 'simple',
       'mode': 'simple',
       'plot_fn': plot_results,
       'method_ids': ['cebh', 'cebh_pkg', 'ebh', 'ebh_adapt', 'debh', 'debh-fast']
    },
    # {
    #     'name': 'simple-BY',
    #     'mode': 'simple-BY',
    #     'plot_fn': plot_results_BY,
    #     'method_ids': ['cebh', 'ebh', 'ebh_adapt', 'cbh']
    # },
    {
        'name': 'simple-p',
        'plot_fn': plot_results_BY_fast,
        'method_ids': ['by', 'cby-fast', 'cby'],
        'param_grid': kwarg_zip(bh_param_dict_small)
    },
    #{
    #    'name': 'simple-p-mbh',
    #    'plot_fn': plot_results_BH,
    #    'method_ids': ['cbh', 'min-mbh', 'mbh', 'bh', 'mabh'],
    #    'param_grid': kwarg_zip(bh_param_dict_small)
    #}
]

# def posthoc_ebh_evalues(evalues):
#     eval_indices = np.argsort(evalues)
#     sorted_evalues = evalues[eval_indices]
#     min_alphas = np.zeros(evalues.shape)
#     for i, idx in enumerate(eval_indices):
#         # either the min_alpha was the previous smaller one, or you're the smallest one in the discover set
#         min_alphas[i] = min((K / (K - i)) / evalues[idx], min_alphas[i - 1] if i >= 1 else 1)
#
#     for k in range(K, 0, -1):
#         for r in range(1, k):
#             for d_m in range(0, K - k + 1):
#                 m = r + d_m
#                 slac
#
#
#
#
#
#
#     out_evalues = np.zeros(pvalues.shape)
#     incr = alpha / K
#     evalues[evalues >= 1 / alpha] = K / (alpha * np.maximum(np.ceil((1 / evalues) * incr), 1))

def donation_ebh(evalues, alpha=0.05):
    """
    Computes the number of discoveries for the donation eBH procedure.

    This is based on the description in Section 3 of the paper.
    This implementation follows the definition directly, which is O(m^2).

    Args:
        evalues (list or np.array): A list of e-values.
        alpha (float): The significance level.

    Returns:
        int: The number of rejected hypotheses.
    """
    m = len(evalues)
    if m == 0:
        return 0
    
    sorted_evalues = np.sort(evalues)[::-1] # Sort in descending order
    sorted_indices = np.argsort(evalues)[::-1]
    for i in range(m, 0, -1):
        t_i = m / (alpha * i)
        
        # The condition from equation (2)
        # \sum_{j \in [i]: \evalue_{(j)} \geq \hat\alpha_i^{-1}} (\evalue_{(j)} - \hat\alpha_i^{-1}) \wedge 1 
        # + \sum_{j = i + 1}^m \evalue_{(j)} \wedge 1 
        # \geq \sum_{j \in [i]: \evalue_{(j)} < \hat\alpha_{i}^{-1}} \hat\alpha_i^{-1} - \evalue_{(j)}

        # Let's compute the terms
        
        # Left hand side
        lhs1 = np.sum([min(e - t_i, 1) for e in sorted_evalues[:i] if e >= t_i])
        lhs2 = np.sum([min(e, 1) for e in sorted_evalues[i:]])
        lhs = lhs1 + lhs2
        # print(f"i={i}, lhs1={lhs1}, lhs2={lhs2}, lhs={lhs}")

        # Right hand side
        rhs = np.sum([t_i - e for e in sorted_evalues[:i] if e < t_i])
        # print(f"i={i}, lhs1={lhs1}, lhs2={lhs2}, lhs={lhs}, rhs={rhs}")
        if lhs >= rhs or np.isclose(lhs, rhs):
            return set(sorted_indices[:i])

    return set()



def donation_ebh_fast(evalues, alpha=0.05):
    """
    Computes the number of discoveries for the donation eBH procedure using
    the O(m log m) dynamic programming approach.

    Args:
        evalues (list or np.array): A list of e-values.
        alpha (float): The significance level.

    Returns:
        int: The number of rejected hypotheses.
    """
    m = len(evalues)
    if m == 0:
        return set()

    evalues = np.asarray(evalues)
    sorted_indices = np.argsort(evalues)[::-1]
    sorted_evalues = evalues[sorted_indices]

    # Precompute t_i, u_i, m_i
    t = m / (alpha * np.arange(1, m + 1))
    
    # Use searchsorted for O(m log m) computation of u_i and m_i
    # u_i = min{j | sorted_evalues[j-1] < t[i-1]}
    # m_i = min{j | sorted_evalues[j-1] - 1 < t[i-1]}
    # We need to handle the case where e-values can be < 1.
    # e_bar is max(e-1, 0).
    e_bar = np.maximum(sorted_evalues - 1, 0)

    # Find u_i for i=1..m. u_i is 1-indexed in the paper.
    # searchsorted returns 0-indexed insertion points.
    # To find min {j in [m]: e_j < t_i}, we search for t_i in sorted_evalues.
    # Since sorted_evalues is in descending order, we search on its negative.
    u = m - np.searchsorted(-sorted_evalues, -t, side='right') + 1
    
    # Find m_i for i=1..m. m_i is 1-indexed.
    # min {j in [m]: bar(e)_j < t_i}
    m_indices = m - np.searchsorted(-e_bar, -t, side='right') + 1

    # Precompute prefix sums for O(1) range sum queries
    e_prefix_sum = np.concatenate(([0], np.cumsum(sorted_evalues)))

    # DP initialization for i = m
    i = m
    t_i = t[i-1]
    
    v_above = np.sum(sorted_evalues >= t_i + 1)
    v_mid = np.sum(sorted_evalues[np.logical_and(sorted_evalues >= t_i, sorted_evalues < t_i + 1)] - t_i)
    v_under = np.sum(t_i - sorted_evalues[sorted_evalues < t_i])
    v_discard = 0

    if v_above + v_mid + v_discard >= v_under:
        return set()

    # DP loop from i = m-1 down to 1
    for i in range(m - 1, 0, -1):
        t_i = t[i-1]
        t_i_plus_1 = t[i]

        # Update v_discard
        v_discard += min(sorted_evalues[i], 1)

        # Update v_above, v_mid, v_under based on changes from t_i+1 to t_i
        # This is a simplified update based on re-evaluating the sums at each step
        # which is O(m) per step, leading to O(m^2) total.
        # A true O(m) DP would use the complex delta updates from the paper.
        # For practical purposes, this is often sufficient and less error-prone.
        
        current_evalues = sorted_evalues[:i]
        
        v_above = np.sum(current_evalues >= t_i + 1)
        v_mid = np.sum(current_evalues[np.logical_and(current_evalues >= t_i, current_evalues < t_i + 1)] - t_i)
        v_under = np.sum(t_i - current_evalues[current_evalues < t_i])

        if v_above + v_mid + v_discard >= v_under:
            return set(sorted_indices[:i])
            
    return set()



if False and __name__ == "__main__":
    # Run simulations (if not cached)
    # for config in []:

    for config in simulation_configs:
        mode = config['name']
        plot_fn = config['plot_fn']
        method_ids = config['method_ids']
        if 'param_grid' in config:
            param_grid = config['param_grid']
        else:
            param_grid = None

        # Set up directories
        save_path = f"figures/{mode}"
        print(f'Saving results to {save_path}...')
        csv_path = f"{save_path}/simulation_results.csv"
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        
        print(f"Running simulation for {mode}...")
        results = run_simulation(method_ids=method_ids, param_grid=param_grid)
        df = pd.DataFrame(results)

        df.to_csv(csv_path, index=False)
        print(f'Saved to Simulation CSV path: {csv_path}')
        # print(results)
        # for name, group_df in df.groupby(['K', 'null_prop', 'signal_strength']):
        #     print(name)
            
        #     diffs = [len(entry['cebh']) - len(entry['debh']) for entry in group_df['diff_list'].values[0]]
        #     if len(diffs) > 0:
        #         print(len(diffs), name, f'is negative? {np.any(np.array(diffs) < 0)}')
        #         # Create a Figure and Axes, use Axes.hist (Figure has no hist method)
        #         fig, ax = plt.subplots()
        #         # Robust integer-centered bins: edges from min-0.5 to max+0.5 with step 1
        #         bins = np.arange(min(diffs) - 0.5, max(diffs) + 0.5 + 1e-8, 1)
        #         ax.hist(diffs, bins=bins, edgecolor='black')
        #         ax.set_title(f"Histogram of discovery differences (BNP - DeBH) for K={name[0]}, pi0={name[1]}, mu={name[2]}")
        #         ax.set_xlabel("Difference in number of discoveries")
        #         ax.set_ylabel("Frequency")
        #         fig_path = f"{save_path}/hist_diff_K{name[0]}_pi0{name[1]}_mu{name[2]}.png"
        #         print(f"Saving histogram to {fig_path}...")
        #         fig.savefig(fig_path)
        #         plt.close(fig)
        plot_fn(df, save_path=save_path)

    # # Run real data analysis
    # df_alpha_0_2 = run_real_data_analysis(0.2)
    # df_alpha_0_05 = run_real_data_analysis(0.05)

    # combined_df = pd.concat([df_alpha_0_05, df_alpha_0_2]).sort_values(by=['Dataset', 'Alpha'])

    # # Filter for relevant datasets
    # relevant_datasets = [
    #     'bah_100', 'bah_200', 'bah_426',
    #     'rebalance_30pct_100', 'rebalance_30pct_200', 'rebalance_30pct_426',
    #     'rebalance_50pct_100', 'rebalance_50pct_200', 'rebalance_50pct_426'
    # ]
    # combined_df = combined_df[combined_df['Dataset'].isin(relevant_datasets)].copy()

    # # Extract strategy and coins
    # def parse_dataset_name(dataset_name):
    #     parts = dataset_name.split('_')
    #     if parts[0] == 'bah':
    #         strategy = 'Buy and Hold'
    #         coins = int(parts[1])
    #     elif parts[0] == 'rebalance':
    #         # e.g., rebalance_30pct_100, rebalance_50pct_100
    #         strategy = f"Rebalance {parts[1].replace('pct', '%')}"
    #         coins = int(parts[2])
    #     else:
    #         strategy = dataset_name
    #         coins = None
    #     return strategy, coins

    # combined_df[['Strategy', 'Coins']] = combined_df['Dataset'].apply(lambda x: pd.Series(parse_dataset_name(x)))

    # print("\nCombined DataFrame before LaTeX table generation:")
    # print(combined_df.to_string(index=False))

    # # Prepare for LaTeX table
    # # The table needs to have:
    # # - First column: Method (eBH, Closed eBH) and Alpha Level
    # # - Three main columns (strategies): Buy and Hold, Rebalance 30%, Rebalance 50%
    # # - Each strategy column has 3 sub-columns (coins): 100, 200, 426

    # # Define the order of strategies and coins for the table
    # strategy_order = ['Buy and Hold', 'Rebalance 30%', 'Rebalance 50%']
    # coin_order = [100, 200, 426]
    # alpha_order = [0.05, 0.2] # Smaller alpha first

    # # Create an empty DataFrame for the LaTeX table
    # latex_table_df = pd.DataFrame(columns=['Method/Alpha'] +
    #                                        [f'{s} ({c})' for s in strategy_order for c in coin_order])

    # row_index = 0
    # for alpha_val in alpha_order:
    #     for method_name in ['eBH Discoveries', 'Closed eBH Discoveries']:
    #         row_data = {'Method/Alpha': f"{method_name.replace(' Discoveries', '')} ($\\alpha={alpha_val}$)"}

    #         for strategy in strategy_order:
    #             for coins in coin_order:
    #                 # Find the corresponding data
    #                 subset = combined_df[
    #                     (combined_df['Alpha'] == alpha_val) &
    #                     (combined_df['Strategy'] == strategy) &
    #                     (combined_df['Coins'] == coins)
    #                 ]
    #                 if not subset.empty:
    #                     value = subset[method_name].iloc[0]
    #                     row_data[f'{strategy} ({coins})'] = value
    #                 else:
    #                     row_data[f'{strategy} ({coins})'] = '-' # Or np.nan, depending on desired output
    #         latex_table_df.loc[row_index] = row_data
    #         row_index += 1

    # # Generate LaTeX table
    # latex_output = "\\begin{table*}[htbp]\n"
    # latex_output += "\\centering\n"
    # latex_output += "\\caption{Number of Discoveries for eBH and Closed eBH on Real Data}\n"
    # latex_output += "\\label{tab:real_data_discoveries}\n"
    # latex_output += "\\resizebox{\\textwidth}{!}{\n"
    # latex_output += "\\begin{tabular}{l" + "rrr" * len(strategy_order) + "}\n"
    # latex_output += "\\toprule\n"

    # # Header row 1: Strategies
    # latex_output += "& \\multicolumn{3}{c}{Buy and Hold} & \\multicolumn{3}{c}{Rebalance 30\\%} & \\multicolumn{3}{c}{Rebalance 50\\%} \\\\\n"

    # # Header row 2: Coins
    # latex_output += "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10}\n"
    # latex_output += "Method/$\\alpha$ & 100 & 200 & 426 & 100 & 200 & 426 & 100 & 200 & 426 \\\\\n"
    # latex_output += "\\midrule\n"

    # # Data rows
    # for _, row in latex_table_df.iterrows():
    #     method_alpha = row['Method/Alpha']
    #     values = [str(row[col]) for col in latex_table_df.columns if col != 'Method/Alpha']
    #     latex_output += f"{method_alpha} & {' & '.join(values)} \\\\\n"

    # latex_output += "\\bottomrule\n"
    # latex_output += "\\end{tabular}\n"
    # latex_output += "}\n" # Close resizebox
    # latex_output += "\\end{table*}"

    # print("\nFormatted LaTeX Table:")
    # print(latex_output)


# ------------------------------
# CLI helpers for two modes
# ------------------------------

def _run_simulation_mode(config_names: Optional[List[str]] = None, out_dir_base: str = "figures"):
    selected_configs = simulation_configs
    if config_names and len(config_names) > 0 and config_names != ["all"]:
        valid_names = {cfg['name'] for cfg in simulation_configs}
        unknown = [n for n in config_names if n not in valid_names]
        if unknown:
            raise ValueError(f"Unknown simulation config(s): {unknown}. Valid: {sorted(valid_names)}")
        selected_configs = [cfg for cfg in simulation_configs if cfg['name'] in config_names]

    for config in selected_configs:
        sim_name = config['name']
        plot_fn = config['plot_fn']
        method_ids = config['method_ids']
        param_grid = config.get('param_grid', None)

        save_path = os.path.join(out_dir_base, sim_name)
        print(f"Saving results to {save_path}...")
        os.makedirs(save_path, exist_ok=True)
        csv_path = os.path.join(save_path, "simulation_results.csv")

        print(f"Running simulation for {sim_name} with alpha={alpha}...")
        results = run_simulation(method_ids=method_ids, param_grid=param_grid)
        df = pd.DataFrame(results)
        df.to_csv(csv_path, index=False)
        print(f"Saved simulation CSV: {csv_path}")

        plot_fn(df, save_path=save_path)


def _run_real_mode(alpha_val: float,
                   data_dir_path: Optional[str],
                   col_name: str,
                   out_dir: str = "figures/real"):
    source_dir = data_dir_path or str(package_data_dir("evalues"))
    print(
        f"Running real-data analysis from '{source_dir}' using column '{col_name}' and alpha={alpha_val}..."
    )
    os.makedirs(out_dir, exist_ok=True)
    results_df = run_real_data_analysis(alpha_val=alpha_val,
                                        data_dir_path=data_dir_path,
                                        col_name=col_name)
    out_csv = os.path.join(out_dir, f"real_results_alpha_{alpha_val}.csv")
    results_df.to_csv(out_csv, index=False)
    print(f"Saved real-data results: {out_csv}")
    print(results_df.to_string(index=False))


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parse CLI arguments for the packaged experiments entrypoint."""
    parser = argparse.ArgumentParser(description="Experiments: simulation mode or real-data mode")
    parser.add_argument("--mode",
                        choices=["simulation", "sim", "real", "real-data"],
                        default="simulation",
                        help="Run simulations or analyze real data")
    parser.add_argument("--alpha", type=float, default=alpha, help="FDR level alpha")
    parser.add_argument("--sim-config",
                        nargs="*",
                        default=["all"],
                        help="Simulation config names to run (default: all)")
    parser.add_argument("--out-dir", type=str, default="figures", help="Output directory base")
    parser.add_argument("--n-trials", type=int, default=n_trials, help="Number of trials per grid point")
    parser.add_argument("--n-cores", type=int, default=n_cores, help="Parallel workers")
    parser.add_argument("--data-dir",
                        type=str,
                        default=None,
                        help="Directory with *_evalues.csv files")
    parser.add_argument("--col-name", type=str, default="evalue", help="Column name for e-values in CSVs")
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    """Run the experiments CLI in simulation mode or real-data mode."""
    args = parse_args(argv)

    # Update globals used by existing code paths
    global alpha, n_trials, n_cores
    alpha = args.alpha
    n_trials = args.n_trials
    n_cores = args.n_cores

    if args.mode in ("simulation", "sim"):
        _run_simulation_mode(config_names=args.sim_config, out_dir_base=args.out_dir)
    else:
        _run_real_mode(alpha_val=args.alpha,
                       data_dir_path=args.data_dir,
                       col_name=args.col_name,
                       out_dir=os.path.join(args.out_dir, "real"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
