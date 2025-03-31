import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
import os
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
import pandas as pd

# Parameters
alpha = 0.1
Ks = [20, 50, 100, 200]
null_props = [0.5, 0.7, 0.9]
signal_strengths = [2, 3, 4]
n_trials = 1000
n_cores = 12


def simulate_e_values(K, null_prop, signal_strength):
    lambda_e = np.sqrt(np.log(2 * K / alpha))
    n_null = int(K * null_prop)
    n_alt = K - n_null
    Z = np.zeros(K)
    Z[:n_null] = np.random.normal(0, 1, n_null)
    Z[n_null:] = np.random.normal(signal_strength, 1, n_alt)
    E = np.maximum(
        np.exp(lambda_e * Z - (lambda_e ** 2) / 2),
        np.exp(-lambda_e * Z - (lambda_e ** 2) / 2)
    ) / 2
    return E, set(range(n_null))

def compute_bnp_discovery_set(E, alpha):
    K = len(E)
    R_ebh = compute_ebh(E, alpha)
    ebh_k = len(R_ebh)
    start_k = len(R_ebh) + 1

    sorted_idx = np.argsort(E)[::-1]
    E_sorted = E[sorted_idx]
    best_R = set(R_ebh)

    # Precompute full cumulative sum for e-values sorted ascending
    E_sorted_asc = np.sort(E)
    cumsum_full = np.cumsum(E_sorted_asc)

    # Cumulative sum for inside-R worst nulls
    asc_cumsum = np.cumsum(E_sorted_asc)
    cumsum_in = np.cumsum(E_sorted_asc[(-1 * ebh_k):])

    for k in range(start_k, K + 1):
        # new_e = E_sorted[k - 1]  # kth largest is at position K - k in sorted ascending
        # if cumsum_in.size > 0:
        #     new_cumsum = np.append(new_e, cumsum_in + new_e)
        # else:
        #     new_cumsum = np.array([new_e])
        cumsum_in = asc_cumsum[(K - k):(K + 1)] - asc_cumsum[K - k - 1]
        assert len(cumsum_in) == k, (len(cumsum_in), k)

        valid = True
        for r in range(1, k + 1):
            sum_in = cumsum_in[r - 1]
            for m in range(k, K + 1):
                if m - r > 0:
                    sum_out = cumsum_full[m - r - 1]
                else:
                    sum_out = 0
                mean_E = (sum_in + sum_out) / m
                fdp = r / k
                if mean_E < fdp / alpha:
                    valid = False
                    break
            if not valid:
                break
        if valid:
            R_idx = sorted_idx[:k]
            best_R = set(R_idx)

    return best_R

# def compute_bnp_discovery_set(E, alpha):
#     K = len(E)
#     R_ebh = compute_ebh(E, alpha)
#     ebh_k = len(R_ebh)
#     start_k = len(R_ebh) + 1
#
#     sorted_idx = np.argsort(E)[::-1]
#     E_sorted = E[sorted_idx]
#     best_R = set(R_ebh)
#
#     # Precompute full cumulative sum for e-values sorted ascending
#     E_sorted_asc = np.sort(E)
#     cumsum_full = np.cumsum(E_sorted_asc)
#
#     # Cumulative sum for inside-R worst nulls
#     cumsum_in = np.cumsum(E_sorted_asc[(-1 * ebh_k):])
#
#     for k in range(start_k, K + 1):
#         new_e = E_sorted[k - 1]  # kth largest is at position K - k in sorted ascending
#         if cumsum_in.size > 0:
#             new_cumsum = np.append(new_e, cumsum_in + new_e)
#         else:
#             new_cumsum = np.array([new_e])
#         cumsum_in = new_cumsum
#
#         valid = True
#         for r in range(1, k + 1):
#             sum_in = cumsum_in[r - 1]
#             for m in range(k, K + 1):
#                 if m - r > 0:
#                     sum_out = cumsum_full[m - r - 1]
#                 else:
#                     sum_out = 0
#                 mean_E = (sum_in + sum_out) / m
#                 fdp = r / k
#                 if mean_E < fdp / alpha:
#                     valid = False
#                     break
#             if not valid:
#                 break
#         if valid:
#             R_idx = sorted_idx[:k]
#             best_R = set(R_idx)
#
#     return best_R


def compute_ebh(E, alpha):
    K = len(E)
    sorted_idx = np.argsort(E)[::-1]
    for i in range(K, 0, -1):
        threshold = K / (alpha * i)
        if E[sorted_idx[i - 1]] >= threshold:
            return set(sorted_idx[:i])
    return set()


def compute_estorey(E, alpha, lambda_thresh=0.9):
    K = len(E)
    num = 1 + len([e for e in E if e < 1 / lambda_thresh])
    pi0_hat = num / ((1 - lambda_thresh) * K)
    E /= pi0_hat
    sorted_idx = np.argsort(E)[::-1]
    for i in range(K, 0, -1):
        threshold =  K / (alpha * i)
        if E[sorted_idx[i - 1]] >= threshold:
            return set(sorted_idx[:i])
    return set()


def single_simulation(args):
    K, null_prop, signal_strength = args
    fdr_bnp, tdp_bnp = [], []
    fdr_ebh, tdp_ebh = [], []
    fdr_estorey, tdp_estorey = [], []

    for _ in range(n_trials):
        E, true_nulls = simulate_e_values(K, null_prop, signal_strength)
        R_bnp = compute_bnp_discovery_set(E, alpha)
        R_ebh = compute_ebh(E, alpha)
        #R_estorey = compute_estorey(E, alpha)

        fdp_bnp = len(R_bnp & true_nulls) / max(len(R_bnp), 1)
        tdp_bnp.append((len(R_bnp) - len(R_bnp & true_nulls)) / (K - len(true_nulls)))
        fdr_bnp.append(fdp_bnp)

        fdp_ebh = len(R_ebh & true_nulls) / max(len(R_ebh), 1)
        tdp_ebh.append((len(R_ebh) - len(R_ebh & true_nulls)) / (K - len(true_nulls)))
        fdr_ebh.append(fdp_ebh)

        #fdp_estorey = len(R_estorey & true_nulls) / max(len(R_estorey), 1)
        #tdp_estorey.append((len(R_estorey) - len(R_estorey & true_nulls)) / (K - len(true_nulls)))
        #fdr_estorey.append(fdp_estorey)

    return {
        'K': K, 'null_prop': null_prop, 'signal': signal_strength,
        'FDR_BNP': np.mean(fdr_bnp), 'TDP_BNP': np.mean(tdp_bnp),
        'FDR_eBH': np.mean(fdr_ebh), 'TDP_eBH': np.mean(tdp_ebh),
#        'FDR_eStorey': np.mean(fdr_estorey), 'TDP_eStorey': np.mean(tdp_estorey),
    }


def run_simulation():
    param_grid = [(K, null_prop, signal_strength)
                  for K in Ks
                  for null_prop in null_props
                  for signal_strength in signal_strengths]

    with Pool(processes=n_cores) as pool:
        results = list(tqdm(pool.imap(single_simulation, param_grid), total=len(param_grid), desc="Running simulations"))
    return results


def plot_results(results, save_path=None):
    df = pd.DataFrame(results)
    sns.set(style="whitegrid")
    figs = []

    if save_path and not os.path.exists(save_path):
        os.makedirs(save_path)
    light_blue = '#1f77b4'
      # usually the light blue
    for null_prop in null_props:
        for signal in signal_strengths:
            subset = df[(df['null_prop'] == null_prop) & (df['signal'] == signal)]

            fig1, ax1 = plt.subplots(figsize=(4.5, 3.3))

            ax1.plot(subset['K'], subset['TDP_BNP'], label='CeBH', marker='o', color='orange')
            ax1.plot(subset['K'], subset['TDP_eBH'], label='eBH', marker='s', color=light_blue)
            ax1.axhline(y=alpha, color='red', linestyle='--', linewidth=1, label=f"$\\alpha={alpha:.1f}$")
            ax1.plot(subset['K'], subset['FDR_BNP'],linestyle='--', marker='o', color='orange')
            ax1.plot(subset['K'], subset['FDR_eBH'],linestyle='--', marker='s', color=light_blue)
            ax1.text(ax1.get_xlim()[1], alpha, 'TDP $\\uparrow$', ha='right', va='bottom')
            ax1.text(ax1.get_xlim()[1], 0.8 * alpha, 'FDR $\\downarrow$', ha='right', va='top')


            ax1.set_title(f'$\\pi_0$={null_prop}, $\\mu$={signal}')
            ax1.set_xlabel('$K$')
            ax1.set_ylabel('Empirical average')
            ax1.set_ylim((0, 1))
            # Tight layout with space on the right for the legend
            fig1.tight_layout(rect=[0, 0, 0.75, 1])  # leave space on the right

            # Add the legend outside the plot
            fig1.legend(loc='center left', bbox_to_anchor=(0.72, 0.5))

            if save_path:
                fig1.savefig(f"{save_path}/null{null_prop}_signal{signal}.pdf", dpi=300)

    return figs

# Run simulations and plot
csv_path = "figures/gaussian/simulation_results.csv"
if not os.path.exists('figures/gaussian'):
    os.makedirs('figures/gaussian')
if os.path.exists(csv_path):
    print("CSV found, loading cached results...")
    df_cached = pd.read_csv(csv_path)
    plot_results(df_cached.to_dict(orient='records'), save_path="figures/gaussian")
else:
    results = run_simulation()
    df = pd.DataFrame(results)
    df.to_csv(csv_path, index=False)
    plot_results(results, save_path="figures/gaussian")

