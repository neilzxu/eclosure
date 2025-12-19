from itertools import product, combinations
import numpy as np
import pandas as pd


def bh_adjust(p: np.ndarray) -> np.ndarray:
    """
    Compute BH-adjusted p-values (Benjamini–Hochberg)."""
    p = np.asarray(p)
    K = p.size
    sorted_idx = np.argsort(p)
    sorted_p = p[sorted_idx]
    adj = sorted_p * K / np.arange(1, K + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    p_adjust = np.empty_like(adj)
    p_adjust[sorted_idx] = adj
    return p_adjust


def bh_reject_count(q: np.ndarray, alpha_level: float) -> int:
    """
    Return BH rejection count on array q at FDR level alpha_level."""
    K = q.size
    sorted_q = np.sort(q)
    r = 0
    for k in range(1, K + 1):
        if sorted_q[k - 1] <= alpha_level * k / K:
            r = k
    return r


def evalue_i(i: int, S: list, p: np.ndarray, alpha: float) -> float:
    """
    Original e-value logic for index i within subset S."""
    # always consider all other p-values except i
    others = np.delete(p, i)
    p2, p3 = np.min(others), np.max(others)
    lS = len(S)
    if lS == 3:
        if p3 <= alpha:
            crit = alpha
        elif p2 <= 2 * alpha / 3:
            crit = 2 * alpha / 3
        else:
            crit = alpha / 3
    elif lS == 2:
        if p2 <= alpha / 3 and p3 <= 3 * alpha / 2:
            crit = 3 * alpha / 2
        elif p2 <= alpha / 3 or p3 <= alpha:
            crit = alpha
        elif p2 <= 2 * alpha / 3:
            crit = 2 * alpha / 3
        else:
            crit = alpha / 3
    elif lS == 1:
        if p3 <= 1.5 * alpha:
            crit = 3 * alpha
        else:
            crit = alpha
    else:
        raise ValueError(f"Unexpected subset size: {lS}")
    return (p[i] <= crit) / crit


def evalue_i_bhcomb(i: int, S: list, p: np.ndarray, alpha: float, m: int = 3) -> float:
    """
    Order-statistic BH-combination e-value: for index i,
    take all other p's, build thresholds j*alpha/|S| for j=1..|S|+1,
    pick the max threshold whose order-statistic condition holds,
    then e_i = I(p_i <= crit)/crit."""
    # consider all p-values except i
    others = np.delete(p, i)
    lS = len(S)
    ps = np.sort(others)
    # thresholds j * alpha / |S| for j=1..|S|+1
    thresholds = [(j * alpha / lS) for j in range(1, lS + 2)]
    valid = []
    for idx, T in enumerate(thresholds, start=1):
        if idx <= lS:
            # require idx-th smallest among others <= T
            if ps[idx - 1] <= T:
                valid.append(T)
        else:
            # require max(other p's) <= T
            if np.max(others) <= T:
                valid.append(T)
    crit = max(valid) if valid else thresholds[0]
    return (p[i] <= crit) / crit


def bh_metrics(p: np.ndarray, alpha: float, m: int = 3):
    """
    Compute BH, step-down BH, MABH, BH2S, and original ePart."""
    p_adjust = bh_adjust(p)
    def r(i):
        return m if i == 0 else int(np.sum(p_adjust <= alpha * m / i))
    bh = r(m)
    i0 = 0
    while r(m - i0) > i0:
        i0 += 1
    sdbh = i0
    mabh = r(m - 1) if bh > 0 else 0
    bh2s = r(m - bh) if bh > 0 else 0
    subsets = [list(c) for r_size in range(1, m + 1)
               for c in combinations(range(m), r_size)]
    ES = np.array([np.mean([evalue_i(i, S, p, alpha) for i in S]) for S in subsets])
    out_js = []
    for j in range(1, m + 1):
        crits = [len([ii for ii in S if ii <= j - 1]) / (j * alpha)
                 for S in subsets]
        if np.all(ES >= np.array(crits) - 1e-10):
            out_js.append(j)
    epart = max([0] + out_js)
    return bh, sdbh, mabh, bh2s, epart


def epart_bhcomb(p: np.ndarray, alpha: float, m: int = 3) -> int:
    """
    Compute ePart using BH-combination e-values for each subset."""
    subsets = [list(c) for r_size in range(1, m + 1)
               for c in combinations(range(m), r_size)]
    ES = np.array([np.mean([evalue_i_bhcomb(i, S, p, alpha, m) for i in S])
                   for S in subsets])
    out_js = []
    for j in range(1, m + 1):
        crits = [len([ii for ii in S if ii <= j - 1]) / (j * alpha)
                 for S in subsets]
        if np.all(ES >= np.array(crits) - 1e-10):
            out_js.append(j)
    return max([0] + out_js)


def generate_p_triplets(alpha: float, m: int = 3, scale: float = 0.99999) -> np.ndarray:
    """
    Generate sorted p-value triplets on a grid."""
    seqp = np.arange(1, 20) / 6 * alpha * scale
    triplets = np.array(list(product(seqp, repeat=m)))
    return np.sort(triplets, axis=1)


def analyze_p_triplets(alpha: float, m: int = 3, scale: float = 0.99999) -> pd.DataFrame:
    """
    Compute BH metrics and compare original vs BH-combination ePart."""
    records = []
    triplets = generate_p_triplets(alpha, m, scale)
    for p in triplets:
        bh, sdbh, mabh, bh2s, epart = bh_metrics(p, alpha, m)
        epart_bc = epart_bhcomb(p, alpha, m)
        rec = {f'p{i+1}': float(p[i]) for i in range(m)}
        rec.update({
            'BH': bh,
            'SD_BH': sdbh,
            'MABH': mabh,
            'BH2S': bh2s,
            'ePart_orig': epart,
            'ePart_bhcomb': epart_bc
        })
        records.append(rec)
        if epart < mabh:
            break
    return pd.DataFrame(records)


def main():
    """
    Run the triplet analysis, generate frequency table of discoveries, and print results."""
    alpha = 0.05
    m = 3
    df = analyze_p_triplets(alpha=alpha, m=m)
    methods = ['BH', 'SD_BH', 'MABH', 'BH2S', 'ePart_orig', 'ePart_bhcomb']
    # Count frequencies for each method and number of discoveries 0..m
    freq = pd.DataFrame(
        {method: df[method].value_counts().reindex(range(m+1), fill_value=0)
         for method in methods},
        index=range(m+1)
    )
    freq.index.name = 'Discoveries'
    print("Frequency of discoveries by method:")
    print(freq)


if __name__ == "__main__":
    main()

