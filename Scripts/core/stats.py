# core/stats.py
#!/usr/bin/env python3
import numpy as np

def bootstrap_mean_ci(data, n_bootstrap=1000, ci=0.95):
    if len(data) == 0:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng()
    means = []
    for _ in range(n_bootstrap):
        sample = rng.choice(data, size=len(data), replace=True)
        means.append(np.mean(sample))
    mean_val = np.mean(data)
    lower = np.percentile(means, (1 - ci) / 2 * 100)
    upper = np.percentile(means, (1 + ci) / 2 * 100)
    return mean_val, lower, upper