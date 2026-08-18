# core/compute.py
#!/usr/bin/env python3
from config import EQUILIBRATION_NS

def trim_equilibration(t, y, t_eq=None):
    if t_eq is None:
        t_eq = EQUILIBRATION_NS
    mask = t >= t_eq
    return t[mask], y[mask]