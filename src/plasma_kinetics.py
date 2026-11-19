"""Magnetron sputtering plasma kinetics calculations (mean free path and collision probability)."""

import numpy as np

K_B = 1.380649e-23  # Boltzmann constant (J/K)
TORR_TO_PA = 133.322  # Conversion factor from Torr to Pascals


def calculate_mean_free_path(
    temperature_k: float | np.ndarray,
    pressure_pa: float | np.ndarray,
    r_mean: float = 1.5e-10,
) -> float | np.ndarray:
    """Compute the mean free path (lambda) of sputtered species in an Ar/N2 plasma.

    Formula:
        lambda = (k_B * T) / (sqrt(2) * 4 * pi * r_mean^2 * P)

    Parameters
    ----------
    temperature_k : float or np.ndarray
        Gas/chamber temperature in Kelvin (K).
    pressure_pa : float or np.ndarray
        Working chamber pressure in Pascals (Pa).
    r_mean : float, optional
        Effective mean atomic/molecular radius in meters (default: 1.5e-10 m).

    Returns
    -------
    float or np.ndarray
        Mean free path in meters (m).
    """
    cross_section = np.sqrt(2.0) * 4.0 * np.pi * (r_mean ** 2)
    return (K_B * temperature_k) / (cross_section * pressure_pa)


def calculate_collision_probability(
    distance_m: float | np.ndarray,
    mean_free_path_m: float | np.ndarray,
) -> float | np.ndarray:
    """Compute the gas-phase collision probability (P_c) across target-to-substrate distance.

    Formula:
        P_c = 1 - exp(-D_TS / lambda)

    Parameters
    ----------
    distance_m : float or np.ndarray
        Target-to-substrate distance D_TS in meters (m).
    mean_free_path_m : float or np.ndarray
        Mean free path lambda in meters (m).

    Returns
    -------
    float or np.ndarray
        Dimensionless collision probability in [0, 1].
    """
    return 1.0 - np.exp(-distance_m / mean_free_path_m)
