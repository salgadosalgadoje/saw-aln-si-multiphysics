"""Rotation matrices, Bond transformations, and Rayleigh wave acoustic confinement for AlN/Si(111)."""

import numpy as np


def euler_rotation_matrix(phi: float, theta: float, psi: float) -> np.ndarray:
    """Compute 3x3 Z-X-Z Euler rotation matrix (angles in radians)."""
    c1, s1 = np.cos(phi), np.sin(phi)
    c2, s2 = np.cos(theta), np.sin(theta)
    c3, s3 = np.cos(psi), np.sin(psi)

    r_z1 = np.array([[c1, s1, 0.0], [-s1, c1, 0.0], [0.0, 0.0, 1.0]])
    r_x = np.array([[1.0, 0.0, 0.0], [0.0, c2, s2], [0.0, -s2, c2]])
    r_z2 = np.array([[c3, s3, 0.0], [-s3, c3, 0.0], [0.0, 0.0, 1.0]])

    return r_z2 @ r_x @ r_z1


def bond_matrix(R: np.ndarray) -> np.ndarray:
    """Construct the 6x6 Bond transformation matrix M from a 3x3 rotation matrix R."""
    M = np.zeros((6, 6), dtype=float)

    # Top-left 3x3: squared elements
    M[:3, :3] = R**2

    # Top-right 3x3: 2 * R_ik * R_il
    M[:3, 3] = 2.0 * R[:, 1] * R[:, 2]
    M[:3, 4] = 2.0 * R[:, 0] * R[:, 2]
    M[:3, 5] = 2.0 * R[:, 0] * R[:, 1]

    # Bottom-left 3x3: R_ik * R_jk
    M[3, :3] = R[1, :] * R[2, :]
    M[4, :3] = R[0, :] * R[2, :]
    M[5, :3] = R[0, :] * R[1, :]

    # Bottom-right 3x3
    idx_pairs = [(1, 2), (0, 2), (0, 1)]
    for row_idx, (i, j) in enumerate(idx_pairs, start=3):
        for col_idx, (k, l) in enumerate(idx_pairs, start=3):
            M[row_idx, col_idx] = R[i, k] * R[j, l] + R[i, l] * R[j, k]

    return M


def rotate_stiffness_tensor(C: np.ndarray, R: np.ndarray) -> np.ndarray:
    """Rotate a 6x6 Voigt stiffness tensor C using Bond transformation: C' = M @ C @ M.T."""
    M = bond_matrix(R)
    return M @ C @ M.T


def get_material_tensors() -> dict:
    """Return elastic stiffness tensors (in Pa) and densities (kg/m^3) for AlN and Si."""
    # Hexagonal wurtzite AlN (c-axis along x3)
    c11_aln, c12_aln, c13_aln, c33_aln, c44_aln = 410e9, 149e9, 99e9, 389e9, 125e9
    c66_aln = 0.5 * (c11_aln - c12_aln)
    C_aln = np.array([
        [c11_aln, c12_aln, c13_aln, 0, 0, 0],
        [c12_aln, c11_aln, c13_aln, 0, 0, 0],
        [c13_aln, c13_aln, c33_aln, 0, 0, 0],
        [0, 0, 0, c44_aln, 0, 0],
        [0, 0, 0, 0, c44_aln, 0],
        [0, 0, 0, 0, 0, c66_aln],
    ], dtype=float)

    # Cubic Si in standard [100] crystallographic axes
    c11_si, c12_si, c44_si = 165.7e9, 63.9e9, 79.6e9
    C_si = np.array([
        [c11_si, c12_si, c12_si, 0, 0, 0],
        [c12_si, c11_si, c12_si, 0, 0, 0],
        [c12_si, c12_si, c11_si, 0, 0, 0],
        [0, 0, 0, c44_si, 0, 0],
        [0, 0, 0, 0, c44_si, 0],
        [0, 0, 0, 0, 0, c44_si],
    ], dtype=float)

    return {
        "AlN": {"C": C_aln, "rho": 3260.0},
        "Si": {"C": C_si, "rho": 2329.0},
    }


def si_111_rotation_matrix(in_plane_angle: float = 0.0) -> np.ndarray:
    """Transform cubic [100] axes to Si(111) surface normal (x3 || [111])."""
    e3 = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)       # [111]
    e1_0 = np.array([1.0, 1.0, -2.0]) / np.sqrt(6.0)    # [11-2]
    e2_0 = np.cross(e3, e1_0)                           # [-110]

    R_base = np.vstack([e1_0, e2_0, e3])
    c, s = np.cos(in_plane_angle), np.sin(in_plane_angle)
    R_az = np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])
    return R_az @ R_base


def rayleigh_depth_profile(
    z_norm: np.ndarray,
    v_r: float,
    v_l: float,
    v_t: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute normalized Rayleigh displacements (u1, u3) and acoustic energy density vs z/lambda."""
    q = np.sqrt(max(1.0 - (v_r / v_l) ** 2, 1e-9))
    s = np.sqrt(max(1.0 - (v_r / v_t) ** 2, 1e-9))
    k_z = 2.0 * np.pi * z_norm

    ratio = (2.0 * q * s) / (1.0 + s**2)
    u1 = np.exp(-q * k_z) - ratio * np.exp(-s * k_z)
    u3 = -q * (np.exp(-q * k_z) - (2.0 / (1.0 + s**2)) * np.exp(-s * k_z))

    # Normalize against surface vertical amplitude |u3(0)|
    norm = np.abs(u3[0]) if np.abs(u3[0]) > 0 else 1.0
    u1_norm = u1 / norm
    u3_norm = u3 / norm

    energy_density = u1_norm**2 + u3_norm**2
    energy_density /= energy_density[0]
    return u1_norm, u3_norm, energy_density


def acoustic_confinement_factor(
    h_over_lambda: np.ndarray,
    v_r: float,
    v_l: float,
    v_t: float,
) -> np.ndarray:
    """Fraction of Rayleigh acoustic energy confined within normalized film thickness h/lambda."""
    z_grid = np.linspace(0.0, 5.0, 5000)
    _, _, energy = rayleigh_depth_profile(z_grid, v_r, v_l, v_t)
    total_energy = np.trapezoid(energy, z_grid)

    gamma = np.zeros_like(h_over_lambda, dtype=float)
    for idx, h_norm in enumerate(h_over_lambda):
        mask = z_grid <= h_norm
        if np.any(mask):
            gamma[idx] = np.trapezoid(energy[mask], z_grid[mask]) / total_energy
    return gamma