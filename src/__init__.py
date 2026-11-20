"""Core physics modules for AlN/Si(111) SAW modeling and sputtering kinetics."""

from .plasma_kinetics import (
    K_B,
    TORR_TO_PA,
    calculate_mean_free_path,
    calculate_collision_probability,
)

__all__ = [
    "K_B",
    "TORR_TO_PA",
    "calculate_mean_free_path",
    "calculate_collision_probability",
]

from .plasma_kinetics import (
    K_B,
    TORR_TO_PA,
    calculate_mean_free_path,
    calculate_collision_probability,
)
from .rayleigh_confinement import (
    euler_rotation_matrix,
    bond_matrix,
    rotate_stiffness_tensor,
    get_material_tensors,
    si_111_rotation_matrix,
    rayleigh_depth_profile,
    acoustic_confinement_factor,
)