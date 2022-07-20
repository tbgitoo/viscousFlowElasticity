"""Provide preconfigured functions, geometrical operations on meshes
"""

__all__ = ["normal_vector_cylinder","normal_vector_torus",
           "get_targetposition_on_torus_surface", "get_targetposition_on_cylinder_surface"
           ]



# Import from the sub files
from .geometry import normal_vector_cylinder
from .geometry import normal_vector_torus
from .geometry import get_targetposition_on_torus_surface
from .geometry import get_targetposition_on_cylinder_surface


