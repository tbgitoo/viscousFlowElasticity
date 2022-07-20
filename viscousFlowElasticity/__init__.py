"""Collection of functions for simulation of viscoelastic flow in a Bingham plastic

Package written by Thomas Braschler (thomas.braschler@gmail.com)
"""

__all__ = ["epsilon","sigma",
           "getElasticBilinearForm","getVanMises",
           "getMax","getAbsMax",
           "getSmoothened","getSmoothenedNormalized", "getSmoothenedNormalizedAxisymmetric",
           "getSmoothenedStressNormalizedAxisymmetric", "getSmoothenedGradient",
           "getSmoothenedLaplacian",
           "sigmaVolumetric","getVolumetricChangeContributionForLinearForm",
           "getPlasticDeformationContributionForLinearForm","get_yielding_deformation",
           "get_deformation_with_yield_strain",
           "normal_vector_cylinder","normal_vector_torus",
           "get_targetposition_on_torus_surface", "get_targetposition_on_cylinder_surface",
           "loadMesh","get_basic_bilinear_form",
           "get_basic_linear_form","stress_surface_force_linear_form_contribution",
           "normal_vector_outlet",
           "get_targetposition_on_tip_outlet", "get_target_displacement_on_tip_outlet",
           "contribution_to_bilinear_form_shape_tip", "contribution_to_linear_form_shape_tip",
           "run_simulation_step", "getEvacuationTimeConstant", "get_simulation"
           ]



# Import from submodules
from linearElasticity import *
from smoothing import *
from volumeChange import *
from plasticDeformation import *
from simulation_1mL_minimal import *
from geometry import *

