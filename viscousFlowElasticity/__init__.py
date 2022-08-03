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
           "tangential_vector_horizontal_plane","vector_product","vector_product_numbers",
           "volumeChangeFromMesh",
           "deformedMeshVolume"
           "loadMesh","get_basic_bilinear_form",
           "get_basic_linear_form","stress_surface_force_linear_form_contribution",
           "normal_vector_outlet",
           "get_targetposition_on_tip_outlet", "get_target_displacement_on_tip_outlet",
           "contribution_to_bilinear_form_shape_tip", "contribution_to_linear_form_shape_tip",
           "run_simulation_step", "getEvacuationTimeConstant", "get_simulation","continue_simulation"
           "getSmoothenedNormalizedAxisymmetricZDeformed",
           "getSmoothenedNormalizedDeformedState",
           "getSmoothenedDeformedState",
           "contribution_to_bilinear_form_shape_tip_fix_coordinates",
           "contribution_to_linear_form_shape_tip_fix_coordinates",
           "contribution_to_linear_form_xy_outlet",
           "contribution_to_bilinear_form_xy_outlet",

           ]



# Import from submodules
from viscousFlowElasticity.linearElasticity import *
from viscousFlowElasticity.smoothing import *
from viscousFlowElasticity.volumeChange import *
from viscousFlowElasticity.plasticDeformation import *
from viscousFlowElasticity.simulation_1mL_minimal import *
from viscousFlowElasticity.geometry import *

