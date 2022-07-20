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
           "get_deformation_with_yield_strain"
           ]


# Import from submodules
from linearElasticity import *
from smoothing import *
from volumeChange import *
from plasticDeformation import *
from simulation_1mL_minimal import *

