"""Provide preconfigured functions, linear and bilinear forms for elasticity simualtion
"""

__all__ = ["getPlasticDeformationContributionForLinearForm","get_yielding_deformation",
           "get_deformation_with_yield_strain"
           ]



# Import from the sub files
from .plasticDeformation import getPlasticDeformationContributionForLinearForm
from .plasticDeformation import get_yielding_deformation
from .plasticDeformation import get_deformation_with_yield_strain



