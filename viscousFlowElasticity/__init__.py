"""Collection of functions for simulation of viscoelastic flow in a Bingham plastic

Package written by Thomas Braschler (thomas.braschler@gmail.com)
"""

__all__ = ["epsilon","sigma",
           "getElasticBilinearForm","getVanMises",
           "getMax","getAbsMax",
           "getSmoothened","getSmoothenedNormalized", "getSmoothenedNormalizedAxisymmetric",
           "getSmoothenedStressNormalizedAxisymmetric", "getSmoothenedGradient",
           "getSmoothenedLaplacian"
           ]


# Import from submodules
from linearElasticity import *
from smoothing import *

