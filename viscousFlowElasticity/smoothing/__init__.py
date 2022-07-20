"""Provide preconfigured functions, linear and bilinear forms for elasticity simualtion
"""

__all__ = ["getMax","getAbsMax",
           "getSmoothened","getSmoothenedNormalized", "getSmoothenedNormalizedAxisymmetric",
           "getSmoothenedStressNormalizedAxisymmetric", "getSmoothenedGradient",
           "getSmoothenedLaplacian"
           ]



# Import from the sub files
from .smoothing import getMax
from .smoothing import getAbsMax
from .smoothing import getSmoothened
from .smoothing import getSmoothenedNormalized
from .smoothing import getSmoothenedNormalizedAxisymmetric
from .smoothing import getSmoothenedStressNormalizedAxisymmetric
from .smoothing import getSmoothenedGradient
from .smoothing import getSmoothenedLaplacian


