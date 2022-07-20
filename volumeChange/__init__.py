"""Provide preconfigured functions, linear and bilinear forms for elasticity simualtion
"""

__all__ = ["sigmaVolumetric","getVolumetricChangeContributionForLinearForm"
           ]



# Import from the sub files
from .volumeChange import sigmaVolumetric
from .volumeChange import getVolumetricChangeContributionForLinearForm



def sigmaVolumetric(volumeChange,E,nu,dim=3):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    return lam*volumeChange*Id(dim)

def getVolumetricChangeContributionForLinearForm(fes,E,nu,targetChange):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    v = fes.TestFunction()
    e=epsilon(v)
    s=lam*targetChange*Id(fes.dim) # Stress associated with the volumetric change
    return (s[0]*e[0]+s[4]*e[4]+s[8]*e[8])*dx # only the volumetric part is of interest

