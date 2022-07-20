# The idea here is to optimize the density of the sphere such that given:
# - the touching area (defined as flat)
# - the mechanical properties and terrestrial acceleration
# we can find a solution such that
# - the stress at the edge of the touching area has an average z-component of zero
#
# The last condition is the natural condition for a flat continuation of the surface, that is, the object at equilibrium
# on its touching area


from ngsolve import *

from viscousFlowElasticity.linearElasticity import epsilon

def sigmaVolumetric(volumeChange,E,nu,dim=3):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    return lam*volumeChange*Id(dim)

def getVolumetricChangeContributionForLinearForm(fes,E,nu,targetChange):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    v = fes.TestFunction()
    e=epsilon(v)
    s=lam*targetChange*Id(fes.dim) # Stress associated with the volumetric change
    return (s[0]*e[0]+s[4]*e[4]+s[8]*e[8])*dx # only the volumetric part is of interest

