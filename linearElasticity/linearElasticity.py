# The idea here is to optimize the density of the sphere such that given:
# - the touching area (defined as flat)
# - the mechanical properties and terrestrial acceleration
# we can find a solution such that
# - the stress at the edge of the touching area has an average z-component of zero
#
# The last condition is the natural condition for a flat continuation of the surface, that is, the object at equilibrium
# on its touching area


from ngsolve import *
from ngsolve.bvp import BVP
import netgen.meshing as nm
import math
from os import chdir
from netgen.csg import *

def epsilon(u):
    """Strain from deformation field, symmetric"""
    return  0.5*(Grad(u)+Grad(u).trans)

def sigma(strain,E,nu):
    mu = E / 2 / (1 + nu)
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    if strain.dims[0]==strain.dims[1]:
        return lam*Trace(strain)*Id(round(math.sqrt(strain.dim))) + 2*mu*strain
    else:
        if strain.dim==4:
            return lam * (strain[0]+strain[4]) * Id(round(math.sqrt(strain.dim))) + 2 * mu * strain
        else:
            return lam * (strain[0] + strain[4]+strain[8]) * Id(round(math.sqrt(strain.dim))) + 2 * mu * strain

def getElasticBilinearForm(fes,E,nu):
    a = BilinearForm(fes, symmetric=True, condense=True)
    v = fes.TestFunction()
    u = fes.TrialFunction()
    e=epsilon(v)
    s=sigma(epsilon(u),E,nu)
    a+=(s[0]*e[0]+s[1]*e[1]+s[2]*e[2]+s[3]*e[3]+s[4]*e[4]+s[5]*e[5]+s[6]*e[6]+s[7]*e[7]+s[8]*e[8])*dx
    return a



def getVanMises(stress_components):
    return sqrt(0.5*(stress_components[0]-stress_components[4])*(stress_components[0]-stress_components[4])+
                0.5 * (stress_components[4] - stress_components[8]) * (stress_components[4] - stress_components[8])+
                0.5 * (stress_components[0] - stress_components[8]) * (stress_components[0] - stress_components[8])+
                3*(stress_components[1]*stress_components[1]+stress_components[2]*stress_components[2]*stress_components[3]*stress_components[3]))


def extrapolate_E(E, nu, nu_incompressible=0.49):
    return E*(1+nu_incompressible)/(1+nu)
