# The idea here is to optimize the density of the sphere such that given:
# - the touching area (defined as flat)
# - the mechanical properties and terrestrial acceleration
# we can find a solution such that
# - the stress at the edge of the touching area has an average z-component of zero
#
# The last condition is the natural condition for a flat continuation of the surface, that is, the object at equilibrium
# on its touching area


from ngsolve import *

from viscousFlowElasticity.linearElasticity import getVanMises
from viscousFlowElasticity.linearElasticity import epsilon
from viscousFlowElasticity.linearElasticity import sigma


def getPlasticDeformationContributionForLinearForm(fes,E,nu,targetChange):
    v = fes.TestFunction()
    e=epsilon(v)
    s=sigma(targetChange,E,nu) # Stress associated with the volumetric change
    return (s[0]*e[0]+s[1]*e[1]+s[2]*e[2]+s[3]*e[3]+s[4]*e[4]+s[5]*e[5]+s[6]*e[6]+s[7]*e[7]+s[8]*e[8])*dx



def get_yielding_deformation(mesh,strain,stress,E,vanMisesRel=1.5):
    fesp = H1(mesh, order=1, dim=9)
    output=GridFunction(fesp)
    ind=0
    for v in mesh.vertices:
        stress_v= stress(mesh(v.point[0],v.point[1],v.point[2]))
        strain_v= strain(mesh(v.point[0],v.point[1],v.point[2]))
        vM=getVanMises(stress_v)
        volume_strain=(strain_v[0]+strain_v[4]+strain_v[8])/3
        if vM > E*vanMisesRel:
            for ind2 in range(9):
                output.vec.data[ind][ind2]=strain_v[ind2]
            output.vec.data[ind][0]=output.vec.data[ind][0]-volume_strain # only deviatory, not compressive strain
            output.vec.data[ind][4] = output.vec.data[ind][4] - volume_strain  # only deviatory, not compressive strain
            output.vec.data[ind][8] = output.vec.data[ind][8] - volume_strain  # only deviatory, not compressive strain
        ind=ind+1
    return output


def get_deformation_with_yield_strain(intended_plastic_deformation,yield_strain):
    return IfPos(getVanMises(intended_plastic_deformation)-yield_strain,(getVanMises(intended_plastic_deformation)-yield_strain) / (
                getVanMises(intended_plastic_deformation)) * intended_plastic_deformation,CoefficientFunction((0, 0, 0, 0, 0, 0, 0, 0, 0)))





