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

from viscousFlowElasticity.linearElasticity import getElasticBilinearForm
from viscousFlowElasticity.linearElasticity import epsilon
from viscousFlowElasticity.linearElasticity import sigma
from viscousFlowElasticity.linearElasticity import extrapolate_E
from viscousFlowElasticity.linearElasticity import getVanMises

from viscousFlowElasticity.smoothing import getSmoothenedNormalizedAxisymmetricZDeformed
from viscousFlowElasticity.smoothing import getSmoothenedNormalizedAxisymmetric
from viscousFlowElasticity.smoothing import getAbsMax

from viscousFlowElasticity.plasticDeformation import get_deformation_with_yield_strain
from viscousFlowElasticity.plasticDeformation import getPlasticDeformationContributionForLinearForm

from viscousFlowElasticity.volumeChange import getVolumetricChangeContributionForLinearForm
from viscousFlowElasticity.volumeChange import sigmaVolumetric

from viscousFlowElasticity.geometry import normal_vector_cylinder
from viscousFlowElasticity.geometry import normal_vector_torus
from viscousFlowElasticity.geometry import get_targetposition_on_torus_surface
from viscousFlowElasticity.geometry import get_targetposition_on_cylinder_surface


# Filename is a .vol file as generated by Netgen, save mesh command
# this is specific to the syringe geometry, where there is a plunger that pushes down,
# an outer barrel surface (cylinder)
# a tip out surface containing the transition of barrel to tip
# and the outlet
def loadMesh(filename):
    m = nm.Mesh(3)
    m.Load(filename)
    m.SetBCName(0, "plunger")
    m.SetBCName(1, "barrel_outer")
    m.SetBCName(2, "tip_outer")
    m.SetBCName(3, "tip_outlet")
    mesh = Mesh(m)
    return mesh



def get_basic_bilinear_form(fes,mesh,E,nu=0.49,plunger_z_penalty_factor=1e10,xy_plunger_penalty_factor=1e10,r_penalty_factor_barrel=1e10):
    a_incompressible=getElasticBilinearForm(fes,E,nu)
    penalty_z_dict = {"plunger": plunger_z_penalty_factor,"barrel_outer": 0,"tip_outer": 0, "tip_outlet": 0
                      }
    penalty_z = CoefficientFunction([penalty_z_dict[bound] for bound in mesh.GetBoundaries()])
    penalty_xy_dict = {"plunger": xy_plunger_penalty_factor,"barrel_outer": r_penalty_factor_barrel,"tip_outer": 0, "tip_outlet": 0
                      }
    penalty_xy = CoefficientFunction([penalty_xy_dict[bound] for bound in mesh.GetBoundaries()])
    v = fes.TestFunction()
    u = fes.TrialFunction()
    a_incompressible += penalty_z * u[2] * v[2] * ds
    a_incompressible += penalty_xy * u[0] * v[0] * ds
    a_incompressible += penalty_xy * u[1] * v[1] * ds
    return a_incompressible



def get_basic_linear_form(fes,mesh,plunger_z_penalty_factor=1e10,plunger_displacement=0):
    penalty_z_dict_linear = {"plunger": plunger_z_penalty_factor*plunger_displacement,"barrel_outer": 0,"tip_outer": 0, "tip_outlet": 0
                      }
    penalty_z_linear = CoefficientFunction([penalty_z_dict_linear[bound] for bound in mesh.GetBoundaries()])
    f = LinearForm(fes)
    v = fes.TestFunction()
    f += penalty_z_linear*v[2]*ds
    return(f)



def stress_surface_force_linear_form_contribution(stress_deformatory,mesh,fes):
    v = fes.TestFunction()
    slip_dict = {"tip_outer": 1, "tip_outlet": 0, "barrel_outer": 1, "plunger": 0,
                      "barrel_front": 1, "transition": 1}
    slip = CoefficientFunction([slip_dict[bound] for bound in mesh.GetBoundaries()])
    return slip*stress_surface_force_tangential(stress_deformatory,mesh)*v*ds




def normal_vector_outlet(u,torus_position_z=-0.001,radius_ring=0.0015):
    r_xy = sqrt((u[0] + x) * (u[0] + x) + (u[1] + y)* (u[1] + y))
    return IfPos(r_xy - radius_ring, CoefficientFunction((0,0,1)),
                 IfPos(torus_position_z - u[2] - z, -normal_vector_cylinder(u),
                       normal_vector_torus(u, torus_position_z=torus_position_z, radius_ring=radius_ring)))







def get_targetposition_on_tip_outlet(u,torus_position_z=-0.001,radius_ring=0.0015,radius_torus=0.001,radius_outlet=0.0005):
    r_xy=sqrt((u[0]+x)*(u[0]+x)+(u[1]+y)*(u[1]+y))
    return IfPos(r_xy-radius_ring,CoefficientFunction((u[0]+x,u[1]+y,torus_position_z+radius_torus)),IfPos(torus_position_z-u[2]-z,get_targetposition_on_cylinder_surface(u,radius=radius_outlet),
          get_targetposition_on_torus_surface(u,torus_position_z=torus_position_z,radius_ring=radius_ring,
                                          radius_torus=radius_torus)))

def get_target_displacement_on_tip_outlet(u,torus_position_z=-0.001,radius_ring=0.0015,radius_torus=0.001,radius_outlet=0.0005):
    return get_targetposition_on_tip_outlet(u,
        torus_position_z=torus_position_z,radius_ring=radius_ring,radius_torus=radius_torus,radius_outlet=radius_outlet)-\
        CoefficientFunction((x,y,z))



def contribution_to_bilinear_form_shape_tip(fes,mesh,u_last_step,torus_position_z=-0.001,radius_ring=0.0015,r_penalty_factor=1e10):
    penalty_r_dict = {"tip_outer": r_penalty_factor, "tip_outlet": 0, "barrel_outer": 0,
                      "plunger": 0
                      }
    penalty_r = CoefficientFunction([penalty_r_dict[bound] for bound in mesh.GetBoundaries()])
    u = fes.TrialFunction()
    v = fes.TestFunction()
    return penalty_r * (
            normal_vector_outlet(u_last_step, torus_position_z=torus_position_z,
                                        radius_ring=radius_ring) * u) * (
            normal_vector_outlet(u_last_step, torus_position_z=torus_position_z,
                                        radius_ring=radius_ring) * v) * ds


def contribution_to_bilinear_form_shape_tip_fix_coordinates(fes,mesh,r_penalty_factor=1e10):
    penalty_r_dict = {"tip_outer": r_penalty_factor, "tip_outlet": 0, "barrel_outer": 0,
                      "plunger": 0
                      }
    penalty_r = CoefficientFunction([penalty_r_dict[bound] for bound in mesh.GetBoundaries()])
    u = fes.TrialFunction()
    v = fes.TestFunction()
    return penalty_r * (u[0]*v[0]+u[1]*v[1]+u[2]*v[2]) * ds






def contribution_to_linear_form_shape_tip(fes,mesh,u_last_step,torus_position_z=-0.001,radius_ring=0.0015,
                                          radius_torus=0.001,r_penalty_factor=1e10):
    penalty_r_dict = {"tip_outer": r_penalty_factor, "tip_outlet": 0, "barrel_outer": 0,
                      "plunger": 0
                      }
    penalty_r = CoefficientFunction([penalty_r_dict[bound] for bound in mesh.GetBoundaries()])
    v = fes.TestFunction()
    return penalty_r * (
            normal_vector_outlet(u_last_step, torus_position_z=torus_position_z,
                                 radius_ring=radius_ring) *
                   get_target_displacement_on_tip_outlet(u_last_step,
                                 torus_position_z=torus_position_z,radius_ring=radius_ring,
                                 radius_torus=radius_torus,radius_outlet=radius_ring-radius_torus)) * (
                   normal_vector_outlet(u_last_step, torus_position_z=torus_position_z,
                                        radius_ring=radius_ring) * v) * ds


def contribution_to_linear_form_shape_tip_fix_coordinates(fes,mesh,u_last_step,torus_position_z=-0.001,radius_ring=0.0015,
                                          radius_torus=0.001,r_penalty_factor=1e10):
    penalty_r_dict = {"tip_outer": r_penalty_factor, "tip_outlet": 0, "barrel_outer": 0,
                      "plunger": 0
                      }
    penalty_r = CoefficientFunction([penalty_r_dict[bound] for bound in mesh.GetBoundaries()])
    v = fes.TestFunction()
    tg=get_target_displacement_on_tip_outlet(u_last_step,
                                 torus_position_z=torus_position_z,radius_ring=radius_ring,
                                 radius_torus=radius_torus,radius_outlet=radius_ring-radius_torus)
    return penalty_r * (tg[0]*v[0]+tg[1]*v[1]+tg[2]*v[2]) * ds

def contribution_to_linear_form_xy_outlet(fes,mesh,u_last_step,torus_position_z=-0.001,radius_ring=0.0015,
                                          radius_torus=0.001,r_penalty_factor=1e10):
    penalty_r_dict = {"tip_outer": r_penalty_factor, "tip_outlet": 0, "barrel_outer": 0,
                      "plunger": 0
                      }
    penalty_r = CoefficientFunction([penalty_r_dict[bound] for bound in mesh.GetBoundaries()])
    v = fes.TestFunction()
    below_tip = IfPos(torus_position_z - (z + u_last_step[2]), 1, 0)
    # The desired displacement, in plane, is such that a radius of radius_ring-radius_torus is reached, from the original x y position
    desirable_displacement_x = IfPos(x * x + y * y, x * (1 - (radius_ring - radius_torus) / sqrt(x * x + y * y)), 0)
    desirable_displacement_y = IfPos(x * x + y * y, y * (1 - (radius_ring - radius_torus) / sqrt(x * x + y * y)), 0)
    return penalty_r * below_tip * desirable_displacement_x * v[
        0] * ds + penalty_r * below_tip * desirable_displacement_y * v[1] * ds

def contribution_to_bilinear_form_xy_outlet(fes,mesh,u_last_step,torus_position_z=-0.001,
                                          r_penalty_factor=1e10):
    penalty_r_dict = {"tip_outer": r_penalty_factor, "tip_outlet": 0, "barrel_outer": 0,
                      "plunger": 0
                      }
    penalty_r = CoefficientFunction([penalty_r_dict[bound] for bound in mesh.GetBoundaries()])
    v = fes.TestFunction()
    u = fes.TrialFunction()
    below_tip = IfPos(torus_position_z - (z + u_last_step[2]), 1, 0)
    return penalty_r * below_tip  * u[0] * v[
        0] * ds + penalty_r * below_tip * u[1] * v[1] * ds




def run_simulation_step(fes,mesh,E,nu,volumeChange,plastic_deformation,u_last_step,
                        torus_position_z=-0.001,radius_ring=0.0015,
                        radius_torus=0.001,
                        plunger_z_penalty_factor=1e10,xy_plunger_penalty_factor=1e10,r_penalty_factor_barrel=1e10,
                        plunger_displacement=0,r_penalty_factor=1e10,fix_coordinates_on_outlet=False,
                        increase_penalty_fix_coordinates=1):
    a_incompressible=get_basic_bilinear_form(fes,mesh,
                                             E,nu=nu,plunger_z_penalty_factor=plunger_z_penalty_factor,
                                             xy_plunger_penalty_factor=xy_plunger_penalty_factor,
                                             r_penalty_factor_barrel=r_penalty_factor_barrel)
    if fix_coordinates_on_outlet:
        a_incompressible += contribution_to_bilinear_form_shape_tip_fix_coordinates(fes,mesh,
                                    r_penalty_factor=r_penalty_factor*increase_penalty_fix_coordinates)
    a_incompressible += contribution_to_bilinear_form_shape_tip(fes,mesh,u_last_step,
            torus_position_z=torus_position_z,radius_ring=radius_ring,r_penalty_factor=r_penalty_factor)
    a_incompressible += contribution_to_bilinear_form_xy_outlet(fes,mesh,u_last_step,
        torus_position_z=torus_position_z,r_penalty_factor=r_penalty_factor)
    c = MultiGridPreconditioner(a_incompressible)
    a_incompressible.Assemble()
    f=get_basic_linear_form(fes,mesh,plunger_z_penalty_factor=plunger_z_penalty_factor,plunger_displacement=plunger_displacement)
    f+=getVolumetricChangeContributionForLinearForm(fes,E,nu,volumeChange)
    f+=getPlasticDeformationContributionForLinearForm(fes,E,nu,plastic_deformation)
    f+=contribution_to_linear_form_shape_tip(fes,mesh,u_last_step,
        torus_position_z=torus_position_z,radius_ring=radius_ring,
        radius_torus=radius_torus,r_penalty_factor=r_penalty_factor)
    if fix_coordinates_on_outlet:
        f += contribution_to_linear_form_shape_tip_fix_coordinates(fes,mesh,u_last_step,
                                        torus_position_z=torus_position_z,radius_ring=radius_ring,
                                        radius_torus=radius_torus,
                                        r_penalty_factor=r_penalty_factor*increase_penalty_fix_coordinates)
    f += contribution_to_linear_form_shape_tip(fes, mesh, u_last_step,
                                                   torus_position_z=torus_position_z, radius_ring=radius_ring,
                                                   radius_torus=radius_torus, r_penalty_factor=r_penalty_factor)
    f+=contribution_to_linear_form_xy_outlet(fes,mesh,u_last_step,
                                          torus_position_z=torus_position_z,radius_ring=radius_ring,
                                          radius_torus=radius_torus,r_penalty_factor=r_penalty_factor)
    f.Assemble()
    u = GridFunction(fes)
    BVP(bf=a_incompressible, lf=f, gf=u, pre=c)
    return(u)



def getEvacuationTimeConstant(mesh,K,E,nu,p_offset=300):
    fesx = H1(mesh, order=1, dim=1)
    kappa = E*3/(1-2*nu)
    time_constants=GridFunction(fesx)
    ind=0
    for v in mesh.vertices:
        fdiv= GridFunction(fesx) # this is our divergence for each element
        fdiv.vec.data[ind]=1
        a = BilinearForm(fesx, symmetric=True, condense=True)
        v = fesx.TestFunction()
        p = fesx.TrialFunction()
        a += Grad(v)*Grad(p)*dx
        f = LinearForm(fesx)
        f += v*K*fdiv*dx
        penalty_dict = {"tip_outer": 0, "tip_outlet": 1e6, "barrel_outer": 0,
                      "plunger": 0,
                      "barrel_front": 0,
                      "transition": 0
                      }
        penalty_dict_coeff = CoefficientFunction([penalty_dict[bound] for bound in mesh.GetBoundaries()])
        a+=p*v*penalty_dict_coeff*ds
        f += p_offset * v * penalty_dict_coeff * ds
        c = MultiGridPreconditioner(a)
        a.Assemble()
        f.Assemble()
        pc = GridFunction(fesx)
        BVP(bf=a, lf=f, gf=pc, pre=c)
        time_constants.vec.data[ind]=pc.vec.data[ind]/kappa
        ind=ind+1
        #return pc
    return time_constants
    #return pc


def get_simulation(mesh_filename,E=2e6, nu=0,nu_i=0.4975,sd=0.1,K=1e5,dt_max=0.01,eta_yield=1e6,n_steps=30,
                   torus_position_z=-0.001,radius_ring=0.0015,
                   radius_torus=0.001,
                   plunger_z_penalty_factor=1e10,xy_plunger_penalty_factor=1e10,r_penalty_factor_barrel=1e10,
                   plunger_displacement_initial=0,
                   plunger_displacement_rate_m_sec=-1e-4,
                   max_plunger_displacement = -5e-4,r_penalty_factor=1e10,safety_factor=1,
                   volume_change_max=0.8,yield_strain=0.2,p_offset=300,
                   increase_penalty_fix_coordinates=1):
    mesh = loadMesh(mesh_filename)
    evacuationTimeConstant=getEvacuationTimeConstant(mesh, K, E, nu,p_offset=p_offset)
    mean_fluid_evacuation_time = Integrate(evacuationTimeConstant,mesh)/Integrate(1,mesh)
    print("Mean fluid evacuation time constant = ",mean_fluid_evacuation_time,"\n")
    Draw(evacuationTimeConstant,mesh,"evacuationTimeConstant")
    fes = H1(mesh, order=2, dim=3)
    fesx = H1(mesh, order=1, dim=1)
    volumeChange = GridFunction(fesx)
    volumeChange.Set(0)
    fesstress = H1(mesh, order=1, dim=9)
    plastic_deformation = GridFunction(fesstress)
    plastic_deformation.Set((0, 0, 0, 0, 0, 0, 0, 0, 0))
    return_dict = {
        "mesh": mesh,
        "evacuationTimeConstant":evacuationTimeConstant,
        "u": [],
        "stress": [],
        "strain": [],
        "delta_plastic": [],
        "p": [],
        "V": [],
        "plasticDeformation": [],
        "t": [],
        "plunger_position": []
    }
    t = 0
    u_old = GridFunction(fes)
    plunger_displacement=plunger_displacement_initial
    for ind in range(n_steps):
        if ind>0:
            u_old.Set(u)
        u = run_simulation_step(fes, mesh, extrapolate_E(E, nu), nu_i, volumeChange, plastic_deformation, u_old,
                                torus_position_z=torus_position_z, radius_ring=radius_ring,
                                radius_torus=radius_torus,
                                plunger_z_penalty_factor=plunger_z_penalty_factor,
                                xy_plunger_penalty_factor=xy_plunger_penalty_factor,
                                r_penalty_factor_barrel=r_penalty_factor_barrel,
                                plunger_displacement=plunger_displacement, r_penalty_factor=r_penalty_factor,
                                fix_coordinates_on_outlet=True,
                                increase_penalty_fix_coordinates=increase_penalty_fix_coordinates
                                )
        u=getSmoothenedNormalizedAxisymmetricZDeformed(mesh, u,u_old, sd)
        u_old.Set(u) # Repetition of the smoothing and simulation step to obtain better agreement
        # between u_old and u and thus better finite geometrical displacements for the tip_outer boundary domain
        u = run_simulation_step(fes, mesh, extrapolate_E(E, nu), nu_i, volumeChange, plastic_deformation, u_old,
                                torus_position_z=torus_position_z, radius_ring=radius_ring,
                                radius_torus=radius_torus,
                                plunger_z_penalty_factor=plunger_z_penalty_factor,
                                xy_plunger_penalty_factor=xy_plunger_penalty_factor,
                                r_penalty_factor_barrel=r_penalty_factor_barrel,
                                plunger_displacement=plunger_displacement, r_penalty_factor=r_penalty_factor,
                                fix_coordinates_on_outlet=False,
                                increase_penalty_fix_coordinates=increase_penalty_fix_coordinates)
        u = getSmoothenedNormalizedAxisymmetric(mesh, u, sd)
        u_scaffold=run_simulation_step(fes, mesh, E, nu, GridFunction(fesx), plastic_deformation, u_old,
                                       torus_position_z=torus_position_z, radius_ring=radius_ring,
                                       radius_torus=radius_torus,
                                       plunger_z_penalty_factor=plunger_z_penalty_factor,
                                       xy_plunger_penalty_factor=xy_plunger_penalty_factor,
                                       plunger_displacement=plunger_displacement, r_penalty_factor=r_penalty_factor,
                                       fix_coordinates_on_outlet = False,
                                       increase_penalty_fix_coordinates=increase_penalty_fix_coordinates)
        u_scaffold = getSmoothenedNormalizedAxisymmetric(mesh, u_scaffold, sd)
        stress = sigma(epsilon(u), extrapolate_E(E, nu), nu_i)  # -sigma(epsilon(gg1), E, nu)
        stress_deformatory = stress - Trace(stress) * Id(3) / 3
        strain_deformatory = (1 + nu) / E * stress_deformatory
        Draw(strain_deformatory, mesh, "strain_deformatory")
        plastic_deformation_delta_all = strain_deformatory - plastic_deformation
        plastic_deformation_delta=GridFunction(fesstress)
        plastic_deformation_delta.Set(get_deformation_with_yield_strain(plastic_deformation_delta_all,yield_strain))
        max_plastic_stress = getAbsMax(getVanMises(sigma(plastic_deformation_delta,extrapolate_E(E, nu), nu_i)) , mesh)
        p_fluid = GridFunction(fesx) # We consider the fluid pressure to be the excess isometric strain develoopped under
        # the incompressible material parameters minus the one that can be sustained by the solid material itself under
        # the same deformation. Signs: if u is compressive, trace(epsilon(u)) is negative, and sigma also; the incompressible
        # sigma is much larger in absolute magnitude and so the stress contribution will be positive
        # This should cause a negative volume shift over time that offsets the positive pressure, and indeed, for negative volumeChange,
        # sigmaVolumetric is negative and so that plays out
        # at equilibrium, p=0 and the fluid loss offsets the difference between compressive and incompressible material scenarios as desired
        p_fluid.Set(Trace(
            sigma(epsilon(u), E, nu)-sigma(epsilon(u),extrapolate_E(E, nu),nu_i)+sigmaVolumetric(volumeChange, extrapolate_E(E, nu), nu_i, 3))/3)
        excessVolume=GridFunction(fesx)
        # For continuity, we limit here the volume change in a single step to 0.5, more corresponds to densification which would prevent further rapid volume change
        excessVolume.Set(Trace(-epsilon(u_scaffold)+epsilon(u)))
        Draw(excessVolume,mesh,"excess_Volume")
        if max_plastic_stress > 0:
            dt_yield = safety_factor * eta_yield / max_plastic_stress
        else:
            dt_yield=dt_max
        print("Plastic dt ", dt_yield)
        dt=dt_yield
        if dt>mean_fluid_evacuation_time:
            dt=mean_fluid_evacuation_time
        if (dt > dt_max):
            dt = dt_max
        if t>5: # For this slow speed, otherwise, it takes too many simulations
            mean_fluid_evacuation_time=0.1
        if t>10:
            mean_fluid_evacuation_time = 0.2
        changeExcessVolume = -excessVolume * (CoefficientFunction(1) - exp(-dt / IfPos(evacuationTimeConstant,evacuationTimeConstant,1e-6)))
        changeExcessVolumeMax = getAbsMax(changeExcessVolume,mesh)
        print("Max projected volume change ",changeExcessVolumeMax)
        cf_fluid=1
        if changeExcessVolumeMax>volume_change_max:
            dt=dt*volume_change_max/changeExcessVolumeMax
            cf_fluid=volume_change_max/changeExcessVolumeMax
        plastic_deformation_new = GridFunction(fesstress)
        plastic_deformation_new.Set(plastic_deformation + safety_factor * plastic_deformation_delta*dt/dt_yield)
        plastic_deformation = plastic_deformation_new
        volumeChangeNew = GridFunction(fesx)
        volumeChangeNew.Set(volumeChange + changeExcessVolume*cf_fluid)
        Draw(plastic_deformation, mesh, "plastic_deformation")
        volumeChange = getSmoothenedNormalizedAxisymmetricZDeformed(mesh,volumeChangeNew,u,sd)
        Draw(volumeChange, mesh, "volume_change")
        return_dict["t"].append(t)
        return_dict["plunger_position"].append(plunger_displacement)
        t = t + dt
        plunger_displacement=plunger_displacement+plunger_displacement_rate_m_sec*dt
        if plunger_displacement<max_plunger_displacement: # Plunger is going down, this is for negative numbers
            plunger_displacement=max_plunger_displacement
        print("current dt", dt, "\n")
        Draw(p_fluid,mesh,"p")
        return_dict["stress"].append(stress)
        return_dict["strain"].append(epsilon(u))
        return_dict["V"].append(volumeChange)
        return_dict["plasticDeformation"].append(plastic_deformation)
        return_dict["delta_plastic"].append(plastic_deformation_delta)
        return_dict["u"].append(u)
        return_dict["p"].append(p_fluid)
        Draw(u, mesh, "u")
        print(ind, " of ", n_steps, " simulation steps done. Present piston position = ", plunger_displacement)
    return return_dict




