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






def loadMesh():
    m = nm.Mesh(3)
    m.Load("/Users/unige/Documents/publications/injectability/1mL_minimal.vol")
    m.SetBCName(0, "plunger")
    m.SetBCName(1, "barrel_outer")
    m.SetBCName(2, "tip_outer")
    m.SetBCName(3, "tip_outlet")
    mesh = Mesh(m)
    return mesh




def getMax(u,mesh):
    output=[]
    started=False
    for v in mesh.vertices:
        pointvalue = u(mesh(v.point[0],v.point[1],v.point[2]))
        if not started:
            started=True
            for ind in range(u.dim):
                output.append(pointvalue[ind])
        for ind in range(u.dim):
            if(pointvalue[ind] > output[ind]):
                output[ind]=pointvalue[ind]
    return output


def getAbsMax(u,mesh):
    if u.dim>1:
        output=[]
        started=False
        for v in mesh.vertices:
            pointvalue = u(mesh(v.point[0],v.point[1],v.point[2]))
            if not started:
                started=True
                for ind in range(u.dim):
                    output.append(abs(pointvalue[ind]))
            for ind in range(u.dim):
                if(abs(pointvalue[ind]) > output[ind]):
                    output[ind]=abs(pointvalue[ind])
    else:
        output=0
        started=False
        for v in mesh.vertices:
            pointvalue = u(mesh(v.point[0],v.point[1],v.point[2]))
            if not started:
                started=True
                output=abs(pointvalue)
            if(abs(pointvalue) > output):
                output=abs(pointvalue)
    return output







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

def isotropic_pressure(u,E,nu_high,nu_low,volumeChange):
    return Trace(
            sigma(epsilon(u), E, nu_low) + sigmaVolumetric(volumeChange, E, nu_high, 3) - sigma(epsilon(u), E, nu_high)) / 3


def sigmaVolumetric(volumeChange,E,nu,dim=3):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    return lam*volumeChange*Id(3)

def getVolumetricChangeContributionForLinearForm(fes,E,nu,targetChange):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    v = fes.TestFunction()
    e=epsilon(v)
    s=lam*targetChange*Id(fes.dim) # Stress associated with the volumetric change
    return (s[0]*e[0]+s[4]*e[4]+s[8]*e[8])*dx # only the volumetric part is of interest


def getPlasticDeformationContributionForLinearForm(fes,E,nu,targetChange):
    v = fes.TestFunction()
    e=epsilon(v)
    s=sigma(targetChange,E,nu) # Stress associated with the volumetric change
    return (s[0]*e[0]+s[1]*e[1]+s[2]*e[2]+s[3]*e[3]+s[4]*e[4]+s[5]*e[5]+s[6]*e[6]+s[7]*e[7]+s[8]*e[8])*dx

def getElasticBilinearForm(fes,E,nu):
    a = BilinearForm(fes, symmetric=True, condense=True)
    v = fes.TestFunction()
    u = fes.TrialFunction()
    e=epsilon(v)
    s=sigma(epsilon(u),E,nu)
    a+=(s[0]*e[0]+s[1]*e[1]+s[2]*e[2]+s[3]*e[3]+s[4]*e[4]+s[5]*e[5]+s[6]*e[6]+s[7]*e[7]+s[8]*e[8])*dx
    return a


# toSmooth is a grid function defined on the mesh
def getSmoothened(mesh,toSmooth,sd):
    fes=H1(mesh, order=1, dim=toSmooth.dim)
    output=GridFunction(fes)
    # loop over the vertices of the mesh, in order 1, those are also the data elements of the grid function
    ind=0
    for v in mesh.vertices:
        posx=v.point[0]
        posy=v.point[1]
        posz=v.point[2]
        sq=math.sqrt(2*math.pi)*math.sqrt(2*math.pi)*math.sqrt(2*math.pi)*sd*sd*sd
        output.vec.data[ind]=Integrate(exp(-((x-posx)*(x-posx)+(y-posy)*(y-posy)+(z-posz)*(z-posz))/2/sd/sd)*toSmooth/sq,mesh)
        ind=ind+1
    return output



# toSmooth is a grid function defined on the mesh
def getSmoothenedNormalized(mesh,toSmooth,sd):
    fes=H1(mesh, order=1, dim=toSmooth.dim)
    input=GridFunction(fes)
    input.Set(toSmooth)
    output=GridFunction(fes)
    # loop over the vertices of the mesh, in order 1, those are also the data elements of the grid function
    ind=0
    for v in mesh.vertices:
        posx=v.point[0]
        posy=v.point[1]
        posz=v.point[2]
        normalization_factor = Integrate(exp(-((x-posx)*(x-posx)+(y-posy)*(y-posy)+(z-posz)*(z-posz))/2/sd/sd),mesh)
        if input.dim==1:
            output.vec.data[ind]=Integrate(exp(-((x-posx)*(x-posx)+(y-posy)*(y-posy)+(z-posz)*(z-posz))/2/sd/sd)*input,mesh)/normalization_factor
        else:
            for ind2 in range(input.dim):
                output.vec.data[ind][ind2] = Integrate(exp(-(
                            (x - posx) * (x - posx) + (y - posy) * (y - posy) + (z - posz) * (
                                z - posz)) / 2 / sd / sd) * input[ind2], mesh) / normalization_factor
        ind=ind+1
        if ind % 100 == 0:
            print(ind, " of ", len(output.vec.data))
    return output


# toSmooth is a grid function defined on the mesh
def getSmoothenedNormalizedAxisymmetric(mesh,toSmooth,sd):
    fes=H1(mesh, order=1, dim=toSmooth.dim)
    input=GridFunction(fes)
    input.Set(toSmooth)
    output=GridFunction(fes)
    # loop over the vertices of the mesh, in order 1, those are also the data elements of the grid function
    ind=0
    for v in mesh.vertices:
        posx=v.point[0]
        posy=v.point[1]
        posz=v.point[2]
        r2=posx*posx+posy*posy
        r=sqrt(r2)
        radial_distance_difference2=IfPos(r2-x*x-y*y,r2-x*x-y*y,-r2+x*x+y*y)
        integration_factor_z = exp(
            -(radial_distance_difference2 + (z - posz) * (z - posz)) / 2 / sd / sd)
        normalization_factor_z = Integrate(
            integration_factor_z,
            mesh)
        if input.dim==1:
            output.vec.data[ind]=Integrate(integration_factor_z*input,mesh)/normalization_factor_z
        else:
            integration_factor_x = IfPos(r, x / r * exp(
                -(radial_distance_difference2 + (z - posz) * (z - posz)) / 2 / sd / sd), CoefficientFunction(0))
            normalization_factor_x = Integrate(
                IfPos(r, x * x / r / r * exp(-(radial_distance_difference2 + (z - posz) * (z - posz)) / 2 / sd / sd),
                      CoefficientFunction(0)),
                mesh)
            integration_factor_y = IfPos(r, y / r * exp(
                -(radial_distance_difference2 + (z - posz) * (z - posz)) / 2 / sd / sd), CoefficientFunction(0))
            normalization_factor_y = Integrate(
                IfPos(r, y * y / r / r * exp(-(radial_distance_difference2 + (z - posz) * (z - posz)) / 2 / sd / sd),
                      CoefficientFunction(0)),
                mesh)
            sx = Integrate(integration_factor_x*input[0],mesh)/normalization_factor_x
            sy = Integrate(integration_factor_y * input[1], mesh) / normalization_factor_y
            # Static situation without rotation around axis, so xy-displacements have to be around radial direction
            if r2==0:
                sx=0
                sy=0
            else:
                if sx > 0:
                    sx_new=posx/r*sqrt(sx*sx+sy*sy)
                    sy_new=posy/r*sqrt(sx*sx+sy*sy)
                else:
                    sx_new = -posx / r * sqrt(sx * sx + sy * sy)
                    sy_new = -posy / r * sqrt(sx * sx + sy * sy)
                sx=sx_new
                sy=sy_new
            output.vec.data[ind][0]=sx
            output.vec.data[ind][1]=sy
            output.vec.data[ind][2] = Integrate(integration_factor_z * input[2], mesh) / normalization_factor_z
        ind=ind+1
        if ind % 100 == 0:
            print(ind, " of ", len(output.vec.data))
    return output


# toSmooth is a grid function defined on the mesh
def getSmoothenedStressNormalizedAxisymmetric(mesh,toSmooth,sd):
    # How does this transform?
    # Stress tensor: there cannot be any rotatory component as we apply no rotatory boundary conditions.
    # we know that utheta=0, also that dur/dtheta=0, duz/dtheta=0
    # and so in the end we only have the components err, erz, and ezz
    # where we further have: err=exx*x^2/r^2+eyy*y^2/r^2; erz=ezx*x/r+ezy*y/r
    fes=H1(mesh, order=1, dim=toSmooth.dim)
    input=GridFunction(fes)
    input.Set(toSmooth)
    output=GridFunction(fes)
    # loop over the vertices of the mesh, in order 1, those are also the data elements of the grid function
    ind=0
    for v in mesh.vertices:
        posx=v.point[0]
        posy=v.point[1]
        posz=v.point[2]
        r2=posx*posx+posy*posy
        r=sqrt(r2)
        radial_distance_difference2 = IfPos(r2 - x * x - y * y, r2 - x * x - y * y, -r2 + x * x + y * y)
        integration_factor_xx = IfPos(r, x*x / r/r * exp(
            -(radial_distance_difference2 + (z - posz) * (z - posz)) / 2 / sd / sd), CoefficientFunction(0))
        normalization_factor_xx = Integrate(
            integration_factor_xx,mesh)
        err=Integrate(toSmooth[0]*integration_factor_xx)/normalization_factor_xx
    return output


# toSmooth is a grid function defined on the mesh
def getSmoothenedGradient(mesh,toSmooth,sd):
    fes=H1(mesh, order=3, dim=3)
    output=GridFunction(fes)
    # loop over the vertices of the mesh, in order 1, those are also the data elements of the grid function
    ind=0
    for v in mesh.vertices:
        posx=v.point[0]
        posy=v.point[1]
        posz=v.point[2]
        sq=math.sqrt(2*math.pi)*math.sqrt(2*math.pi)*math.sqrt(2*math.pi)*sd*sd*sd*sd*sd
        output.vec.data[ind][0]=Integrate((x-posx)*exp(-((x-posx)*(x-posx)+(y-posy)*(y-posy)+(z-posz)*(z-posz))/2/sd/sd)*toSmooth/sq,mesh)
        output.vec.data[ind][1]=Integrate((y-posy)*exp(-((x-posx)*(x-posx)+(y-posy)*(y-posy)+(z-posz)*(z-posz))/2/sd/sd)*toSmooth/sq,mesh)
        output.vec.data[ind][2]=Integrate((z-posz)*exp(-((x-posx)*(x-posx)+(y-posy)*(y-posy)+(z-posz)*(z-posz))/2/sd/sd)*toSmooth/sq,mesh)
        ind=ind+1
    return output


# toSmooth is a grid function defined on the mesh
def getSmoothenedLaplacian(mesh,toSmooth,sd):
    fes=H1(mesh, order=3, dim=toSmooth.dim)
    output=GridFunction(fes)
    # loop over the vertices of the mesh, in order 1, those are also the data elements of the grid function
    ind=0
    for v in mesh.vertices:
        posx=v.point[0]
        posy=v.point[1]
        posz=v.point[2]
        sq=math.sqrt(2*math.pi)*math.sqrt(2*math.pi)*math.sqrt(2*math.pi)*sd*sd*sd*sd*sd
        output.vec.data[ind]=Integrate((((x-posx)*(x-posx)+(y-posy)*(y-posy)+(z-posz)*(z-posz))/sd/sd-3)*exp(-((x-posx)*(x-posx)+(y-posy)*(y-posy)+(z-posz)*(z-posz))/2/sd/sd)*toSmooth/sq,mesh)
        ind=ind+1
    return output



def get_basic_bilinear_form(fes,mesh,E,nu=0.49):
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



def get_basic_linear_form(fes,mesh):
    penalty_z_dict_linear = {"plunger": plunger_z_penalty_factor*plunger_displacement,"barrel_outer": 0,"tip_outer": 0, "tip_outlet": 0
                      }
    penalty_z_linear = CoefficientFunction([penalty_z_dict_linear[bound] for bound in mesh.GetBoundaries()])
    f = LinearForm(fes)
    v = fes.TestFunction()
    f += penalty_z_linear*v[2]*ds
    return(f)


def xy_contribution_at_torus_to_linear_form(fes,mesh,u_old):
    applicable_boundaries_dict = {"plunger": 0,"barrel_outer": 0,"tip_outer": penalty_xy_torus, "tip_outlet": 0}
    applicable_boundaries = CoefficientFunction([applicable_boundaries_dict[bound] for bound in mesh.GetBoundaries()])
    v = fes.TestFunction()
    below_barrel = IfPos(torus_position_z+radius_torus - (z + u_old[2]), 1, 0)
    above_tubing = IfPos(z+u_old[2]-torus_position_z,1,0)
    inside_torus = IfPos(radius_ring-sqrt((x+u_old[0])*(x+u_old[0])+(y+u_old[1])*(y+u_old[1])),1,0)
    applicable = applicable_boundaries*below_barrel*above_tubing*inside_torus
    target_radius=IfPos(radius_torus*radius_torus - (z+u_old[2]-torus_position_z)*(z+u_old[2]-torus_position_z),
                        radius_ring - sqrt(radius_torus*radius_torus - (z+u_old[2]-torus_position_z)*(z+u_old[2]-torus_position_z)),
                        radius_ring)
    initial_radius = sqrt(x*x+y*y)
    desirable_displacement_x = IfPos(initial_radius,(target_radius-initial_radius)*x/initial_radius,0)
    desirable_displacement_y = IfPos(initial_radius,(target_radius-initial_radius)*y/initial_radius,0)
    return applicable * desirable_displacement_x * v[
        0] * ds + applicable * desirable_displacement_y * v[1] * ds


def xy_contribution_at_torus_to_bilinear_form(fes,mesh,u_old):
    applicable_boundaries_dict = {"plunger": 0,"barrel_outer": 0,"tip_outer": penalty_xy_torus, "tip_outlet": 0}
    applicable_boundaries = CoefficientFunction([applicable_boundaries_dict[bound] for bound in mesh.GetBoundaries()])
    v = fes.TestFunction()
    u = fes.TrialFunction()
    below_barrel = IfPos(torus_position_z+radius_torus - (z + u_old[2]), 1, 0)
    above_tubing = IfPos(z+u_old[2]-torus_position_z,1,0)
    inside_torus = IfPos(radius_ring-sqrt((x+u_old[0])*(x+u_old[0])+(y+u_old[1])*(y+u_old[1])),1,0)
    applicable = applicable_boundaries*below_barrel*above_tubing*inside_torus
    return applicable * u[0] * v[
        0] * ds + applicable * u[1] * v[1] * ds


def xy_contribution_at_outlet_to_bilinear_form(fes,mesh,u_old):
    applicable_boundaries_dict = {"plunger": 0,"barrel_outer": 0,"tip_outer": penalty_tubing, "tip_outlet": 0}
    applicable_boundaries = CoefficientFunction([applicable_boundaries_dict[bound] for bound in mesh.GetBoundaries()])
    v = fes.TestFunction()
    u = fes.TrialFunction()
    below_tip = IfPos(torus_position_z-(z+u_old[2]),1,0)
    above_tubing_end = IfPos(z+u_old[2]-lower_end_outlet,1,0)
    return applicable_boundaries*below_tip*above_tubing_end * u[0] * v[0] * ds + applicable_boundaries *below_tip*above_tubing_end* u[1] * v[1] * ds


def xy_contribution_at_outlet_to_linear_form(fes,mesh,u_old):
    applicable_boundaries_dict = {"plunger": 0,"barrel_outer": 0,"tip_outer": penalty_tubing, "tip_outlet": 0}
    applicable_boundaries = CoefficientFunction([applicable_boundaries_dict[bound] for bound in mesh.GetBoundaries()])
    v = fes.TestFunction()
    u = fes.TrialFunction()
    below_tip = IfPos(torus_position_z - (z + u_old[2]), 1, 0)
    above_tubing_end = IfPos(z + u_old[2] - lower_end_outlet, 1, 0)
    # The desired displacement, in plane, is such that a radius of radius_ring-radius_torus is reached, from the original x y position
    desirable_displacement_x = IfPos(x*x+y*y,x*(1-(radius_ring-radius_torus)/sqrt(x*x+y*y)),0)
    desirable_displacement_y = IfPos(x * x + y * y, y * (1 - (radius_ring - radius_torus) / sqrt(x * x + y * y)), 0)
    return applicable_boundaries*below_tip*above_tubing_end * desirable_displacement_x * v[0] * ds + applicable_boundaries *below_tip*above_tubing_end* desirable_displacement_y * v[1] * ds






def normal_vector_torus_segment_in_situ():
    x_f = CoefficientFunction(x )
    y_f = CoefficientFunction(y )
    z_f = CoefficientFunction(z )
    position_vector = CoefficientFunction((x ,y ,z ))
    r_xy = sqrt(x_f*x_f+y_f*y_f)
    torus_center = IfPos(
        r_xy,CoefficientFunction((x_f/r_xy*radius_ring, y_f/r_xy*radius_ring, torus_position_z)),CoefficientFunction((radius_ring, 0, torus_position_z)))
    return (position_vector-torus_center)/sqrt(
        (position_vector[0]-torus_center[0])*(position_vector[0]-torus_center[0])+
        (position_vector[1] - torus_center[1]) * (position_vector[1] - torus_center[1]) +
        (position_vector[2] - torus_center[2]) * (position_vector[2] - torus_center[2]))

def stress_surface_force_linear_form_contribution(stress_deformatory,mesh,fes):
    v = fes.TestFunction()
    slip_dict = {"tip_outer": 1, "tip_outlet": 0, "barrel_outer": 1, "plunger": 0,
                      "barrel_front": 1, "transition": 1}
    slip = CoefficientFunction([slip_dict[bound] for bound in mesh.GetBoundaries()])
    return slip*stress_surface_force_tangential(stress_deformatory,mesh)*v*ds


def penalty_large_scale_deformation_rounded_segments(u):
    # The idea here is that additional normal displacement is added so that the final positions do not go outside the radius
    # geometric positions
    # this is the torus having been cut out to permit the transition between the barrel and the needle without singulariy
    x_f=CoefficientFunction(x+u[0])
    y_f = CoefficientFunction(y + u[1])
    z_f = CoefficientFunction(z + u[2])
    d_torus_center=CoefficientFunction(sqrt((sqrt(x_f*x_f+y_f*y_f)-radius_ring)*(sqrt(x_f*x_f+y_f*y_f)-radius_ring)+(z_f-torus_position_z)*(z_f-torus_position_z)))
    d_torus_inside=IfPos(radius_torus-d_torus_center,radius_torus-d_torus_center,0)
    return d_torus_inside

def getFiniteCurvatureCorrectionTorus_bilinear_form_contribution(fes,mesh):
    penalty_r_dict = {"tip_outer": r_penalty_factor, "tip_outlet": 0, "barrel_outer": 0,
                      "plunger": 0,
                      "barrel_front": 0, "transition": r_penalty_factor
                      }
    penalty_r = CoefficientFunction([penalty_r_dict[bound] for bound in mesh.GetBoundaries()])
    v = fes.TestFunction()
    u = fes.TestFunction()
    return penalty_r * penalty_large_scale_deformation_rounded_segments(u) * (specialcf.normal(3) * v) * ds

def normal_vector_lower_surface(u):
    applicable_dict = {"plunger": 0,"barrel_outer": 0,"tip_outer": 1, "tip_outlet": 0}
    applicable = CoefficientFunction([applicable_dict[bound] for bound in mesh.GetBoundaries()])
    torus_normal=normal_vector_torus_segment(u)
    x_f = CoefficientFunction(x + u[0])
    y_f = CoefficientFunction(y + u[1])
    z_f = CoefficientFunction(z + u[2])
    outwardRadial=IfPos(x*x+y*y,CoefficientFunction((-x/sqrt(x*x+y*y),-y/sqrt(x*x+y*y),0)),CoefficientFunction((0,0,0)))
    return IfPos(torus_position_z-z_f,outwardRadial,IfPos(z_f-(torus_position_z+radius_torus),CoefficientFunction((0,0,1)),torus_normal))

def normal_vector_torus_segment(u):
    x_f = CoefficientFunction(x + u[0])
    y_f = CoefficientFunction(y + u[1])
    z_f = CoefficientFunction(z + u[2])
    position_vector = CoefficientFunction((x + u[0],y + u[1],z + u[2]))
    r_xy = sqrt(x_f*x_f+y_f*y_f)
    torus_center = IfPos(
        r_xy,CoefficientFunction((x_f/r_xy*radius_ring, y_f/r_xy*radius_ring, torus_position_z)),CoefficientFunction((radius_ring, 0, torus_position_z)))
    basic_normvector=(position_vector-torus_center)/sqrt(
        (position_vector[0]-torus_center[0])*(position_vector[0]-torus_center[0])+
        (position_vector[1] - torus_center[1]) * (position_vector[1] - torus_center[1]) +
        (position_vector[2] - torus_center[2]) * (position_vector[2] - torus_center[2]))
    return IfPos(z_f-torus_position_z,basic_normvector,CoefficientFunction((0, 0, 0)))



#    torus_position_z = -0.05
#    radius_torus = 0.05
#    radius_ring = 0.1

#    radius_main_barrel = 0.15
#    z_lower_main_barrel = 0
#    z_upper_main_barrel = 0.2


def run_simulation_step(fes,mesh,E,nu,volumeChange,plastic_deformation,u_last_step,stress_deformatory_old):
    a_incompressible=get_basic_bilinear_form(fes,mesh,E,nu=nu)
    penalty_r_dict = {"tip_outer": r_penalty_factor, "tip_outlet": 0, "barrel_outer": 0,
                      "plunger": 0
                      }
    v = fes.TestFunction()
    u = fes.TrialFunction()
    penalty_r = CoefficientFunction([penalty_r_dict[bound] for bound in mesh.GetBoundaries()])
    a_incompressible += penalty_r * (normal_vector_torus_segment(u_last_step) * u) * (
            normal_vector_torus_segment(u_last_step) * v) * ds
    #a_incompressible += xy_contribution_at_torus_to_bilinear_form(fes, mesh, u_last_step)
    c = MultiGridPreconditioner(a_incompressible)
    a_incompressible.Assemble()
    f=get_basic_linear_form(fes,mesh)
    f+=getVolumetricChangeContributionForLinearForm(fes,E,nu,volumeChange)
    f+=getPlasticDeformationContributionForLinearForm(fes,E,nu,plastic_deformation)
    #f+=xy_contribution_at_outlet_to_linear_form(fes,mesh,u_last_step)
    #f+=xy_contribution_at_torus_to_linear_form(fes, mesh, u_last_step)
    f.Assemble()
    u = GridFunction(fes)
    BVP(bf=a_incompressible, lf=f, gf=u, pre=c)
    return(u)


def getVanMises(stress_components):
    return sqrt(0.5*(stress_components[0]-stress_components[4])*(stress_components[0]-stress_components[4])+
                0.5 * (stress_components[4] - stress_components[8]) * (stress_components[4] - stress_components[8])+
                0.5 * (stress_components[0] - stress_components[8]) * (stress_components[0] - stress_components[8])+
                3*(stress_components[1]*stress_components[1]+stress_components[2]*stress_components[2]*stress_components[3]*stress_components[3]))

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


def extrapolate_E(E, nu, nu_incompressible=0.49):
    return E*(1+nu_incompressible)/(1+nu)

def get_deformation_with_yield_strain(intended_plastic_deformation,yield_strain):
    return IfPos(getVanMises(intended_plastic_deformation)-yield_strain,(getVanMises(intended_plastic_deformation)-yield_strain) / (
                getVanMises(intended_plastic_deformation)) * intended_plastic_deformation,CoefficientFunction((0, 0, 0, 0, 0, 0, 0, 0, 0)))

def getEvacuationTimeConstant(mesh,K,E,nu,nu_i):
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






def get_simulation(E=2e6, nu=0,nu_i=0.4975,sd=0.1,K=1e5,dt_max=0.01,eta_yield=1e6,n_steps=30):
    mesh = loadMesh()
    evacuationTimeConstant=getEvacuationTimeConstant(mesh, K, E, nu, nu_i)
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
    stress_deformatory_old = CoefficientFunction((0, 0, 0, 0, 0, 0, 0, 0, 0))
    for ind in range(n_steps):
        if ind>0:
            u_old.Set(u)
        u = run_simulation_step(fes, mesh, extrapolate_E(E, nu), nu_i, volumeChange, plastic_deformation, u_old,
                                    stress_deformatory_old)
        for gg in range(5):
            u = run_simulation_step(fes, mesh, extrapolate_E(E, nu), nu_i, volumeChange, plastic_deformation, u,
                                stress_deformatory_old)
            Draw(u,mesh,"u")
        u = run_simulation_step(fes, mesh, extrapolate_E(E, nu), nu_i, volumeChange, plastic_deformation, u,
                                stress_deformatory_old)
        u=getSmoothenedNormalizedAxisymmetric(mesh, u, sd)
        u_scaffold=run_simulation_step(fes, mesh, E, nu, GridFunction(fesx), plastic_deformation, u_old,
                                stress_deformatory_old)
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
        kappa=3*E/(1-2*nu)
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
        volumeChange = getSmoothenedNormalizedAxisymmetric(mesh,volumeChangeNew,sd)
        Draw(volumeChange, mesh, "volume_change")
        global plunger_displacement
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






# Penalty factors
penalty_general=1e10
plunger_z_penalty_factor = 1e10 # This one can be seen as an external force. If we have 5N on 20^2mm^^*pi, that's 12500Pa
barrel_penalty_z_factor=1e10
plunger_displacement=0 # Starting value, will be updated
plunger_displacement_rate_m_sec = -1e-4 # to be used with imposed displacement. The area ratio is 400, so this can't be too big
max_plunger_displacement=-5e-4
r_penalty_factor = penalty_general
r_penalty_factor_barrel=penalty_general
xy_plunger_penalty_factor=penalty_general
penalty_transition=penalty_general
penalty_tubing=penalty_general
penalty_xy_torus=penalty_general/1000


# geometric positions
# this is the torus having been cut out to permit the transition between the barrel and the needle without singulariy
torus_position_z=-0.001
radius_torus=0.001
radius_ring=0.0015

radius_main_barrel=0.002
z_lower_main_barrel=0
z_upper_main_barrel=0.01


lower_end_outlet=-0.1

# For the step size
safety_factor=1

volume_change_max=0.8


dt_max=1


n_steps=300


E=2e3

nu=0.1
nu_i=0.4999

sd=0.2e-3
sd_laplace=2e-3

K=7.1e7



yield_strain=0.2

# In practice, this is eta_water/Ks, Ks=r_pore^2/64 so K=64*eta_water/r_pore^2
# so for water, and a pore radius of 3e-5m

# Eta yield: At some 100Pa stress, the material yields and deformation rates are on the order 10/s
#eta_yield=1e3
eta_yield=1e3

p_offset=300 # Negative pressure in material, must be offset by fluid to leave as drop




ret=get_simulation(E=E, nu=nu,nu_i=nu_i, sd=sd,K=K,eta_yield=eta_yield,dt_max=dt_max,n_steps=n_steps)


mesh=ret["mesh"]
u=ret["u"][0]
#us=getSmoothenedNormalized(mesh,u,sd)
#fesp = H1(mesh, order=1, dim=9)
#stress = GridFunction(fesp)
#stress.Set(sigma(epsilon(us), E, nu))

chdir("/Users/unige/Documents/publications/injectability/solution")
ret["evacuationTimeConstant"].Save("evacuationTimeConstant.sol")
for ind in range(len(ret["u"])):
    ret["u"][ind].Save("u_"+str(ind)+".sol")
    fesp = H1(mesh, order=1, dim=9)
    stress = GridFunction(fesp)
    stress.Set(ret["stress"][ind])
    stress.Save("stress_" + str(ind) + ".sol")
    strain = GridFunction(fesp)
    strain.Set(ret["strain"][ind])
    strain.Save("strain_" + str(ind) + ".sol")
    ret["V"][ind].Save("V_" + str(ind) + ".sol")
    ret["plasticDeformation"][ind].Save("plasticDeformation_" + str(ind) + ".sol")
    ret["delta_plastic"][ind].Save("delta_plastic" + str(ind) + ".sol")

file = open("t.txt", "w")
file.write(str(ret["t"]))
file.close()

file = open("plunger_position.txt", "w")
file.write(str(ret["plunger_position"]))
file.close()



file = open("parameters.txt", "w")
file.write("E="+str(E)+"\n"+
           "nu="+str(nu)+"\n"+
           "nu_i"+str(nu_i)+"\n"+
           "sd="+str(sd)+"\n"+
           "K="+str(K)+"\n"+
           "eta_yield="+str(eta_yield)+"\n"+
           "yield_strain"+str(yield_strain)+"\n"+
           "dt_max="+str(dt_max)+"\n"+
           "penalty_general="+str(penalty_general)+"\n"+
           "plunger_displacement_initial=0\n"+
           "plunger_displacement_rate_m_sec="+str(plunger_displacement_rate_m_sec)+"\n"+
           "max_plunger_displacement="+str(max_plunger_displacement)+"\n"
           )
file.close()



