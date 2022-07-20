# The idea here is to optimize the density of the sphere such that given:
# - the touching area (defined as flat)
# - the mechanical properties and terrestrial acceleration
# we can find a solution such that
# - the stress at the edge of the touching area has an average z-component of zero
#
# The last condition is the natural condition for a flat continuation of the surface, that is, the object at equilibrium
# on its touching area


from ngsolve import *
import math

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





