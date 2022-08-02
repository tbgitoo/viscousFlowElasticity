# The idea here is to optimize the density of the sphere such that given:
# - the touching area (defined as flat)
# - the mechanical properties and terrestrial acceleration
# we can find a solution such that
# - the stress at the edge of the touching area has an average z-component of zero
#
# The last condition is the natural condition for a flat continuation of the surface, that is, the object at equilibrium
# on its touching area


from ngsolve import *
import ngsolve

# outward pointing normal for cylinder
def normal_vector_cylinder(u):
    x_f = CoefficientFunction(x + u[0])
    y_f = CoefficientFunction(y + u[1])
    r_xy = sqrt(x_f * x_f + y_f * y_f)
    return IfPos(r_xy, CoefficientFunction((x_f / r_xy , y_f / r_xy , 0)),
        CoefficientFunction((1, 0, 0)))


def normal_vector_torus(u,torus_position_z=-0.001,radius_ring=0.0015):
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
    return basic_normvector



def get_targetposition_on_torus_surface(u,torus_position_z=-0.001,radius_ring=0.0015,
                                          radius_torus=0.001):
    n=normal_vector_torus(u,torus_position_z=torus_position_z,radius_ring=radius_ring)
    current_position = CoefficientFunction((ngsolve.x + u[0], ngsolve.y + u[1], ngsolve.z + u[2]))
    # We have the surface normal. To go to the torus center, we need to get the z position right
    d_to_torus_center = (current_position[2]-torus_position_z)/n[2]
    return current_position +(radius_torus-d_to_torus_center)*n

def get_targetposition_on_cylinder_surface(u,radius=0.0005):
    n = normal_vector_cylinder(u)
    return CoefficientFunction((n[0]*radius,n[1]*radius,u[2]+z))

def tangential_vector_horizontal_plane(position):
    r=sqrt(position[0]*position[0]+position[1]*position[1])
    return IfPos(r,
          CoefficientFunction((-position[1]/r,position[0]/r,0)),
          CoefficientFunction((0,0,0)))

def vector_product(a,b):
    return CoefficientFunction((a[1]*b[2]-a[2]*b[1],a[0]*b[2]-a[2]*b[0],a[0]*b[1]-a[1]*b[0]))

def vector_product_numbers(a,b):
    return [a[1] * b[2] - a[2] * b[1], a[0] * b[2] - a[2] * b[0], a[0] * b[1] - a[1] * b[0]]


def totalMeshVolume(mesh,u):
    fes3 = H1(mesh, order=1, dim=3)
    fes1 = H1(mesh, order=1, dim=1)
    u1 = GridFunction(fes3)
    u1.Set(u)  # to be sure that we have a displacement per node
    x_coords = []
    y_coords = []
    z_coords = []
    V_u_tot=0
    for v in mesh.vertices:
        x_coords.append(v.point[0])
        y_coords.append(v.point[1])
        z_coords.append(v.point[2])
    for tetraeder in mesh.Elements():
        v1 = [x_coords[tetraeder.vertices[1].nr] - x_coords[tetraeder.vertices[0].nr],
              y_coords[tetraeder.vertices[1].nr] - y_coords[tetraeder.vertices[0].nr],
              z_coords[tetraeder.vertices[1].nr] - z_coords[tetraeder.vertices[0].nr]]
        v1_u = [v1[0] + u1.vec.data[tetraeder.vertices[1].nr][0] - u1.vec.data[tetraeder.vertices[0].nr][0],
                v1[1] + u1.vec.data[tetraeder.vertices[1].nr][1] - u1.vec.data[tetraeder.vertices[0].nr][1],
                v1[2] + u1.vec.data[tetraeder.vertices[1].nr][2] - u1.vec.data[tetraeder.vertices[0].nr][2]]
        v2 = [x_coords[tetraeder.vertices[2].nr] - x_coords[tetraeder.vertices[0].nr],
              y_coords[tetraeder.vertices[2].nr] - y_coords[tetraeder.vertices[0].nr],
              z_coords[tetraeder.vertices[2].nr] - z_coords[tetraeder.vertices[0].nr]]
        v2_u = [v2[0] + u1.vec.data[tetraeder.vertices[2].nr][0] - u1.vec.data[tetraeder.vertices[0].nr][0],
                v2[1] + u1.vec.data[tetraeder.vertices[2].nr][1] - u1.vec.data[tetraeder.vertices[0].nr][1],
                v2[2] + u1.vec.data[tetraeder.vertices[2].nr][2] - u1.vec.data[tetraeder.vertices[0].nr][2]]
        v3 = [x_coords[tetraeder.vertices[3].nr] - x_coords[tetraeder.vertices[0].nr],
              y_coords[tetraeder.vertices[3].nr] - y_coords[tetraeder.vertices[0].nr],
              z_coords[tetraeder.vertices[3].nr] - z_coords[tetraeder.vertices[0].nr]]
        v3_u = [v3[0] + u1.vec.data[tetraeder.vertices[3].nr][0] - u1.vec.data[tetraeder.vertices[0].nr][0],
                v3[1] + u1.vec.data[tetraeder.vertices[3].nr][1] - u1.vec.data[tetraeder.vertices[0].nr][1],
                v3[2] + u1.vec.data[tetraeder.vertices[3].nr][2] - u1.vec.data[tetraeder.vertices[0].nr][2]]
        v1v2_u = vector_product_numbers(v1_u, v2_u)
        V_u = abs(v3_u[0] * v1v2_u[0] + v3_u[1] * v1v2_u[1] + v3_u[2] * v1v2_u[2]) / 6
        V_u_tot=V_u_tot+V_u
    return V_u_tot


def deformedMeshVolume(mesh,u):
    fes3 = H1(mesh, order=1, dim=3)
    fes1 = H1(mesh, order=1, dim=1)
    u1=GridFunction(fes3)
    u1.Set(u) # to be sure that we have a displacement per node
    n=GridFunction(fes1) # Number of elements in which a given node is placed, will be used to normalize at the end
    dv=GridFunction(fes1)
    x_coords=[]
    y_coords=[]
    z_coords=[]
    vol=[]
    for v in mesh.vertices:
        x_coords.append(v.point[0])
        y_coords.append(v.point[1])
        z_coords.append(v.point[2])
        vol.append(0)
    for tetraeder in mesh.Elements():
        v1=[x_coords[tetraeder.vertices[1].nr]-x_coords[tetraeder.vertices[0].nr],
            y_coords[tetraeder.vertices[1].nr]-y_coords[tetraeder.vertices[0].nr],
            z_coords[tetraeder.vertices[1].nr]-z_coords[tetraeder.vertices[0].nr]]
        v1_u = [v1[0]+u1.vec.data[tetraeder.vertices[1].nr][0]-u1.vec.data[tetraeder.vertices[0].nr][0],
                v1[1] + u1.vec.data[tetraeder.vertices[1].nr][1] - u1.vec.data[tetraeder.vertices[0].nr][1],
                v1[2] + u1.vec.data[tetraeder.vertices[1].nr][2] - u1.vec.data[tetraeder.vertices[0].nr][2]]
        v2 = [x_coords[tetraeder.vertices[2].nr] - x_coords[tetraeder.vertices[0].nr],
              y_coords[tetraeder.vertices[2].nr] - y_coords[tetraeder.vertices[0].nr],
              z_coords[tetraeder.vertices[2].nr] - z_coords[tetraeder.vertices[0].nr]]
        v2_u = [v2[0] + u1.vec.data[tetraeder.vertices[2].nr][0] - u1.vec.data[tetraeder.vertices[0].nr][0],
                v2[1] + u1.vec.data[tetraeder.vertices[2].nr][1] - u1.vec.data[tetraeder.vertices[0].nr][1],
                v2[2] + u1.vec.data[tetraeder.vertices[2].nr][2] - u1.vec.data[tetraeder.vertices[0].nr][2]]
        v3 = [x_coords[tetraeder.vertices[3].nr] - x_coords[tetraeder.vertices[0].nr],
              y_coords[tetraeder.vertices[3].nr] - y_coords[tetraeder.vertices[0].nr],
              z_coords[tetraeder.vertices[3].nr] - z_coords[tetraeder.vertices[0].nr]]
        v3_u = [v3[0] + u1.vec.data[tetraeder.vertices[3].nr][0] - u1.vec.data[tetraeder.vertices[0].nr][0],
                v3[1] + u1.vec.data[tetraeder.vertices[3].nr][1] - u1.vec.data[tetraeder.vertices[0].nr][1],
                v3[2] + u1.vec.data[tetraeder.vertices[3].nr][2] - u1.vec.data[tetraeder.vertices[0].nr][2]]
        v1v2 = vector_product_numbers(v1,v2)
        v1v2_u = vector_product_numbers(v1_u,v2_u)
        V=abs(v3[0]*v1v2[0]+v3[1]*v1v2[1]+v3[2]*v1v2[2])/6
        V_u = abs(v3_u[0] * v1v2_u[0] + v3_u[1] * v1v2_u[1] + v3_u[2] * v1v2_u[2]) / 6
        vol[tetraeder.vertices[0].nr]=vol[tetraeder.vertices[0].nr]+V_u*V
        vol[tetraeder.vertices[1].nr] = vol[tetraeder.vertices[1].nr] + V_u*V
        vol[tetraeder.vertices[2].nr] = vol[tetraeder.vertices[2].nr] + V_u*V
        vol[tetraeder.vertices[3].nr] = vol[tetraeder.vertices[3].nr] + V_u*V
        n.vec.data[tetraeder.vertices[0].nr]=n.vec.data[tetraeder.vertices[0].nr]+V
        n.vec.data[tetraeder.vertices[1].nr] = n.vec.data[tetraeder.vertices[1].nr] + V
        n.vec.data[tetraeder.vertices[2].nr] = n.vec.data[tetraeder.vertices[2].nr] + V
        n.vec.data[tetraeder.vertices[3].nr] = n.vec.data[tetraeder.vertices[3].nr] + V
    for v in mesh.vertices:
        dv.vec.data[v.nr]=vol[v.nr]/n.vec.data[v.nr]
    return dv



def volumeChangeFromMesh(mesh,u):
    fes3 = H1(mesh, order=1, dim=3)
    fes1 = H1(mesh, order=1, dim=1)
    u1=GridFunction(fes3)
    u1.Set(u) # to be sure that we have a displacement per node
    n=GridFunction(fes1) # Number of elements in which a given node is placed, will be used to normalize at the end
    dv=GridFunction(fes1)
    x_coords=[]
    y_coords=[]
    z_coords=[]
    vol=[]
    for v in mesh.vertices:
        x_coords.append(v.point[0])
        y_coords.append(v.point[1])
        z_coords.append(v.point[2])
        vol.append(0)
    for tetraeder in mesh.Elements():
        v1=[x_coords[tetraeder.vertices[1].nr]-x_coords[tetraeder.vertices[0].nr],
            y_coords[tetraeder.vertices[1].nr]-y_coords[tetraeder.vertices[0].nr],
            z_coords[tetraeder.vertices[1].nr]-z_coords[tetraeder.vertices[0].nr]]
        v1_u = [v1[0]+u1.vec.data[tetraeder.vertices[1].nr][0]-u1.vec.data[tetraeder.vertices[0].nr][0],
                v1[1] + u1.vec.data[tetraeder.vertices[1].nr][1] - u1.vec.data[tetraeder.vertices[0].nr][1],
                v1[2] + u1.vec.data[tetraeder.vertices[1].nr][2] - u1.vec.data[tetraeder.vertices[0].nr][2]]
        v2 = [x_coords[tetraeder.vertices[2].nr] - x_coords[tetraeder.vertices[0].nr],
              y_coords[tetraeder.vertices[2].nr] - y_coords[tetraeder.vertices[0].nr],
              z_coords[tetraeder.vertices[2].nr] - z_coords[tetraeder.vertices[0].nr]]
        v2_u = [v2[0] + u1.vec.data[tetraeder.vertices[2].nr][0] - u1.vec.data[tetraeder.vertices[0].nr][0],
                v2[1] + u1.vec.data[tetraeder.vertices[2].nr][1] - u1.vec.data[tetraeder.vertices[0].nr][1],
                v2[2] + u1.vec.data[tetraeder.vertices[2].nr][2] - u1.vec.data[tetraeder.vertices[0].nr][2]]
        v3 = [x_coords[tetraeder.vertices[3].nr] - x_coords[tetraeder.vertices[0].nr],
              y_coords[tetraeder.vertices[3].nr] - y_coords[tetraeder.vertices[0].nr],
              z_coords[tetraeder.vertices[3].nr] - z_coords[tetraeder.vertices[0].nr]]
        v3_u = [v3[0] + u1.vec.data[tetraeder.vertices[3].nr][0] - u1.vec.data[tetraeder.vertices[0].nr][0],
                v3[1] + u1.vec.data[tetraeder.vertices[3].nr][1] - u1.vec.data[tetraeder.vertices[0].nr][1],
                v3[2] + u1.vec.data[tetraeder.vertices[3].nr][2] - u1.vec.data[tetraeder.vertices[0].nr][2]]
        v1v2 = vector_product_numbers(v1,v2)
        v1v2_u = vector_product_numbers(v1_u,v2_u)
        V=abs(v3[0]*v1v2[0]+v3[1]*v1v2[1]+v3[2]*v1v2[2])/6
        V_u = abs(v3_u[0] * v1v2_u[0] + v3_u[1] * v1v2_u[1] + v3_u[2] * v1v2_u[2]) / 6
        vol[tetraeder.vertices[0].nr]=vol[tetraeder.vertices[0].nr]+V_u-V
        vol[tetraeder.vertices[1].nr] = vol[tetraeder.vertices[1].nr] + V_u - V
        vol[tetraeder.vertices[2].nr] = vol[tetraeder.vertices[2].nr] + V_u - V
        vol[tetraeder.vertices[3].nr] = vol[tetraeder.vertices[3].nr] + V_u - V
        n.vec.data[tetraeder.vertices[0].nr]=n.vec.data[tetraeder.vertices[0].nr]+V
        n.vec.data[tetraeder.vertices[1].nr] = n.vec.data[tetraeder.vertices[1].nr] + V
        n.vec.data[tetraeder.vertices[2].nr] = n.vec.data[tetraeder.vertices[2].nr] + V
        n.vec.data[tetraeder.vertices[3].nr] = n.vec.data[tetraeder.vertices[3].nr] + V
    for v in mesh.vertices:
        dv.vec.data[v.nr]=vol[v.nr]/n.vec.data[v.nr]
    return dv












