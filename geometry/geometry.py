# The idea here is to optimize the density of the sphere such that given:
# - the touching area (defined as flat)
# - the mechanical properties and terrestrial acceleration
# we can find a solution such that
# - the stress at the edge of the touching area has an average z-component of zero
#
# The last condition is the natural condition for a flat continuation of the surface, that is, the object at equilibrium
# on its touching area


from ngsolve import *
import netgen.meshing as nm


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
    n=normal_vector_torus_segment(u,torus_position_z=torus_position_z,radius_ring=radius_ring)
    current_position = CoefficientFunction((x + u[0], y + u[1], z + u[2]))
    # We have the surface normal. To go to the torus center, we need to get the z position right
    d_to_torus_center = IfPos(n[2],(torus_position_z-current_position[2])/n[2],
                              IfPos(n[0]*current_position[0]+n[1]*current_position[1],
                                    radius_ring-(n[0]*current_position[0]+n[1]*current_position[1]),
                                    n[0] * current_position[0] + n[1] * current_position[1]-radius_ring))
    return current_position + IfPos(d_to_torus_center,(d_to_torus_center-radius_torus)*n,(radius_torus+d_to_torus_center)*n)

def get_targetposition_on_cylinder_surface(u,radius=0.0005):
    n = normal_vector_cylinder(u)
    return CoefficientFunction((n[0]*radius,n[1]*radius,u[2]+z))



