import math, mathutils as mu

def my_reflect(u_inc : mu.Vector, normal: mu.Vector) -> mu.Vector:
    """Calculates the reflection of an incident vector.

    Args:
        u_inc (Vector): The incident direction vector
        normal (Vector): The normal vector of the intersected surface.

    Returns:
        Vector: The reflected direction vector
    """
    normal.normalize()
    u_inc.normalize()
    u_inc_norm = (normal.dot(u_inc))*normal
    u_inc_tang = u_inc - u_inc_norm
    return u_inc_norm - u_inc_tang

def phong(uo:mu.Vector,
          ui: mu.Vector,
          n: mu.Vector,
          Ci: mu.Vector,
          Cs: mu.Vector,
          sh : float) -> mu.Vector:
    """
    The specular component computed using the Phong model and characterized by a specular reflection
    color and shininess.

    Args:
        uo (Vector): The view direction
        ui (Vector): The incident direction vector of the light
        n (Vector): The normal of the intersected surface
        Ci (Vector): The color of the incident light
        Cs (Vector): The specular color
        sh (float): The shininess factor

    Returns:
        Vector: RGB color of the pixel
    """
    
    Ci = Ci / 255.0
    Cs = Cs / 255.0
 
    ur = my_reflect(ui, n)
    scalaire = max(0.0 ,uo.dot(ur))
    return Ci * Cs * (scalaire ** sh)

def lambert(ui : mu.Vector,
            normal: mu.Vector,
            Ci : mu.Vector,
            Cd : mu.Vector = mu.Vector((1,1,1)),
            darkness : float = 0.1
            ) -> mu.Vector:
    """    
    The diffusion component is computed using the Lambert model. It calculates the 
    diffuse intensity using the angle between light direction vector and the 
    surface normal.

    Args:
        ui (mu.Vector): The incident direction vector of the light
        normal (mu.Vector): The normal direction
        Ci (mu.Vector): The color of the incident light
        Cd (mu.Vector, optional): The diffuse color. Defaults to Vector((1,1,1)).
        darkness (float, optional): The ambient factor. Defaults to 0.1.

    Returns:
        Vector: RGB color of the pixel
    """
    # Normalize all vectors to get coherent calculation
    Ci = Ci / 255.0
    Cd = Cd / 255.0
    

    dot = max(0.0, ui.dot(normal))
    

    Ca = (Cd * darkness)
    Cdiff = (Ci * Cd * dot) #dot = scalaire
    
    # Apply the ambient component to the diffuse components
    return Ca + Cdiff

def calc_refr_vector(d_vec : mu.Vector, n_face : mu.Vector, n1 : float, n2: float):
    """
    Computes the refracted direction vector.

    Args:
        d_vec (Vector): Incident direction vector
        n_face (Vector): Surface normal vector
        n1 (float): Refractive index of the first environment
        n2 (float): Refractive index of the second environment

    Returns:
        Vector: Refracted direction vector vector
    """
    cos_i = -d_vec.dot(n_face)
    eta = n1/n2
    
    k = 1.0 - eta*eta * (1.0 - cos_i*cos_i)
    
    if k <= 0.0:
        return None
    
    return (eta * d_vec + (eta * cos_i - math.sqrt(k)) * n_face).normalized()
 
def calc_refl_vector(u_ray : mu.Vector, n: mu.Vector) -> mu.Vector:
    """Computes the reflection of an incident direction vector
    about a surface normal.

    Args:
        u_ray (Vector): Normalized incident direction vector.
        n (Vector): Normalized surface normal vector.

    Returns:
        Vector: Normalized reflected direction vector.
    """
    udotn = u_ray.dot(n)
    return (u_ray - 2.0 * udotn * n).normalized()

def build_local_frame(v_normal: mu.Vector) -> mu.Matrix:
    """    
    Builds a local orthonormal frame from a given direction.
    In other words, it converts a single direction into a complete local coordinate system.
    
    A normal vector alone does not define a complete orientation,
    which can hinder computations.
    

    Args:
        v_normal (mu.Vector): Surface normal direction vector

    Returns:
        Matrix: 3x3 matrix whose rows are the tangent, bitangent, and normal vectors forming an orthonormal basis.
    """
 
    up = mu.Vector((0.0, 0.0, 1.0))
    if abs(v_normal.dot(up)) > 0.99:
        up = mu.Vector((0.0, 1.0, 0.0))
    t = v_normal.cross(up).normalized()      # tangent
    b = v_normal.cross(t).normalized()       # bitangent

    return mu.Matrix((t, b, v_normal))
