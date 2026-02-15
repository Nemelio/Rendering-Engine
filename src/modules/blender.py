
import bpy
import mathutils as mu
from helpers import vec

class ObjectShader:
    """
    Class representing an object material. 
    It is composed of (Cd, Cs, sh, ka, reflec, transparence, idr), respectively: 
    - The diffuse color
    - The specular reflection color
    - The shininess
    - The ambiant shading factor
    - Reflection level
    - Transparancy level
    - Index of refraction
    """
    def __init__(self,  Cd, Cs, sh, ka, reflec, transparence, idr):
        self.Cd = Cd
        self.Cs = Cs
        self.sh = sh
        self.ka = ka
        self.reflec = reflec
        self.transparence = transparence
        self.idr = idr
        
    
    def get(self):
        return self.Cd, self.Cs, self.sh, self.ka, self.reflec, self.transparence, self.idr

class Shaders:
    """Collections of shader settings associated with objects of the scene.
    """
    def __init__(self):
        self.data = {}
        
    def set_shader(self, name : str,
                   Cd: mu.Vector= vec((50,50,50)),
                   Cs: mu.Vector= vec((1,1,1)),
                   sh: float=0.0,
                   ka: float=0.2,
                   reflec : float=0.0,
                   transparence : float=0.0,
                   idr : float=1.0
                   ) -> None:
        self.data[name] = ObjectShader(Cd, Cs,sh,ka, reflec, transparence, idr)
    
    def get_shader(self, name: str) -> tuple:
        if self.data.get(name):
            return self.data.get(name).get()
        return ObjectShader(
            vec((50,50,50)),
            vec((1,1,1)),
            0.0, 0.2, 0.0, 0.0, 1.0).get()

def enclosing_box(obj : bpy.types.Object) -> tuple: 
    """  
    Calculates the minimal and maximal value of all vertices' coordinates for each axis.

    Args:
        obj (bpy.types.Object): An object existing in the blender scene

    Returns:
        (`min`, `max`) (tuple) : Two vectors each representing the positions of the minimum `(xmin, ymin, zmin)` and maximum `(xmax, ymax, zmax)` point in the object.
    """
    M = obj.matrix_world
    v0 = obj.data.vertices[0].co
    xmin, ymin, zmin = M @ v0
    xmax, ymax, zmax = M @ v0
    
    for vertex in obj.data.vertices:
        x, y, z = M @ vertex.co
        
        if x > xmax : xmax = x
        if y > ymax : ymax = y
        if z > zmax : zmax = z
        
        if x < xmin : xmin = x
        if y < ymin : ymin = y
        if z < zmin : zmin = z
    
    
    pmin = mu.Vector((xmin, ymin, zmin))
    pmax = mu.Vector((xmax, ymax, zmax))

    return (pmin, pmax)

def vector_to_world(obj: bpy.types.Object ,vector: mu.Vector) -> mu.Vector:
    return (obj.matrix_world.to_3x3() @ vector).normalized()

def coord_to_world(obj: bpy.types.Object ,coord: mu.Vector) -> mu.Vector:
    return (obj.matrix_world @ coord)

def objects_by_type(type : str) -> list[bpy.types.Object]:
    """Get a list of all the objects of a given type from the scene.

    Args:
        type (str): The type of objects to retrieve

    Returns:
        list[bpy.types.Object]: List of objects matching the specified type.
    """
    objects = []
    for obj in bpy.data.objects:
        if obj.type == type:
            objects.append(obj)
    return objects 

def get_vertex_by_index(obj : bpy.types.Object, index: int) -> mu.Vector:
    """
    Calculates the coordinates of a vertex using its index.
    
    Args:
        obj (bpy.types.Object): The base object of the vertex.
        n (int): The index of the vertex inside the object.

    Returns:
        Vector: The coordinates of the vertex transformed by the world matrix
    """
    return obj.matrix_world @ obj.data.vertices[index].co

def get_hitboxes(meshes : list, hide_viewport = True) -> list: 
    """
    Generates a list of hitboxes for each meshes specified in the list from argument.
    
    Args:
        meshes (list): List of meshes of the scene.
        hide_viewport (bool, optional): If the hitboxes need to be hidden from the viewport. Defaults to True.

    Returns:
        list: A list of blender objects representing the hitboxes of each meshes.
    """
    hitboxes = []
    for mesh in meshes: 
        pmin, pmax = enclosing_box(mesh)
        hitbox = create_enclosing_box(pmin, pmax)
        hitbox.name = "hitbox"
        hitbox.hide_viewport = hide_viewport

        hitboxes.append(hitbox)
        
    return hitboxes

def create_enclosing_box(pmin: mu.Vector, pmax: mu.Vector) -> bpy.types.Object:
    """
    Creates a triangulated cube with `pmin` and `pmax` as it bounding box.
    
    Args:
        pmin (mu.Vector): The coordinates of the lower point of the cube
        pmax (mu.Vector): The coordinates of the upper point of the cube

    Returns:
        bpy.types.Object: The triangulated cube generated.
    """
    
    # Calculates the average position of the two points
    location = (pmin + pmax) / 2
    
    
    # Calculates size of the cube using the delta of both positions
    delta = pmax - pmin
    dx, dy, dz = delta
    size = max(dx,dy,dz)
    
    bpy.ops.mesh.primitive_cube_add(size=size, enter_editmode=True, location=location, scale=(1,1,1))
    obj = bpy.context.active_object
    bpy.ops.mesh.quads_convert_to_tris()
    bpy.ops.object.editmode_toggle()
    
    return obj

def delete_hitboxes(hitboxes: list) -> None:
    """Delete all the hitboxes specified in argument.

    Args:
        hitboxes (list): List of the hitboxes of the scene.
    """
    for hitbox in list(hitboxes):
        if hitbox is not None and hitbox.name in bpy.data.objects:
            bpy.data.objects.remove(hitbox, do_unlink=True)
    hitboxes.clear()
