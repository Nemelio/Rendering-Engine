
import bpy, math
import mathutils as mu
from blender import get_vertex_by_index, Shaders, vector_to_world
from helpers import lencol, vec, direction_vector, clamp_color
from maths import calc_refl_vector, calc_refr_vector, build_local_frame, phong, lambert
from ressources import PPM
from ressources.Barycentre import coords3

def hit_object(nbobj : int, p_int : mu.Vector, objs : list, hitboxes: list,  light : bpy.types.Object, eps=1e-4) -> tuple[bool, bpy.types.Object]: 
    """
    Determines whether an object lies between the intersection point and the light source 
    by casting a ray from the point toward the light. 

    Args:
        nbobj (int): Index of the original object
        p_int (Vector): Intersection point
        objs (list): List of scene objects
        hitboxes (list): Bounding boxes corresponding to each object
        light (Vector): Light object
        eps (float, optional): Epsilon offset applied to the ray origin. Defaults to 1e-4.

    Returns:
        tuple : A tuple (is_hit, object). If an object blocks the light, is_hit is True and object is that object; otherwise is_hit is False and object is None.
    """
    # Ray from intersection toward light
    ray = (light.location - p_int)
    
    distance_ray_light = ray.length
    u_ray = ray.normalized()
    p_origin = p_int + u_ray * eps
    
    for idx, obj in enumerate(objs):
        # Ignore the source object to avoid self intersection
        if idx != nbobj:
            hit,_ , _= rayon_intersect_opt(p_origin, u_ray, [hitboxes[idx]], [obj], distance_ray_light)     
            if hit is not None:
                return True, obj
            
    return False, None

def rayon_intersect1(p_source : mu.Vector, u_ray : mu.Vector, obj : bpy.types.Object) -> tuple:
    """
    Computes the intersection point between a ray and a given object. 

    Args:
        p_source (mu.Vector): Ray origin coordinates
        u_ray (mu.Vector): Ray direction vector
        obj (bpy.types.Object): Object intersected

    Returns:
        tuple: Tuple (p_int, nf) where `p_int` is the location of the nearest intersected point and `nf` is the index of the intersected face of the object.
    """
   
    
    u = u_ray.normalized()
    p_int = None 
    nf_best = None
    best_d2 = None
    
    
    for polygons in obj.data.polygons:
        buffer = []
    
        
        for idx_vertex in polygons.vertices:
            # Retrivies vertices coordinates
            buffer.append(get_vertex_by_index(obj, idx_vertex))
            
            # Ignore if the buffer has less than 3 vertex 
            if len(buffer) < 3: 
                continue
            
            # Get the 3 vertices and calculates the intersection point
            A,B,C = buffer
            hit = mu.geometry.intersect_ray_tri(A,B,C,u,p_source)
            
            # Computes distances and retrieves the nearest intersection point
            if hit is not None:
                d2 = (hit - p_source).length_squared
                if best_d2 is None or d2 < best_d2:
                    best_d2 = d2
                    p_int = hit
                    nf_best = polygons.index 
                    
            buffer.pop(0)
            
            

            
    return (p_int, nf_best)

def rayon_intersect_opt(p_source: mu.Vector, u_ray: mu.Vector, hitboxes : list , meshs : list, l_max : float = 0.0 ) -> tuple[mu.Vector, int, int]:
    """   
    Computes ray intersections with a set of meshes.
    Returns the nearest intersection point.
    Uses the bounding box of each mesh to reduce the number of tested faces.

    Args:
        p_source (mu.Vector): Ray origin coordinates
        u_ray (mu.Vector): Ray direction vector
        objs (list[bpy.types.Object]): List of scene meshes
        l_max (float, optional): Maximum distance between the ray origin and the object. Defaults to 0.

    Returns:
        tuple: A tuple (`p_int`: Vector, `nobj`: int, `nf`: int). `p_int` is the nearest intersection point, `nobj` is the index of the intersected object in the list and `nf` is the index of the intersected face.
    """
    
    obj_nf = None
    p_int = None
    nobj = None
    best_d2 = None

    for idx, box in enumerate(hitboxes):
        hit, _ = rayon_intersect1(p_source, u_ray, box)
        
        if hit is None:
            continue
        
        objet = meshs[idx]
        # Retrives the intersection point and the index of the intersected face
        tmp_int, tmp_nf = rayon_intersect1(p_source, u_ray, objet)
        
        if tmp_int is None:
            continue
        
        ## Computes the distance and ignores if the distance is above `l_max`
        # If two objects are intersected, keep the nearest point
        d2 = ( tmp_int - p_source).length_squared
        d = ( tmp_int - p_source).length
        if best_d2 is None or d2 < best_d2 :
            

            if l_max > 0.0 and d > l_max:
                continue
            
            obj_nf = tmp_nf
            p_int = tmp_int
            nobj = idx
            best_d2 = d2
        
    return (p_int, nobj, obj_nf)
        
def render(path : str, pixels : list[list], height : int, width: int) -> None:
    """    
    Generates a PPM image from a matrix of pixels.

    Args:
        path (str): Output file path
        pixels (list): Matrix HxW populated with RGB vectors
        height (int): Height of the image
        width (int): Width of the image
    """
    img = PPM.creer(path,width, height)
    
    pixels_height = len(pixels)
    pixels_width = lencol(pixels)
    inversed_height = 1.0 / height
    inversed_width = 1.0 / width
    for lgn in range(height):
        for col in range(width):
            px = int(lgn * inversed_height * pixels_height)
            py = int(col * inversed_width * pixels_width)
        
            color = pixels[px][py]
            PPM.writePixel(img, color)
            

    PPM.fini(img)
    
def calc_pixel_color(
    p_source : mu.Vector,
    u_ray : mu.Vector,
    meshes : list,
    hitboxes: list,
    shaders : Shaders,
    light : bpy.types.Object,
    camera : bpy.types.Object,
    nb_max_reflec: int = 0,
    nb_max_refr: int = 0
) -> mu.Vector:
    """    
    Recursively computes the color of a pixel based on shaders and scene meshes.
    
    If reflectivity is non-zero, it calls `calc_color` and 
    then itself to compute reflection.
    
    If transparency is non-zero, it calls `calc_color` and 
    then itself to compute refraction. 


    Args:
        p_source (mu.Vector): Ray origin
        u_ray (mu.Vector): Ray direction vector
        meshes (list): List of scene meshes
        hitboxes (list): Bounding boxes corresponding of each mesh
        shaders (Shaders): Shader parameters 
        light (bpy.types.Object): Light object 
        camera (bpy.types.Object): Camera object
        nb_max_reflec (int, optional): Maximum reflection recursion depth. Defaults to 0.
        nb_max_refr (int, optional): Maximum refraction recursion depth. Defaults to 0.

    Returns:
        mu.Vector: RGB vector
    """
    
    p_int,nbobj,nf = rayon_intersect_opt(p_source, u_ray, hitboxes, meshes)
    
    
    if p_int is None:
        return vec((1,1,1))
   
    c_local = calc_color(p_source, u_ray, meshes, hitboxes, shaders, light, camera)
    mesh = meshes[nbobj]
    _,_,_,_,reflectivite,transparence,idr = shaders.get_shader(mesh.name)
    
    # If reflectivity and transparency of the object is zero, return the color
    if reflectivite <= 0.0 and transparence <= 0.0:
        return c_local

    
    n = vector_to_world(mesh, mesh.data.polygons[nf].normal)
    
    # Computes the normal orientation relative to the object
    sens = u_ray.dot(n) > 0 
    if sens: n = -n
    
    eps = 1e-3 # Avoid self intersection
    c_refl  = vec((0,0,0))
    c_trans = vec((0,0,0))
    
    # Reflection
    if reflectivite > 0.0 and nb_max_reflec > 0.0:
        urefl = calc_refl_vector(u_ray,n)
        p_refl = p_int + eps * urefl
        c_refl = calc_pixel_color(p_refl, urefl, meshes, hitboxes, shaders, light, camera, nb_max_reflec-1, nb_max_refr)
        
    # Refraction / Transparency
    if transparence > 0.0 and nb_max_refr > 0.0:
        n1, n2 = (idr,1.0) if sens else (1.0, idr)
        u_refr = calc_refr_vector(u_ray, n, n1, n2)
        if u_refr is not None:
            p_refr = p_int + eps * u_refr
            c_trans = calc_pixel_color(p_refr, u_refr, meshes, hitboxes, shaders, light, camera, nb_max_reflec, nb_max_refr-1)
       
    return reflectivite * c_refl + transparence * c_trans + ( 1 - reflectivite - transparence) * c_local

def generate_camera_ray_grid(cam: bpy.types.Object, dpx: float = 0.01) -> list:
    """
    Creates a set of vectors representing points in the camera frame. 
    Each point is then used to define a ray origin and compute a ray direction
    vector.

    Args:
        cam (bpy.types.Object): Camera object
        dpx (float, optional): Spacing between sample points. Defaults to 0.01.

    Returns:
        list: List of vectors.
    """

    alpha = cam.data.angle_x / 2
    Dy = 1.39 
    Dx = Dy * math.tan(alpha) # L'opposé = tan(theta) * adjacent
    Dz = 9/16*Dx # On utilise un rapport 16/9 pour déterminer la résolution Dz de l'écran

    # Calcul les vecteurs en fonction de la position de la caméra
    rays = []
    sx,sy,sz = cam.location
    x_max, z_max = sx+Dx, sz+Dz
    z = sz-Dz
    while z < z_max:
        x = sx-Dx
        line = []
        while x < x_max:
            line.append(mu.Vector((x,(sy+Dy),z)))
            x += dpx
        rays.append(line)
        z += dpx

    return rays

def calc_color(
    p_source : mu.Vector,
    u_ray : mu.Vector,
    meshes : list,
    hitboxes: list,
    shaders : Shaders,
    light : bpy.types.Object,
    camera : bpy.types.Object,
    ) -> mu.Vector:
    """   
    Casts a ray toward the scene, finds the intersection point 
    with an object, casts a shadow ray and if the point is 
    illuminated, computes and returns its color.
    
    Args:
        p_source (mu.Vector): Ray origin
        u_ray (mu.Vector): Ray direction vector
        meshes (list): List of scene meshes
        hitboxes (list): Bounding boxes corresponding to each mesh
        shaders (Shaders): Shaders parameters
        light (bpy.types.Object): Light object
        camera (bpy.types.Object): Camera object

    Returns:
        mu.Vector: RGB vector
    """
    
    color = vec((0,0,0))
    Ci, _,_,_,_,_,_ = shaders.get_shader(light.name)
    intensity = light.data.energy / 1000.0 

    p_int,nbobj,nf = rayon_intersect_opt(p_source, u_ray, hitboxes, meshes)
    if p_int is None or nbobj is None or nf is None:
        return color
    
    
    
    cam_ray = direction_vector(p_int, camera.location)
    mesh = meshes[nbobj]
   
    MW = mesh.matrix_world
    
    Cd, Cs, sh, ka, _, _, _ = shaders.get_shader(mesh.name)
    poly = mesh.data.polygons[nf]
    

    if len(poly.vertices) == 3:
        
        iv0, iv1, iv2 = poly.vertices
        v0 = MW @ mesh.data.vertices[iv0].co
        v1 = MW @ mesh.data.vertices[iv1].co
        v2 = MW @ mesh.data.vertices[iv2].co

    
        face_normal_world = vector_to_world(mesh, poly.normal) 

        # construction du repère local
        M_local = build_local_frame(face_normal_world)
        
        # on travaille dans un repère centré sur v0 pour éviter les grosses valeurs
        origin = v0

        P_local  = M_local @ (p_int - origin)
        v0_local = M_local @ (v0 - origin)
        v1_local = M_local @ (v1 - origin)
        v2_local = M_local @ (v2 - origin)
        


        # Local vectors must not be normalized
        a,b,c = coords3(P_local, v0_local, v1_local, v2_local)
        
        
        nv0 = vector_to_world(mesh, mesh.data.vertices[iv0].normal)
        nv1 = vector_to_world(mesh, mesh.data.vertices[iv1].normal)
        nv2 = vector_to_world(mesh, mesh.data.vertices[iv2].normal)
        
        
        nINT = a * nv0 + b * nv1 + c * nv2
        nINT = nINT.normalized() 
        
        
    
        u_light = direction_vector(p_int, light.location)
        diff_amb = lambert(u_light, nINT, Ci, Cd, ka)
        spec = phong(cam_ray, u_light, nINT, Ci, Cs, sh)
   
        diff_amb_spec = diff_amb + spec 
        
        hit,hit_obj = hit_object(nbobj, p_int,hitboxes, meshes, light)
        
        if hit:
            _, _, _, _, _, hit_trans, _ = shaders.get_shader(hit_obj.name)
            obscurite = max(0.0, min(1.0, (ka + hit_trans)))
            color = diff_amb_spec * obscurite * 255 
            
        else:
            color =  diff_amb_spec * intensity * 255
    

    return color
    
def trace_rays(
    p_source : mu.Vector,
    rays : mu.Vector,
    meshs : list,
    hitboxes: list,
    shaders : Shaders,
    light : bpy.types.Object,
    camera : bpy.types.Object,
    nb_max_reflec: int = 0,
    nb_max_refr: int = 0
): 
    """Algorithme de lancé de rayon. \n
    
    **Objectif** \: A partir d'une source, crée des vecteurs "d'observation"
    et calcul la couleur de chaque pixel en fonction de différents paramètres (Shaders). 
    
    **Particularités** \: La fonction utilise une fonction récursive pour calculer 
    les couleurs et des hitboxes pour optimiser les calculs. 
    
    :param p_source: Source de l'observateur (origin du rayon)
    :type p_source: mu.Vector
    
    :param rays: Liste comportant tous les points utilisés pour calculer les vecteurs
    d'observation.
    :type rays: mu.Vector
    
    :param meshs: Ensemble des objets concernés par l'observation (générallement des polygones)
    :type meshs: list
    
    :param hitboxes: Ensemble des hitboxes associés aux objets dans `meshs`
    :type hitboxes: list
    
    :param shaders: Objet/class comportant toutes les textures de tous les objets
    :type shaders: Shaders
    
    :param light: Objet blender représentant la source de lumière
    :type light: bpy.types.Object
    
    :param camera: Objet blender représentant la caméra (observateur)
    :type camera: bpy.types.Object
    
    :param nb_max_reflec: Nombre maximum de reflexion
    :type nb_max_reflec: int
    
    :param nb_max_refr: Nombre maximum de refraction
    :type nb_max_refr: int

    Args:
        p_source (mu.Vector): _description_
        rays (mu.Vector): _description_
        meshs (list): _description_
        hitboxes (list): _description_
        shaders (Shaders): _description_
        light (bpy.types.Object): _description_
        camera (bpy.types.Object): _description_
        nb_max_reflec (int, optional): _description_. Defaults to 0.
        nb_max_refr (int, optional): _description_. Defaults to 0.

    Returns:
        _type_: _description_
    """
    
    # Lance un rayon en fonction des vecteurs de la camera
    # Récupère les points d'intersection avec un objet 
    h = len(rays)
    pixels = [None] * h
    nbr_of_pixels = h * lencol(rays)

    pixel_count = 0
    for idx, line in enumerate(rays):
        print(f"[{(pixel_count/nbr_of_pixels)*100}/100%]", end="\r")
        
        
        pixels_line = []
        for ray in line:
            pixel_count+=1
            # Regarde si un objet est intersecté 
            u_ray = direction_vector(p_source, ray)
            
            # Fonction récursive refraction
            color = calc_pixel_color(p_source, u_ray, meshs, hitboxes, shaders, light, camera, nb_max_reflec, nb_max_refr)
            
            # Fonction récursive reflexion
            ## color = calc_color_refl(p_source, u_ray, meshs, hitboxes, shaders, light, camera, nb_max_reflec)
           
               
            pixels_line.append(clamp_color(color))
         
        # Remplissage inverse
        # Sinon l'image est a l'envers
        pixels[h - 1 - idx] = pixels_line
        
    
    return pixels
