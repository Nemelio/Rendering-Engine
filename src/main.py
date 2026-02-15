import bpy, mathutils as mu
from modules.blender import Shaders
from modules.core import generate_camera_ray_grid, trace_rays, render
from modules.blender import objects_by_type, get_hitboxes, delete_hitboxes
import time, sys



if __name__ == "__main__":
    
    
    print("[Start Program]")
    start = time.time()   
    
    # Essentials objects
    cam = bpy.data.objects["Camera"]
    light = bpy.data.objects["Light"]
    
    # Shaders definitions
    shaders = Shaders()
    shaders.set_shader("Sphere", mu.Vector((255,0,0)),Cs=mu.Vector((255,255,255)), ka=0.1, sh=64, reflec=0.5, transparence=0.4)
    shaders.set_shader("Sphere.001", mu.Vector((100,100,0)),Cs=mu.Vector((255,255,255)), ka=0.2, sh=64, reflec=0.7,)
    shaders.set_shader("Light", mu.Vector((255,255,255)))
    shaders.set_shader("Plane", mu.Vector((0,0,255)), ka=0.1)
    shaders.set_shader("Plane.001", mu.Vector((255,0,255)))
    shaders.set_shader("Plane.002", mu.Vector((0,255,255)))
    shaders.set_shader("Plane.003", mu.Vector((255,100,255)))

    # Rays precomputing
    rays = generate_camera_ray_grid(cam, dpx=0.01)
   
    meshs = objects_by_type("MESH")
    hitboxes = get_hitboxes(meshs)
    

    
    # Ray tracing + rendering
    try:
        root = sys.path[-1]
        path = f"{root}/ppm"
        file_name = "/save.ppm"
     
        nb_max_reflec = 3
        nb_max_refr = 3

        pixels = trace_rays(cam.location, rays, meshs, hitboxes, shaders, light, cam, nb_max_reflec, nb_max_refr)
        im_largeur = bpy.context.scene.render.resolution_x
        im_hauteur = bpy.context.scene.render.resolution_y
        
        render(path=f"{path}{file_name}",
                pixels=pixels,
                height=im_hauteur,
                width=im_largeur)
    finally:
        # Delete all hitboxes
        delete_hitboxes(hitboxes)
        pass 
    
    
    end = time.time() - start
    print("[End Program]")
    print(f"Executing Time : {end}")