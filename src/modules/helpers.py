import mathutils as mu

def minmax_color(n : float) -> float:
    """Clamps the numeric value between min (0) and max (255)

    Args:
        n (float): The value to clamp

    Returns:
        float: The clampled value 
    """
    return min(255.0,max(0.0, n))

def clamp_color(color : mu.Vector) -> mu.Vector:
    """Clamps the specified color passed as an argument.

    Args:
        color (mu.Vector): Triplet RGB representing the color.

    Returns:
        mu.Vector: Triplet RGB representing the clamped color
    """
    new_color = [minmax_color(elt) for elt in tuple(color)]
    return mu.Vector(tuple(new_color)) 

def direction_vector(_from:  mu.Vector , _to :  mu.Vector) -> mu.Vector:
    """
    Calculates the direction vector between two points.

    Args:
        _from (mu.Vector): The position of the start of the vector.
        _to (mu.Vector): The position the vector points to.

    Returns:
        mu.Vector: The normalized direction vector
    """
    return (_to - _from).normalized()

def lencol(lst: list) -> int:
    """
    Args:
        lst (list): 2 by 2 matrix

    Returns:
        int: Width of the matrix
    """
    if len(lst) > 0:
        return len(lst[0])
    return 0
      
def vec(vec : tuple):
    """Converts a tuple into a Vector"""
    return mu.Vector(vec)
