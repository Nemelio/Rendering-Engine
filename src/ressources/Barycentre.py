# Auteur : Arash HABIBI


'''
Dans cette version : 
La fonction intersect retourne aussi le numéro de la face intersectée.
Dans le main, on calcule les couleurs de chaque sommet de cette face 
et on fait une interpolation pour trouver la couleur courante. 
Ombrage de Gouraud.
'''


import bpy
import mathutils as mu
import random

#======================================================

def solutionSystEquationLin2D(vcol1, vcol2, v):
    '''
    vcol1, vcol2 et v : Vecteurs à 2 dimensions
    valeur de retour : un couple de flottants (x,y)
    Soit un système d'équations : 
    (vcol1, vcol2) * (x,y) = V
    la valeur de retour est la solution de cette 
    équation. Si ce système n'admet pas d'équation, 
    alors la valeur de retour est None.
    '''
    a = vcol1[0]
    b = vcol2[0]
    c = vcol1[1]
    d = vcol2[1]
    X = v[0]
    Y = v[1]
    det = a*d - b*c
    if det==0:
        return None
    else:
        A = d/det
        B = -b/det
        C = -c/det
        D = a/det
        x = A*X + B*Y
        y = C*X + D*Y
        return (x,y)
    
#======================================================

def coords2(M, A, B, C):
    '''
    M, A, B et C : un point (Vector)
    valeur de retour : un couple de flottants (a,b)
    de façon à ce que AM = a*AB + b*AC
    '''
    AM = M - A
    AB = B - A
    AC = C - A
    return solutionSystEquationLin2D(AB, AC, AM)
    
    
#======================================================

def coords3(M, A, B, C):
    '''
    M, A, B et C : un point (Vector)
    valeur de retour : un triplet de flottants (a,b,c)
    de façon à ce que M = a*A + b*B + c*C
    '''
    (alpha,beta) = coords2(M, A, B, C)
    return mu.Vector((1-alpha-beta, alpha, beta))

