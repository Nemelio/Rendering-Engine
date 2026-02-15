# Mardi 04 février 2025
# Arash HABIBI
# Fonction sur les images ppm

def creer(nom_fichier, largeur_pix, hauteur_pix):
    '''
    nom_fichier : string
    largeur_pix entier largeur de l'image en pixels 
    hauteur_pix : entier hauteur de l'image en pixels
    Cette fonction créé un fichier appelé nom_fichier
    et écrit son entete. 
    '''
    image = open(nom_fichier,"w+")
    image.write("P3\n")
    image.write(str(largeur_pix) + " " + str(hauteur_pix)+"\n")
    image.write("255\n")
    return image

#----------------------------

def writePixel(image, couleur):
    '''
    image : descripteur de fichier ppm
    couleur : Vector représentant la couleur à écrire.
    valeur de retour : aucune
    precondition : le fichier doit avoir ete ouvert.
    écrit ces trois valeurs dans le fichier
    '''
    image.write(str(int(couleur[0])) + " " + str(int(couleur[1])) + " " + str(int(couleur[2])) + "\n")
    
#----------------------------

def fini(image):
    image.close()
