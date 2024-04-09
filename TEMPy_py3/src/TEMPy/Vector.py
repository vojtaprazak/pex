##===============================================================================
#     This file is part of TEMPy.
#     
#     TEMPy is a software designed to help the user in the manipulation 
#     and analyses of macromolecular assemblies using 3D electron microscopy maps. 
#     
#     Copyright 2010-2014 TEMPy Inventors and Birkbeck College University of London.
#                          The TEMPy Inventors are: Maya Topf, Daven Vasishtan, 
#                           Arun Prasad Pandurangan, Irene Farabella, Agnel-Praveen Joseph,
#                          Harpal Sahota
# 
# 
#     TEMPy is available under Public Licence.
#     
#     Please cite your use of TEMPy in published work:
#     
#     Vasishtan D, Topf M. (2011) J Struct Biol 174:333-343. Scoring functions for cryoEM density fitting.
#
#===============================================================================

from numpy import sqrt, matrix, random, array
from math import cos, sin, pi, acos, asin, atan2

class Vector:
    """A class representing Cartesian 3-dimensonal vectors."""
    
    def __init__(self, x,y,z):
        """x, y, z = Cartesian co-ordinates of vector."""
        self.x = x
        self.y = y
        self.z = z

    @classmethod
    def fromlist(cls, vec):
        """
        Create a vector from a Python iterable type.

        *vec*
            A Python iterable type of length 3, with numerical values.
        
        Return:
            A Vector instance with values taken from vec.

        """
        if len(vec) != 3:
            raise IndexError("Input value must have length 3!")
        return cls(float(vec[0]), float(vec[1]), float(vec[2]))
    
    def __repr__(self):
        return "(%.3f,%.3f,%.3f)" %(self.x, self.y, self.z)

    def __getitem__(self, index):
        l = [self.x, self.y, self.z]
        return l[index]

    def __iter__(self):
        l = [self.x, self.y, self.z]
        return l.__iter__()

    def copy(self):
        """
        Return:
            A copy of Vector instance
        """
        return Vector(self.x, self.y, self.z)
    
    def mod(self):
        """
        Return:
            The modulus (length) of the vector.
        """
        return sqrt(self.x**2 + self.y**2 + self.z**2)
    
    def reverse(self):
        """
        Flip the direction of a Vector instance.
        
        Return:
          A Vector instance  
        """
        return Vector(-self.x,-self.y,-self.z)
    
    def arg(self, vector):
        """Return the argument (angle) between this and another vector.RAD"""
        top  = self.dot(vector)
        bottom = self.mod()*vector.mod()
        #print 'top/bottom, top, bottom ', top/bottom, top, bottom
        if abs(top-bottom) < 0.00001:
            return 0.0
        else:
            #print 'top/bottom, top, bottom ', top/bottom, top, bottom
            #return acos(top/bottom)
            return acos(min(max(top/bottom,-1.0),1.0))
    
    def times(self, factor):
        """
        Multiplies a Vector instance by a scalar factor.
        
        Return:
          A Vector instance
        """
        return Vector(factor*self.x, factor*self.y, factor*self.z)
    
    def dot(self, vector):
        """
        Return:
            The dot product of this and another vector specified as input parameter. 
        """
        return vector.x * self.x + vector.y * self.y + vector.z * self.z
    
    def cross(self, vector):
        """
        Return:
            A Vector instance of the cross product of this and another vector specified as input parameter
        """
        newX = self.y*vector.z - self.z*vector.y
        newY = self.z*vector.x - self.x*vector.z
        newZ = self.x*vector.y - self.y*vector.x
        return Vector(newX, newY, newZ)
    
    def __sub__(self, vector):
        """Return a Vector instance of the subtraction of a vector from this one."""
        newX = self.x - vector.x
        newY = self.y - vector.y
        newZ = self.z - vector.z
        return Vector(newX, newY, newZ)

    def dist(self, vector):
        """
        Return:
            The distance between this and another vector specified as input parameter. 
        """
        return (self-vector).mod()
    
    def __add__(self, vector):
        """Return a Vector instance of the addition of a vector from this one."""
        newX = self.x + vector.x
        newY = self.y + vector.y
        newZ = self.z + vector.z
        return Vector(newX, newY, newZ)

    def __mul__(self, prod):
        newX = self.x*prod
        newY = self.y*prod
        newZ = self.z*prod
        return Vector(newX, newY, newZ)

    def __div__(self, divisor):
        newX = self.x/float(divisor)
        newY = self.y/float(divisor)
        newZ = self.z/float(divisor)
        return Vector(newX, newY, newZ)

    def __neg__(self):
        return Vector(-self.x, -self.y, -self.z)

    def translate(self, x, y, z):
        """
        Translate a Vector instance.
        
        Arguments:
            *x, y, z*
                distance in Angstroms in respective Cartesian directions to translate vector.
                
        Return:
            Vector instance.    
        """
        newX = self.x + x
        newY = self.y + y
        newZ = self.z + z
        return Vector(newX, newY, newZ)

    def matrix_transform(self, rot_mat):
        """
        Transform the vector using a transformation matrix.
        
        Arguments:
            *rot_mat*
                a 3x3 Python matrix instance.
        Return:
            A vector instance
        """
        vec_mat = matrix([[self.x],[self.y],[self.z]])
        new_pos = rot_mat*vec_mat
        x = float(new_pos[0])
        y = float(new_pos[1])
        z = float(new_pos[2])
        return Vector(x,y,z)

    def to_atom(self):
        """
        Create an Atom instance based on Vector instance.
        
        Return:
            Atom instance
        """
        from ProtRep_Biopy import BioPyAtom
        template = 'ATOM      1  C   NOR A   1      23.161  39.732 -25.038  1.00 10.00             C'
        a = BioPyAtom(template)
        a.x = self.x
        a.y = self.y
        a.z = self.z
        return a

    def to_array(self):
        return array([self.x,self.y,self.z])
    
    def unit(self):
        """
        Return:
            Vector instance of a unit vector.
        """
        mod = self.mod()
        if mod==0:
            return Vector(0,0,0)
        return Vector(self.x/mod, self.y/mod, self.z/mod)

###########################################################################
###########################################################################
###########################################################################
#### def out of the class . 
#### better have them separate as these definition are an adaptation of 
#### Transformations Python Module from Christoph Gohlke
#### http://www.lfd.uci.edu/~gohlke/
#### http://www.lfd.uci.edu/~gohlke/code/transformations.py.html
###########################################################################
###########################################################################
###########################################################################
   
def random_vector(min_v, max_v):
    """
    Generate a random vector.
    The values for the vector component x, y, and z are randomly sampled between minimum and maximum values specified.
    
    Argument:
        *min_v, max_v*
            minimum and maximum value
    Return:
        A Vector instance.
    """
    
    x = random.uniform(min_v, max_v)
    y = random.uniform(min_v, max_v)
    z = random.uniform(min_v, max_v)
    return Vector(x,y,z)

def axis_angle_to_matrix(x, y, z, turn, rad=False):
    """
    Converts the axis angle rotation to a matrix form.
    
    Arguments:
       *x, y, z*
           axis of rotation (does not need to be normalised).
       *turn*
           angle of rotation, in radians if rad=True, else in degrees.
    Return:
        A 3X3 transformation matrix.   
    """
    if not rad:
        turn = turn*pi/180
    c_a = cos(turn)
    s_a = sin(turn)
    v = Vector(x,y,z).unit()
    x = v.x
    y = v.y
    z = v.z

    rot_mat = matrix([[x**2+(1-x**2)*c_a, x*y*(1-c_a)-z*s_a, x*z*(1-c_a)+y*s_a],
                          [ x*y*(1-c_a)+z*s_a, y**2+(1-y**2)*c_a, y*z*(1-c_a)-x*s_a],
                          [x*z*(1-c_a)-y*s_a, y*z*(1-c_a)+x*s_a, z**2+(1-z**2)*c_a]])
    return rot_mat

def euler_to_matrix(x_turn, y_turn, z_turn, rad=False):
    """
    Converts an euler rotation to a matrix form.
    
    Arguments:
       *x_turn, y_turn, z_turn*
           rotation angles around respective axis, in radians if rad=True, else in degrees.
    Return:
        A 3X3 transformation matrix.   
    """
    if not rad:
        x_turn = x_turn*pi/180
        y_turn = y_turn*pi/180
        z_turn = z_turn*pi/180
    x_mat = axis_angle_to_matrix(0,0,1,x_turn, rad=True)
    y_mat = axis_angle_to_matrix(0,1,0,y_turn, rad=True)
    z_mat = axis_angle_to_matrix(1,0,0,z_turn, rad=True)
    return x_mat*y_mat*z_mat

def axis_angle_to_euler(x,y,z, turn, rad=False):
    """
    Converts the axis angle rotation to an Euler form.
    
    Arguments:
       *x, y, z*
           axis of rotation (does not need to be normalised).
       *turn*
           angle of rotation, in radians if rad=True, else in degrees.
    Returns:
        A 3-tuple (x,y,z) containing the Euler angles. .
    """
    if not rad:
        turn = turn*pi/180
    z_rot = atan2(y*sin(turn)-x*z*(1-cos(turn)), 1-(y**2+z**2)*(1-cos(turn)))
    x_rot = asin(x*y*(1-cos(turn))+z*sin(turn))
    y_rot = atan2(x*sin(turn)-y*z*(1-cos(turn)), 1-(x**2+z**2)*(1-cos(turn)))
    return (x_rot, y_rot, z_rot)

# -- Vector methods for torsion angle geometry -- #



def torsion(a, b, c):
    """
    Find the torsion angle between planes ab and bc.
    
    Arguments:
        
        *a,b,c*
            Vector instances.
    
    Returns:
        The torsion angle in radians
    """
    n1 = a.cross(b)
    n2 = b.cross(c)
    return n1.arg(n2)

def altTorsion(a,b,c):
    """
    An alternate and better way to find the torsion angle between planes ab and bc.
    
    Arguments:
        *a,b,c*
            Vector instances.
    Return:
        The torsion angle (radians)
    """
    A = a.dot(b.cross(c))*b.mod()
    B = (a.cross(b)).dot(b.cross(c))
    return atan2(A,B)
