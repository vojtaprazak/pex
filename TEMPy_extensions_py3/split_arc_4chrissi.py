from Vector import Vector, axis_angle_to_matrix
from transformations import euler_from_matrix
from numpy import array


# Uses vector maths to create inermediate points between two give points on an arc.
# v1 and v2 are Vector objects of your two points
# centre is the centre of your circle.
# n is the number of new intermediate points to create.
def split_arc_4chrissi(v1, v2, centre, n): 
    h = []
    angs = []
    # Subtract the centre so that vectors start from origin, for axis and angle calculations
    v1_c = v1-centre
    v2_c = v2-centre
    # The cross product of two vectors gives the line perpendicular to both (ie. the rotation axis between the two)
    axis = v1_c.cross(v2_c)
    # The normalised dot product (ie. the 'arg' function in the TEMPy Vector class) gives the angle between two vectors.
    angle = v1_c.arg(v2_c)/(n+1)
    for a in range(1, n+1):
        # Use the angle and axis to create a rotation matrix, because that is how TEMPy vectors are rotated.
        mat = axis_angle_to_matrix(axis[0], axis[1], axis[2], angle*a, rad=True)
        # Rotate v1_c by this matrix and add back the centre
        h.append(v1_c.matrix_transform(mat)+centre)
        # Convert the matrix to ZXZ euler rotation angles
        angs.append(euler_from_matrix(mat, axes='rzxz'))
    return h, angs

# An example
a = Vector(5,15,7)
b = Vector(5,3,19)
centre = Vector(5,3,7)
i, angs = split_arc_4chrissi(a,b,centre,10)
print(array(i))
print(array(angs))
