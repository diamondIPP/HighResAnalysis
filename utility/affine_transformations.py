# --------------------------------------------------------
#       methods to apply 2D affine transformations
# created on April 5th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from numpy import array, append, sin, cos, ones, dot
from numpy.linalg import inv


def scale_matrix(sx=1, sy=None):
    return array([[sx, 0, 0], [0, sx if sy is None else sy, 0], [0, 0, 1]])


def transition_matrix(ox=0, oy=0):
    return array([[1, 0, ox], [0, 1, oy], [0, 0, 1]])


def rotation_matrix(rx=0, ry=None):
    x = [cos(rx), -sin(rx), 0] if ry is None else append(rx[:2], 0)
    y = [sin(rx), cos(rx), 0] if ry is None else append(ry[:2], 0)
    return array([x, y, [0, 0, 1]])


def matrix_order(order):
    d = {'s': 0, 't': 1, 'r': 2}
    return array([d[i] for i in order])


def multiply(m0, m1, m2):
    """returns: M0 * M1 * M2 """
    return dot(m0, dot(m1, m2))


def matrix(sx=1, sy=1, ox=0, oy=0, rx=0, ry=None, order='srt'):
    s, t, r = scale_matrix(sx, sy), transition_matrix(ox, oy), rotation_matrix(rx, ry)
    return multiply(*array([s, t, r])[matrix_order(order)])


def transform(x=None, y=None, sx=1, sy=1, ox=0, oy=0, rx=0, ry=None, order='srt', invert=False):
    x = ones(y.size, 'i') if x is None else x
    y = ones(x.size, 'i') if y is None else y
    z = ones(x.size, 'i')
    m = matrix(sx, sy, ox, oy, rx, ry, order)
    return dot(inv(m) if invert else m, array([x, y, z]))[:2]
