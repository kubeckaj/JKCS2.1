from numpy import array, asarray, clip, dot, arccos, degrees, zeros_like
from numpy.linalg import norm
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation


def calculate_vector(coord1, coord2):
    return array(coord1) - array(coord2)


def atom_distance(atom1, atom2):
    return norm(array(atom2) - array(atom1))

def vector_length(vector):
    return norm(vector)

def normalize_vector(vector):
    norm_ = norm(vector)
    if norm_ < 1e-8:
        return zeros_like(vector)
    return vector / norm_


def rotate_vector(vector, axis, angle):
    # axis must be a unit vector (same contract as the Rodrigues formula this replaces)
    return Rotation.from_rotvec(angle * asarray(axis, dtype=float)).apply(vector)


def calculate_angle(coord1, coord2, coord3):
    vector1 = calculate_vector(coord2, coord1)
    vector2 = calculate_vector(coord2, coord3)
    cos_angle = dot(vector1, vector2) / (norm(vector1) * norm(vector2))
    return degrees(arccos(clip(cos_angle, -1.0, 1.0)))


def distance_matrix(coordinates):
    coords = asarray(coordinates, dtype=float)
    return cdist(coords, coords)
