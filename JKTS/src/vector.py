import numpy as np

class Vector:
    @staticmethod
    def calculate_vector(coord1, coord2):
        return np.array(coord1) - np.array(coord2)

    @staticmethod
    def vector_length(vector):
        return np.linalg.norm(vector)

    @staticmethod
    def atom_distance(atom1, atom2):
        return np.linalg.norm(np.array(atom2) - np.array(atom1))

    @staticmethod
    def normalize_vector(vector):
        norm_ = np.linalg.norm(vector)
        if norm_ < 1e-8:
            return np.zeros_like(vector)
        return vector / norm_

    @staticmethod
    def rotate_vector(vector, axis, angle):
        cos_theta = np.cos(angle)
        sin_theta = np.sin(angle)
        cross_product = np.cross(axis, vector)
        return (vector * cos_theta +
                cross_product * sin_theta +
                axis * np.dot(axis, vector) * (1 - cos_theta))

    @staticmethod
    def calculate_angle(coord1, coord2, coord3):
        vector1 = Vector.calculate_vector(coord2, coord1)
        vector2 = Vector.calculate_vector(coord2, coord3)

        dot_product = np.dot(vector1, vector2)
        magnitude1 = Vector.vector_length(vector1)
        magnitude2 = Vector.vector_length(vector2)

        angle_rad = np.arccos(dot_product / (magnitude1 * magnitude2))
        angle_deg = np.degrees(angle_rad)

        return angle_deg

    @staticmethod
    def get_upwards_perpendicular_axis(self, direction):
        # direction is the vector lying between 2 carbons in C=C
        pass