import Torque_calculation_inputs
import Torque_calculation
import numpy as np
import math


Distance_between_two_parts = Torque_calculation_inputs.tube_diameter_1/2 + Torque_calculation_inputs.tube_diameter_2/2 + 0.001
print("Distance between two parts:", Distance_between_two_parts)

alpha = Torque_calculation_inputs.alpha
beta = Torque_calculation_inputs.beta

x = -420
y = -700
z = 0

Forces_c3 = np.array([x, y, z])
added_Moment_c3 = np.array([80, 0, 0])


def add_corrections_for_the_moment(z_moment_vector, force_vector_c3, alpha, beta):
    added_Moment_c3 = np.array([80, 0, 0])
    relative_vector_only_z_direction = np.array([0, 0, -Distance_between_two_parts])

    angle_between_alpha_and_beta_in_arc_below = alpha + beta

    shift_x2_to_x1 = Torque_calculation.create_z_axis_rotation_matrix(math.radians(angle_between_alpha_and_beta_in_arc_below))

    Joint_3_Moment_c3_in_x1 = np.matmul(shift_x2_to_x1, added_Moment_c3)
    Forces_c3_in_x1 = np.matmul(shift_x2_to_x1, force_vector_c3)
    Moments_due_to_forces_in_x1 = np.cross(relative_vector_only_z_direction, Forces_c3_in_x1)

    # print(angle_between_alpha_and_beta_in_arc_below)

    rotation_matrix_alpha = Torque_calculation.create_z_axis_rotation_matrix(alpha)
    x2_axis = np.array([1, 0, 0])

    # print(Joint_3_Moment_c3_in_x1)

    # print(relative_vector_only_z_direction)
    # print(Forces_c3_in_x1)
    # print("Moment vector:", Moments_due_to_forces_in_x1)

    Added_moment_in_x1 = Joint_3_Moment_c3_in_x1 + Moments_due_to_forces_in_x1


    # print("\n\n\n")
    # print(Added_moment_in_x1, "Nm")

    return z_moment_vector + Added_moment_in_x1

