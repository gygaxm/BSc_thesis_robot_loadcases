import numpy as np
import matplotlib.pyplot as plt
import math
import scipy
import Torque_calculation_inputs

alpha_rad = Torque_calculation_inputs.alpha/180 * math.pi
beta_rad = Torque_calculation_inputs.beta/180 * math.pi


"""
# Constants
lenth_per_segment = 0.14
radius_1 = 0.1
radius_2 = 0.1
joint_length_2 = 0.2
radius_3 = 0.1


number_of_segments_in_link_1 = 8
number_of_segments_in_link_2 = 8


# mass of the joint connection which encompasses the joint 2 & 3 in [kg]
mass_joint_2_3 = 10

# mass of the joint connection which encompasses the joint 2 & 3 in [kg]
mass_joint_4_5 = 5

# the safety factor for the holding torque without overheating the motors
safety_factor = 0.9



tube_diameter = 0.2
tube_wall_thickness = 0.002

end_effector_distance = 0.3
end_effector_mass = 6
"""

class Motor:
    def __init__(self, id, continuous, peak, mass=0):
        # the id shows which joint the motor is
        self.number = id

        # the continous attribute is the continuous torque [Nm] the motor can hold
        self.continuous = continuous

        # the peak torque [Nm]
        self.peak = peak

        # Not relevant, as the whole joint needs to be considered
        self.mass = mass
class Material:
    def __init__(self, name, e_modul, density):
        # the name of the material (like iron)
        self.name = name

        # e-modul in [Gpa]
        self.e_modul = e_modul

        # density in [kg/m^3]
        self.density = density

    def __str__(self):
        return f"{self.name}:\n  E-modul: {self.e_modul} GPa\n  density: {self.density} kg/m^3"
class Linkage:
    def __init__(self, index, material, mass=0, number_of_segments=0, length_per_segment = 0, tube_diameter=0.2, tube_wall_thickness = 0.002, length=0, safety_weight=0.25):
        self.index = index

        # the number of segments in this linkage
        # self.number_of_segments = number_of_segments

        # is the material (is a Material object)
        self.material = material

        # length in [m]
        self.length = length

        self.diameter = tube_diameter

        # mass in [kg]
        self.mass = math.pi * ((tube_diameter/2)**2 - (tube_diameter/2 - tube_wall_thickness)**2) * self.length * self.material.density + safety_weight
        # self.mass = mass

        self.gravity_force = self.mass*9.81

    def __str__(self):
        return f"Link {self.index}:\n  {self.number_of_segments} number of segments\n  material: {self.material.name}\n  length: {self.length} m\n  mass: {self.mass} kg"

gravity_constant = 9.81



# Motors
joint_0 = Motor(id=0, continuous=Torque_calculation_inputs.Motor_0['continuous'], peak=Torque_calculation_inputs.Motor_0['peak'])
joint_1 = Motor(id=1, continuous=Torque_calculation_inputs.Motor_1['continuous'], peak=Torque_calculation_inputs.Motor_1['peak'])
joint_2 = Motor(id=2, continuous=Torque_calculation_inputs.Motor_2['continuous'], peak=Torque_calculation_inputs.Motor_2['peak'])
joint_3 = Motor(id=3, continuous=Torque_calculation_inputs.Motor_3['continuous'], peak=Torque_calculation_inputs.Motor_3['peak'])
joint_4 = Motor(id=4, continuous=Torque_calculation_inputs.Motor_4['continuous'], peak=Torque_calculation_inputs.Motor_4['peak'])
joint_5 = Motor(id=5, continuous=Torque_calculation_inputs.Motor_5['continuous'], peak=Torque_calculation_inputs.Motor_5['peak'])


Steel = Material(name='Steel', e_modul=210, density=7850)

# inputs


link_1 = Linkage(1, length=Torque_calculation_inputs.length_link_1, material=Steel, tube_diameter=Torque_calculation_inputs.tube_diameter_1, tube_wall_thickness=Torque_calculation_inputs.tube_wall_thickness_1)
link_2 = Linkage(2, length=Torque_calculation_inputs.length_link_2, material=Steel, tube_diameter=Torque_calculation_inputs.tube_diameter_2, tube_wall_thickness=Torque_calculation_inputs.tube_wall_thickness_2)


"""
print(Steel)
print()
print(link_1)
print()
print(link_2)
"""


def add_corrections_for_the_moment(z_moment_vector, force_vector_c3, alpha, beta, relative_vector):
    # Distance_between_two_parts = Torque_calculation_inputs.tube_diameter_1 / 2 + Torque_calculation_inputs.tube_diameter_2 / 2 + 0.001
    added_Moment_c3 = np.array([80, 0, 0])

    # relative_vector_only_z_direction = np.array([0, 0, -Distance_between_two_parts])

    angle_between_alpha_and_beta_in_arc_below = alpha + beta

    shift_x2_to_x1 = create_z_axis_rotation_matrix(math.radians(angle_between_alpha_and_beta_in_arc_below))

    Joint_3_Moment_c3_in_x1 = np.matmul(shift_x2_to_x1, added_Moment_c3)
    Forces_c3_in_x1 = np.matmul(shift_x2_to_x1, force_vector_c3)
    Moments_due_to_forces_in_x1 = np.cross(relative_vector, Forces_c3_in_x1)

    # print(angle_between_alpha_and_beta_in_arc_below)

    """rotation_matrix_alpha = create_z_axis_rotation_matrix(alpha)
    x2_axis = np.array([1, 0, 0])

    # print(Joint_3_Moment_c3_in_x1)

    # print(relative_vector_only_z_direction)
    # print(Forces_c3_in_x1)
    # print("Moment vector:", Moments_due_to_forces_in_x1)"""

    Added_moment_in_x1 = Joint_3_Moment_c3_in_x1 + Moments_due_to_forces_in_x1

    # print("\n\n\n")
    # print(Added_moment_in_x1, "Nm")

    return z_moment_vector + Added_moment_in_x1
def calculate_the_position_of_the_center_of_mass(angle_alpha=0.0, angle_beta=0.0):
    # Positions of the center of mass of the objects from the base
    position_link_1 = np.array([math.cos(angle_alpha) * (link_1.length/2 + Torque_calculation_inputs.radius_1),
                                math.sin(angle_alpha) * (link_1.length/2 + Torque_calculation_inputs.radius_1),
                                0])
    # print(position_link_1)
    position_joint_2_3 = np.array([math.cos(angle_alpha) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2) + math.cos(angle_beta) * Torque_calculation_inputs.joint_length_2/2,
                                   math.sin(angle_alpha) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2) + math.sin(angle_beta) * Torque_calculation_inputs.joint_length_2/2,
                                   0])

    position_link_2 = np.array([math.cos(angle_alpha) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2) + math.cos(angle_beta) * (Torque_calculation_inputs.joint_length_2 + link_2.length/2),
                                math.sin(angle_alpha) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2) + math.sin(angle_beta) * (Torque_calculation_inputs.joint_length_2 + link_2.length/2),
                                0])

    position_joint_4_5 = np.array([math.cos(angle_alpha) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2) + math.cos(angle_beta) * (Torque_calculation_inputs.joint_length_2 + link_2.length + Torque_calculation_inputs.radius_3),
                                   math.sin(angle_alpha) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2) + math.sin(angle_beta) * (Torque_calculation_inputs.joint_length_2 + link_2.length + Torque_calculation_inputs.radius_3),
                                   0])

    position_end_effector = position_joint_4_5 + np.array([Torque_calculation_inputs.end_effector_distance, 0, 0])

    position_payload = position_end_effector + np.array([0, 0, 0])

    position_counterweight = np.array([math.cos(angle_alpha + math.pi)*Torque_calculation_inputs.counterweight_distance,
                                       math.sin(angle_alpha + math.pi) * Torque_calculation_inputs.counterweight_distance,
                                       0])

    return position_link_1, position_joint_2_3, position_link_2, position_joint_4_5, position_end_effector, position_payload, position_counterweight
def calculate_the_position_of_the_cutting_points(angle_alpha=0.0, angle_beta=0.0):
    # Positions of the center of mass of the objects from the base
    position_cut_1 = np.array([math.cos(angle_alpha) * (Torque_calculation_inputs.radius_1),
                                math.sin(angle_alpha) * (Torque_calculation_inputs.radius_1),
                                0])
    # print(position_link_1)
    position_cut_2 = np.array([math.cos(angle_alpha) * (link_1.length + Torque_calculation_inputs.radius_1),
                                math.sin(angle_alpha) * (link_1.length + Torque_calculation_inputs.radius_1),
                                0])

    position_cut_3 = np.array([math.cos(angle_alpha) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2) + math.cos(angle_beta) * (Torque_calculation_inputs.joint_length_2),
                                math.sin(angle_alpha) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2) + math.sin(angle_beta) * (Torque_calculation_inputs.joint_length_2),
                                0])

    position_cut_4 = np.array([math.cos(angle_alpha) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2) + math.cos(angle_beta) * (Torque_calculation_inputs.joint_length_2 + link_2.length),
                                math.sin(angle_alpha) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2) + math.sin(angle_beta) * (Torque_calculation_inputs.joint_length_2 + link_2.length),
                                0])

    return position_cut_1, position_cut_2, position_cut_3, position_cut_4
def calculate_relative_position_vector(new_origin, vector):
    return np.array([vector[0] - new_origin[0], vector[1]-new_origin[1], vector[2]-new_origin[2]])
def create_x_axis_cylinder_moment_of_Inertia(radius, height, mass):
    moment_of_inertia_tensor = np.array([[6 * (radius**2) * mass / 12, 0, 0],
                                         [0, mass/12 * (3 * radius**3 + height**2), 0],
                                         [0, 0, mass/12 * (3 * radius**3 + height**2)]])
    return moment_of_inertia_tensor
def create_x_axis_tube_moment_of_Inertia(radius, height, mass):
    moment_of_inertia_tensor = np.array([[(radius**2) * mass, 0, 0],
                                         [0, mass/12 * (6 * radius**2 + height**2), 0],
                                         [0, 0, mass/12 * (6 * radius**2 + height**2)]])
    return moment_of_inertia_tensor
def create_rectangular_prism_moment_of_Inertia(x, y, z, mass):
    moment_of_inertia_tensor = np.array([[mass/12 * (y**2 + z**2), 0, 0],
                                         [0, mass/12 * (x**2 + z**2), 0],
                                         [0, 0, mass/12 * (x**2 + y**2)]])
    return moment_of_inertia_tensor
def create_z_axis_rotation_matrix(angle):
    rotation_array = np.array([[math.cos(angle), -math.sin(angle), 0],
                              [math.sin(angle), math.cos(angle), 0],
                              [0, 0, 1]])
    return rotation_array
def calculate_shifting_tensor_steiner(shifting_vector, mass):
    x = shifting_vector[0]
    y = shifting_vector[1]
    z = shifting_vector[2]
    shifting_tensor = np.array([[y**2 + z**2, -x * y, - x * z],
                                [-x*y, x**2 + z**2, -y * z],
                                [-x * z, - y * z, x**2 + y**2]])
    shifting_tensor = shifting_tensor * mass
    return shifting_tensor
def check_if_static_torque_in_safety(static_torque, motor:Motor, safety_factor):
    if static_torque < (motor.continuous * safety_factor):
        return True
    return False
def calculate_forces_for_joint_1_acceleration(positions, masses, omega, alpha=0.0, beta=0.0):
    if omega[2] > Torque_calculation_inputs.Max_angular_acceleration_joint_1:
        omega = np.array([0, 0, Torque_calculation_inputs.Max_angular_acceleration_joint_1])

    position_cut_1, position_cut_2, position_cut_3, position_cut_4 = calculate_the_position_of_the_cutting_points(alpha, beta)
    force_c1 = np.zeros(3)
    z_momentum_c1 = np.zeros(3)
    for i in range(len(positions)):
        linear_acceleration = np.cross(omega, positions[i])
        force = linear_acceleration * masses[i]

        gravity_force = np.array([0, -masses[i] * gravity_constant, 0])
        # print(force, gravity_force)
        total_force_for_part = force - gravity_force
        force_c1 += total_force_for_part
        relative_cut_position = calculate_relative_position_vector(position_cut_1, positions[i])
        z_momentum_c1 += np.cross(relative_cut_position, total_force_for_part)
        # z_momentum_c1 += np.cross(positions[i], total_force_for_part)

    # print("C1:", force_c1, "   C1-momentum:", z_momentum_c1)

    del positions[0]
    del masses[0]

    force_c2 = np.zeros(3)
    z_momentum_c2 = np.zeros(3)
    for i in range(len(positions)):
        linear_acceleration = np.cross(omega, positions[i])
        force = linear_acceleration * masses[i]
        gravity_force = np.array([0, -masses[i] * gravity_constant, 0])
        total_force_for_part = force - gravity_force
        force_c2 += total_force_for_part
        relative_cut_position = calculate_relative_position_vector(position_cut_1, positions[i])
        z_momentum_c2 += np.cross(relative_cut_position, total_force_for_part)
        # z_momentum_c2 += np.cross(positions[i], total_force_for_part)
    # print("C2:", force_c2, "   C2-momentum:", z_momentum_c2)

    del positions[0]
    del masses[0]

    force_c3 = np.zeros(3)
    z_momentum_c3 = np.zeros(3)
    for i in range(len(positions)):
        linear_acceleration = np.cross(omega, positions[i])
        force = linear_acceleration * masses[i]
        gravity_force = np.array([0, -masses[i] * gravity_constant, 0])
        total_force_for_part = force - gravity_force
        force_c3 += total_force_for_part
        relative_cut_position = calculate_relative_position_vector(position_cut_3, positions[i])
        z_momentum_c3 += np.cross(relative_cut_position, total_force_for_part)
        # z_momentum_c3 += np.cross(positions[i], total_force_for_part)
    # print("C3:", force_c3)
    # print("C3:", force_c3, "   C3-momentum:", z_momentum_c3)

    del positions[0]
    del masses[0]

    force_c4 = np.zeros(3)
    z_momentum_c4 = np.zeros(3)
    for i in range(len(positions)):
        linear_acceleration = np.cross(omega, positions[i])
        force = linear_acceleration * masses[i]
        gravity_force = np.array([0, -masses[i] * gravity_constant, 0])
        total_force_for_part = force - gravity_force
        force_c4 += total_force_for_part
        relative_cut_position = calculate_relative_position_vector(position_cut_4, positions[i])
        z_momentum_c4 += np.cross(relative_cut_position, total_force_for_part)
        # z_momentum_c4 += np.cross(positions[i], total_force_for_part)
    # print("C4:", force_c4)
    """print("C1:", force_c1, "   C1-momentum:", z_momentum_c1)
    print("C2:", force_c2, "   C2-momentum:", z_momentum_c2)
    print("C3:", force_c3, "   C3-momentum:", z_momentum_c3)
    print("C4:", force_c4, "   C4-momentum:", z_momentum_c4)"""

    # return {'C1': force_c1, 'C2': force_c2, 'C3': force_c3, 'C4': force_c4}
    return {'C1': force_c1, 'C1-momentum': z_momentum_c1, 'C2': force_c2, 'C2-momentum': z_momentum_c2,
            'C3': force_c3, 'C3-momentum': z_momentum_c3, 'C4': force_c4, 'C4-momentum': z_momentum_c4}
def calculate_forces_for_joint_2_acceleration(positions, positions_joint_2, masses, omega, alpha=0.0, beta=0.0):
    if omega[2] > Torque_calculation_inputs.Max_angular_acceleration_joint_2:
        omega = np.array([0, 0, Torque_calculation_inputs.Max_angular_acceleration_joint_2])

    position_cut_1, position_cut_2, position_cut_3, position_cut_4 = calculate_the_position_of_the_cutting_points(alpha, beta)
    z_momentum_c1 = np.zeros(3)
    force_c1 = np.zeros(3)
    force_c1 += np.array([0, masses[0] * gravity_constant, 0])
    z_momentum_c1 += np.cross(calculate_relative_position_vector(position_cut_1, positions[0]), np.array([0, masses[0] * gravity_constant, 0]))
    for i in range(1, len(positions)):
        relative_position = calculate_relative_position_vector(positions_joint_2, positions[i])
        linear_acceleration = np.cross(omega, relative_position)
        force = linear_acceleration * masses[i]

        gravity_force = np.array([0, -masses[i] * gravity_constant, 0])
        # print(force, gravity_force)
        total_force_for_part = force - gravity_force
        force_c1 += total_force_for_part
        relative_cut_position = calculate_relative_position_vector(position_cut_1, positions[i])
        z_momentum_c1 += np.cross(relative_cut_position, total_force_for_part)

    # print("C1:", force_c1)

    del positions[0]
    del masses[0]

    force_c2 = np.zeros(3)
    z_momentum_c2 = np.zeros(3)
    for i in range(len(positions)):
        relative_position = calculate_relative_position_vector(positions_joint_2, positions[i])
        linear_acceleration = np.cross(omega, relative_position)
        force = linear_acceleration * masses[i]
        gravity_force = np.array([0, -masses[i] * gravity_constant, 0])
        total_force_for_part = force - gravity_force
        force_c2 += total_force_for_part
        relative_cut_position = calculate_relative_position_vector(position_cut_2, positions[i])
        z_momentum_c2 += np.cross(relative_cut_position, total_force_for_part)
        # z_momentum_c2 += np.cross(relative_position, total_force_for_part)
    # print("C2:", force_c2)

    del positions[0]
    del masses[0]

    force_c3 = np.zeros(3)
    z_momentum_c3 = np.zeros(3)
    for i in range(len(positions)):
        relative_position = calculate_relative_position_vector(positions_joint_2, positions[i])
        linear_acceleration = np.cross(omega, relative_position)
        force = linear_acceleration * masses[i]
        gravity_force = np.array([0, -masses[i] * gravity_constant, 0])
        total_force_for_part = force - gravity_force
        force_c3 += total_force_for_part
        relative_cut_position = calculate_relative_position_vector(position_cut_3, positions[i])
        z_momentum_c3 += np.cross(relative_cut_position, total_force_for_part)
        # z_momentum_c3 += np.cross(relative_position, total_force_for_part)
    # print("C3:", force_c3)

    del positions[0]
    del masses[0]

    force_c4 = np.zeros(3)
    z_momentum_c4 = np.zeros(3)
    for i in range(len(positions)):
        relative_position = calculate_relative_position_vector(positions_joint_2, positions[i])
        linear_acceleration = np.cross(omega, relative_position)
        force = linear_acceleration * masses[i]
        gravity_force = np.array([0, -masses[i] * gravity_constant, 0])
        total_force_for_part = force - gravity_force
        force_c4 += total_force_for_part
        relative_cut_position = calculate_relative_position_vector(position_cut_4, positions[i])
        z_momentum_c4 += np.cross(relative_cut_position, total_force_for_part)
        # z_momentum_c4 += np.cross(relative_position, total_force_for_part)
    # print("C4:", force_c4)

    return {'C1': force_c1, 'C1-momentum': z_momentum_c1, 'C2': force_c2, 'C2-momentum': z_momentum_c2,
            'C3': force_c3, 'C3-momentum': z_momentum_c3, 'C4': force_c4, 'C4-momentum': z_momentum_c4}
    # return {'C1': force_c1, 'C2': force_c2, 'C3': force_c3, 'C4': force_c4}
def calculate_forces_for_joint_4_acceleration(positions, positions_joint_4, masses, omega, alpha=0.0, beta=0.0):
    if omega[2] > Torque_calculation_inputs.Max_angular_acceleration_joint_4:
        omega = np.array([0, 0, Torque_calculation_inputs.Max_angular_acceleration_joint_4])

    position_cut_1, position_cut_2, position_cut_3, position_cut_4 = calculate_the_position_of_the_cutting_points(alpha, beta)
    force_c1 = np.zeros(3)
    z_momentum_c1 = np.zeros(3)
    for i in range(len(positions)-3):
        total_force = np.array([0, masses[i] * gravity_constant, 0])
        relative_position = calculate_relative_position_vector(positions_joint_4, positions[i])
        force_c1 += total_force
        relative_cut_position = calculate_relative_position_vector(position_cut_1, positions[i])
        z_momentum_c1 += np.cross(relative_cut_position, total_force)
        # z_momentum_c1 += np.cross(relative_position, total_force)
    for i in range(len(positions)-2, len(positions)):
        relative_position = calculate_relative_position_vector(positions_joint_4, positions[i])
        linear_acceleration = np.cross(omega, relative_position)
        force = linear_acceleration * masses[i]

        gravity_force = np.array([0, -masses[i] * gravity_constant, 0])
        # print(force, gravity_force)
        total_force_for_part = force - gravity_force
        force_c1 += total_force_for_part

        relative_cut_position = calculate_relative_position_vector(position_cut_1, positions[i])
        z_momentum_c1 += np.cross(relative_cut_position, total_force_for_part)
        # z_momentum_c1 += np.cross(relative_position, total_force_for_part)


    # print("C1:", force_c1)

    del positions[0]
    del masses[0]

    force_c2 = np.zeros(3)
    z_momentum_c2 = np.zeros(3)
    for i in range(len(positions)-3):
        total_force = np.array([0, masses[i] * gravity_constant, 0])
        relative_position = calculate_relative_position_vector(positions_joint_4, positions[i])
        force_c2 += total_force
        relative_cut_position = calculate_relative_position_vector(position_cut_2, positions[i])
        z_momentum_c2 += np.cross(relative_cut_position, total_force)
        # z_momentum_c2 += np.cross(relative_position, total_force)
    for i in range(len(positions)-2, len(positions)):
        relative_position = calculate_relative_position_vector(positions_joint_4, positions[i])
        linear_acceleration = np.cross(omega, relative_position)
        force = linear_acceleration * masses[i]
        gravity_force = np.array([0, -masses[i] * gravity_constant, 0])
        total_force_for_part = force - gravity_force
        force_c2 += total_force_for_part

        relative_cut_position = calculate_relative_position_vector(position_cut_2, positions[i])
        z_momentum_c2 += np.cross(relative_cut_position, total_force_for_part)
        # z_momentum_c2 += np.cross(relative_position, total_force_for_part)
    # print("C2:", force_c2)

    del positions[0]
    del masses[0]

    force_c3 = np.zeros(3)
    z_momentum_c3 = np.zeros(3)
    for i in range(len(positions)-3):
        total_force = np.array([0, masses[i] * gravity_constant, 0])
        relative_position = calculate_relative_position_vector(positions_joint_4, positions[i])
        force_c3 += total_force
        relative_cut_position = calculate_relative_position_vector(position_cut_3, positions[i])
        z_momentum_c3 += np.cross(relative_cut_position, total_force)
        # z_momentum_c3 += np.cross(relative_position, total_force)
    for i in range(len(positions)-2, len(positions)):
        relative_position = calculate_relative_position_vector(positions_joint_4, positions[i])
        linear_acceleration = np.cross(omega, relative_position)
        force = linear_acceleration * masses[i]
        gravity_force = np.array([0, -masses[i] * gravity_constant, 0])
        total_force_for_part = force - gravity_force
        force_c3 += total_force_for_part

        relative_cut_position = calculate_relative_position_vector(position_cut_3, positions[i])
        z_momentum_c3 += np.cross(relative_cut_position, total_force_for_part)
        # z_momentum_c3 += np.cross(relative_position, total_force_for_part)
    # print("C3:", force_c3)

    del positions[0]
    del masses[0]

    force_c4 = np.zeros(3)
    z_momentum_c4 = np.zeros(3)
    for i in range(len(positions)-3):
        total_force = np.array([0, masses[i] * gravity_constant, 0])
        relative_position = calculate_relative_position_vector(positions_joint_4, positions[i])
        force_c4 += total_force

        relative_cut_position = calculate_relative_position_vector(position_cut_4, positions[i])
        z_momentum_c4 += np.cross(relative_cut_position, total_force)
        # z_momentum_c4 += np.cross(relative_position, total_force)
    for i in range(len(positions)-2, len(positions)):
        relative_position = calculate_relative_position_vector(positions_joint_4, positions[i])
        linear_acceleration = np.cross(omega, relative_position)
        force = linear_acceleration * masses[i]
        gravity_force = np.array([0, -masses[i] * gravity_constant, 0])
        total_force_for_part = force - gravity_force
        force_c4 += total_force_for_part

        relative_cut_position = calculate_relative_position_vector(position_cut_4, positions[i])
        z_momentum_c4 += np.cross(relative_cut_position, total_force_for_part)
        # z_momentum_c4 += np.cross(relative_position, total_force_for_part)
    # print("C4:", force_c4)

    return {'C1': force_c1, 'C1-momentum': z_momentum_c1, 'C2': force_c2, 'C2-momentum': z_momentum_c2,
            'C3': force_c3, 'C3-momentum': z_momentum_c3, 'C4': force_c4, 'C4-momentum': z_momentum_c4}


# POST: return dictionary with keys [1,2,4] for the specific joints
def calculate_the_static_torque(angle_alpha=0.0, angle_beta=0.0):
    gravity_force_joint_2_3 = Torque_calculation_inputs.mass_joint_2_3 * 9.81
    gravity_force_joint_4_5 = Torque_calculation_inputs.mass_joint_4_5 * 9.81
    gravity_force_end_effector = Torque_calculation_inputs.end_effector_mass * 9.81
    gravity_force_payload = Torque_calculation_inputs.payload * 9.81
    gravity_force_counterweight = Torque_calculation_inputs.counterweight * 9.81


    # Positions of the center of mass of the objects from the base
    position_link_1, position_joint_2_3, position_link_2, position_joint_4_5, position_end_effector, position_payload, position_counterweight = calculate_the_position_of_the_center_of_mass(angle_alpha, angle_beta)


    static_torque_joint_1 = link_1.gravity_force * position_link_1[0] + gravity_force_joint_2_3 * position_joint_2_3[0] +\
                            link_2.gravity_force * position_link_2[0] + gravity_force_joint_4_5 * position_joint_4_5[0] +\
                            gravity_force_end_effector * position_end_effector[0] + gravity_force_payload * position_payload[0] +\
                            gravity_force_counterweight * position_counterweight[0]

    # print(static_torque_joint_1)


    position_joint_2 = np.array([math.cos(angle_alpha) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2),
                                 math.sin(angle_alpha) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2),
                                 0])

    relative_position_joint_2_to_link_2 = calculate_relative_position_vector(position_joint_2, position_link_2)
    relative_position_joint_2_to_joint_4_5 = calculate_relative_position_vector(position_joint_2, position_joint_4_5)
    relative_position_joint_2_to_end_effector = calculate_relative_position_vector(position_joint_2, position_end_effector)
    relative_position_joint_2_to_payload = calculate_relative_position_vector(position_joint_2, position_payload)

    static_torque_joint_2 = link_2.gravity_force * relative_position_joint_2_to_link_2[0] + gravity_force_joint_4_5 * relative_position_joint_2_to_joint_4_5[0] +\
                            gravity_force_end_effector * relative_position_joint_2_to_end_effector[0] + gravity_force_payload * relative_position_joint_2_to_payload[0]

    # print(static_torque_joint_2)


    relative_position_joint_4_to_end_effector = calculate_relative_position_vector(position_joint_4_5, position_end_effector)
    relative_position_joint_4_to_payload = calculate_relative_position_vector(position_joint_4_5, position_payload)

    static_torque_joint_4 = gravity_force_end_effector * relative_position_joint_4_to_end_effector[0] +\
                            gravity_force_payload * relative_position_joint_4_to_payload[0]

    return {1: static_torque_joint_1, 2: static_torque_joint_2, 4: static_torque_joint_4}

# POST: return dictionary with keys [1,2,4] for the specific origin joint
def calculate_moment_of_inertia_tensors(angle_alpha=0.0, angle_beta=0.0):
    position_link_1, position_joint_2_3, position_link_2, position_joint_4_5, position_end_effector, position_payload, position_counterweight = calculate_the_position_of_the_center_of_mass(angle_alpha, angle_beta)

    positions = [position_counterweight, position_link_1, position_joint_2_3, position_link_2, position_joint_4_5, position_end_effector, position_payload]
    masses = [Torque_calculation_inputs.counterweight, link_1.mass, Torque_calculation_inputs.mass_joint_2_3, link_2.mass, Torque_calculation_inputs.mass_joint_4_5, Torque_calculation_inputs.end_effector_mass, Torque_calculation_inputs.payload]

    tensors = [np.zeros((3, 3)),
               create_x_axis_tube_moment_of_Inertia(radius=0.1, height=link_1.length, mass=link_1.mass),
               create_x_axis_cylinder_moment_of_Inertia(radius=0.1, height=Torque_calculation_inputs.joint_length_2, mass=Torque_calculation_inputs.mass_joint_2_3),
               create_x_axis_tube_moment_of_Inertia(radius=0.1, height=link_2.length, mass=link_2.mass),
               create_x_axis_cylinder_moment_of_Inertia(radius=0.1, height=Torque_calculation_inputs.end_effector_distance/2, mass=Torque_calculation_inputs.mass_joint_4_5),
               create_x_axis_tube_moment_of_Inertia(radius=0.1, height=Torque_calculation_inputs.end_effector_distance/2, mass=Torque_calculation_inputs.end_effector_mass),
               np.zeros((3, 3))]

    rotation_matrices = [create_z_axis_rotation_matrix(angle_alpha + math.pi),
                         create_z_axis_rotation_matrix(angle_alpha),
                         create_z_axis_rotation_matrix(angle_beta),
                         create_z_axis_rotation_matrix(angle_beta),
                         create_z_axis_rotation_matrix(0),
                         create_z_axis_rotation_matrix(0),
                         create_z_axis_rotation_matrix(0)]



    # here the rotated tensors are created
    rotated_tensors = []
    for i in range(len(tensors)):
        rotated_tensors.append(np.matmul(tensors[i], rotation_matrices[i]))


    # these are purely to shift the tensor (steiger theorem part)
    shifting_tensors_for_base = []
    for i in range(len(positions)):
        shifting_tensors_for_base.append(calculate_shifting_tensor_steiner(positions[i], masses[i]))

    # these are the rotated & shifted tensors
    moment_of_inertia_tensors_from_base = []
    for i in range(len(rotated_tensors)):
        moment_of_inertia_tensors_from_base.append(rotated_tensors[i] + shifting_tensors_for_base[i])

    # here is the complete tensor sum of all part tensors
    moment_of_inertia_from_base = 0
    for i in range(0, len(moment_of_inertia_tensors_from_base)):
        # print(moment_of_inertia_tensors_from_base[i])
        # print()
        moment_of_inertia_from_base = moment_of_inertia_from_base + moment_of_inertia_tensors_from_base[i]
    # print()
    # print(moment_of_inertia_from_base)

    """for i in range(5):
        print()"""


    # take out the part that is no longer important
    del positions[0]
    del masses[0]
    del rotated_tensors[0]

    del positions[0]
    del masses[0]
    del rotated_tensors[0]

    position_joint_2 = np.array([math.cos(angle_alpha) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2),
                                 math.sin(angle_alpha) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2),
                                 0])

    # print(positions)
    shifting_tensors_for_joint_2 = []
    for i in range(len(positions)):
        shifting_tensors_for_joint_2.append(calculate_shifting_tensor_steiner(calculate_relative_position_vector(position_joint_2, positions[i]), masses[i]))

    # these are the rotated & shifted tensors for joint 2
    moment_of_inertia_tensors_from_joint_2 = []
    for i in range(len(rotated_tensors)):
        moment_of_inertia_tensors_from_joint_2.append(rotated_tensors[i] + shifting_tensors_for_joint_2[i])

    # here is the complete tensor sum of all part tensors
    moment_of_inertia_from_joint_2 = 0
    for i in range(0, len(moment_of_inertia_tensors_from_joint_2)):
        moment_of_inertia_from_joint_2 = moment_of_inertia_from_joint_2 + moment_of_inertia_tensors_from_joint_2[i]
    # print()
    # print(moment_of_inertia_from_joint_2)



    # take out the part that is no longer important
    del positions[0]
    del masses[0]
    del rotated_tensors[0]

    del positions[0]
    del masses[0]
    del rotated_tensors[0]

    position_joint_4 = position_joint_4_5

    shifting_tensors_for_joint_4 = []
    for i in range(len(positions)):
        shifting_tensors_for_joint_4.append(calculate_shifting_tensor_steiner(calculate_relative_position_vector(position_joint_4, positions[i]), masses[i]))

    # these are the rotated & shifted tensors for joint 2
    moment_of_inertia_tensors_from_joint_4 = []
    for i in range(len(rotated_tensors)):
        moment_of_inertia_tensors_from_joint_4.append(rotated_tensors[i] + shifting_tensors_for_joint_4[i])

    # here is the complete tensor sum of all part tensors
    moment_of_inertia_from_joint_4 = 0
    for i in range(0, len(moment_of_inertia_tensors_from_joint_4)):
        moment_of_inertia_from_joint_4 = moment_of_inertia_from_joint_4 + moment_of_inertia_tensors_from_joint_4[i]
    # print()
    # print(moment_of_inertia_from_joint_4)

    return {1: moment_of_inertia_from_base, 2: moment_of_inertia_from_joint_2, 4: moment_of_inertia_from_joint_4}


def calculate_the_maximum_angular_acceleration(angle_alpha=0.0, angle_beta=0.0, peak=False):
    inertia_tensors = calculate_moment_of_inertia_tensors(angle_alpha, angle_beta)
    static_torques = calculate_the_static_torque(angle_alpha, angle_beta)
    if peak==True:
        motor_torques = {1:joint_1.peak, 2:joint_2.peak, 4:joint_4.peak}
    else:
        motor_torques = {1:joint_1.continuous, 2:joint_2.continuous, 4:joint_4.continuous}

    free_torques = {1: np.array([0, 0, (motor_torques[1] - static_torques[1])]), 2: np.array([0, 0, (motor_torques[2] - static_torques[2])]), 4: np.array([0, 0, (motor_torques[4] - static_torques[4])])}
    # print(free_torques)
    angular_acceleration = {1: (np.matmul(np.linalg.inv(inertia_tensors[1]), free_torques[1])),
                            2: (np.matmul(np.linalg.inv(inertia_tensors[2]), free_torques[2])),
                            4: (np.matmul(np.linalg.inv(inertia_tensors[4]), free_torques[4]))}

    return angular_acceleration

# The x- & y-moment are not yet added
def calculate_force_on_linkage_in_global_coordinate_system(angle_alpha=0.0, angle_beta=0.0):
    position_link_1, position_joint_2_3, position_link_2, position_joint_4_5, position_end_effector, position_payload, position_counterweight = calculate_the_position_of_the_center_of_mass(angle_alpha, angle_beta)
    positions = [position_link_1, position_joint_2_3, position_link_2, position_joint_4_5, position_end_effector, position_payload]
    masses = [link_1.mass, Torque_calculation_inputs.mass_joint_2_3, link_2.mass, Torque_calculation_inputs.mass_joint_4_5, Torque_calculation_inputs.end_effector_mass, Torque_calculation_inputs.payload]
    position_joint_2 = np.array([math.cos(angle_alpha) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2),
                                 math.sin(angle_alpha) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2),
                                 0])

    angular_accelerations = calculate_the_maximum_angular_acceleration(angle_alpha, angle_beta, peak=True)

    force_on_linkage_ends_joint_1_acceleration = calculate_forces_for_joint_1_acceleration(positions.copy(), masses.copy(), angular_accelerations[1])
    # print(f"Forces for joint 1 accelerations\nC1: {force_on_linkage_ends_joint_1_acceleration['C1']}\nC2: {force_on_linkage_ends_joint_1_acceleration['C2']}\nC3: {force_on_linkage_ends_joint_1_acceleration['C3']}\nC4: {force_on_linkage_ends_joint_1_acceleration['C4']}\n")
    # print()
    force_on_linkage_ends_joint_2_acceleration = calculate_forces_for_joint_2_acceleration(positions.copy(), position_joint_2, masses.copy(), angular_accelerations[2])
    # print(f"Forces for joint 2 accelerations\nC1: {force_on_linkage_ends_joint_2_acceleration['C1']}\nC2: {force_on_linkage_ends_joint_2_acceleration['C2']}\nC3: {force_on_linkage_ends_joint_2_acceleration['C3']}\nC4: {force_on_linkage_ends_joint_2_acceleration['C4']}\n")
    # print()
    force_on_linkage_ends_joint_4_acceleration = calculate_forces_for_joint_4_acceleration(positions.copy(), position_joint_4_5, masses.copy(), angular_accelerations[2])
    # print(f"Forces for joint 4 accelerations\nC1: {force_on_linkage_ends_joint_4_acceleration['C1']}\nC2: {force_on_linkage_ends_joint_4_acceleration['C2']}\nC3: {force_on_linkage_ends_joint_4_acceleration['C3']}\nC4: {force_on_linkage_ends_joint_4_acceleration['C4']}\n")



    force_on_linkage_ends = {1: force_on_linkage_ends_joint_1_acceleration, 2: force_on_linkage_ends_joint_2_acceleration, 4: force_on_linkage_ends_joint_4_acceleration}

    # force_for_acceleration = {1: force_on_linkage_ends, 2: force_on_linkage_ends, 4: force_on_linkage_ends}

    return force_on_linkage_ends

# The y- & y-moment are added in this function
def calculate_force_on_linkage_in_linkage_coordinate_systems(angle_alpha=0.0, angle_beta=0.0):
    forces = calculate_force_on_linkage_in_global_coordinate_system(angle_alpha, angle_beta)

    """test_rotatation = create_z_axis_rotation_matrix(-0.1)
    test_vector = np.array([1,0,0])
    print(np.matmul(test_rotatation, test_vector))"""

    alpha_rotation = create_z_axis_rotation_matrix(angle_alpha)
    beta_rotation = create_z_axis_rotation_matrix(angle_beta)

    force = np.zeros(3)
    adjusted_forces = {1: {'C1': force, 'C1-momentum': force, 'C2': force, 'C2-momentum': force, 'C3': force, 'C3-momentum': force, 'C4': force, 'C4-momentum': force},
                       2: {'C1': force, 'C1-momentum': force, 'C2': force, 'C2-momentum': force, 'C3': force, 'C3-momentum': force, 'C4': force, 'C4-momentum': force},
                       4: {'C1': force, 'C1-momentum': force, 'C2': force, 'C2-momentum': force, 'C3': force, 'C3-momentum': force, 'C4': force, 'C4-momentum': force}}
    # print(forces)


    for k in forces:
        v = forces[k]
        adjusted_forces[k]['C4'] = np.matmul(beta_rotation, forces[k]["C4"])
        adjusted_forces[k]["C4-momentum"] = np.matmul(alpha_rotation, forces[k]["C4-momentum"])
        adjusted_forces[k]['C3'] = np.matmul(beta_rotation, forces[k]["C3"])
        adjusted_forces[k]["C3-momentum"] = np.matmul(alpha_rotation, forces[k]["C3-momentum"])


        relative_vector_1 = np.array([0, 0, Torque_calculation_inputs.tube_diameter_1 / 2 + Torque_calculation_inputs.tube_diameter_2 / 2 + 0.001])
        adjusted_forces[k]['C2'] = np.matmul(alpha_rotation, forces[k]["C2"])
        adjusted_forces[k]["C2-momentum"] = add_corrections_for_the_moment(np.matmul(alpha_rotation, forces[k]["C2-momentum"]), adjusted_forces[k]['C3'], math.degrees(angle_alpha), math.degrees(-angle_beta), relative_vector_1)

        relative_vector_2 = np.array([Torque_calculation_inputs.length_link_1, 0, Torque_calculation_inputs.tube_diameter_1 / 2 + Torque_calculation_inputs.tube_diameter_2 / 2 + 0.001])
        adjusted_forces[k]['C1'] = np.matmul(alpha_rotation, forces[k]["C1"])
        adjusted_forces[k]["C1-momentum"] = add_corrections_for_the_moment(np.matmul(alpha_rotation, forces[k]["C1-momentum"]), adjusted_forces[k]['C3'], math.degrees(angle_alpha), math.degrees(-angle_beta), relative_vector_2)


        # print(v)

    return adjusted_forces


def calculate_max_necessary_emergency_brake_energy(angular_velocity):
    moment_of_inertia_tensors = calculate_moment_of_inertia_tensors()
    moment_from_base = moment_of_inertia_tensors[1]
    # moment_from_joint_2 = moment_of_inertia_tensors[2]
    # moment_from_joint_4 = moment_of_inertia_tensors[4]


    energy = 1/2 * np.matmul(np.matmul(moment_from_base, angular_velocity), angular_velocity)

    return energy


def create_display_depending_on_angles(alpha_rad, beta_rad, axis, title):
    # print(alpha_rad)
    position_joint_2 = np.array([math.cos(alpha_rad) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2),
                                     math.sin(alpha_rad) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2),
                                     0])
    connection_points = [np.array([0, 0, 0]),
                         np.array([math.cos(alpha_rad) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2),
                                     math.sin(alpha_rad) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2),
                                     0]),
                         np.array([math.cos(alpha_rad) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2) + math.cos(beta_rad) * (Torque_calculation_inputs.joint_length_2 + link_2.length + Torque_calculation_inputs.radius_3),
                                       math.sin(alpha_rad) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2) + math.sin(beta_rad) * (Torque_calculation_inputs.joint_length_2 + link_2.length + Torque_calculation_inputs.radius_3),
                                       0]),
                         np.array([math.cos(alpha_rad) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2) + math.cos(beta_rad) * (Torque_calculation_inputs.joint_length_2 + link_2.length + Torque_calculation_inputs.radius_3),
                                       math.sin(alpha_rad) * (Torque_calculation_inputs.radius_1 + link_1.length + Torque_calculation_inputs.radius_2) + math.sin(beta_rad) * (Torque_calculation_inputs.joint_length_2 + link_2.length + Torque_calculation_inputs.radius_3),
                                       0]) + np.array([Torque_calculation_inputs.end_effector_distance, 0, 0])
                         ]
    border = 0.25
    xlimits = (np.min([x[0] for x in connection_points]) - border, np.max([x[0] for x in connection_points]) + border)
    ylimits = (np.min([x[1] for x in connection_points]) - border, np.max([x[1] for x in connection_points]) + border)
    axis.set_xlim(np.min([x[0] for x in connection_points]) - border, np.max([x[0] for x in connection_points]) + border)
    axis.set_ylim(np.min([x[1] for x in connection_points]) - border, np.max([x[1] for x in connection_points]) + border)
    # axis.axis('equal')
    # ax = axis.gca()
    # axis.gca()
    axis.set_aspect('equal', adjustable='box')
    # ax.set_aspect('equal', adjustable='box')
    # print(f"Limits {alpha_rad} {beta_rad}: {xlimits}; {ylimits}")
    axis.autoscale(False)
    for i in range(len(connection_points)-1):
        axis.plot((connection_points[i][0], connection_points[i+1][0]), (connection_points[i][1], connection_points[i+1][1]), 'ro-')
    axis.set_title(title)






static_torques = calculate_the_static_torque(alpha_rad, -beta_rad)
print(f"Holding Torque Joint 1: {static_torques[1]} Nm\nHolding Torque Joint 2: {static_torques[2]} Nm\nHolding Torque Joint 4: {static_torques[4]} Nm")
print()
print()

momenta_of_inertia = calculate_moment_of_inertia_tensors(alpha_rad, -beta_rad)
# print(f"Moment of inertia from joint 1:\n{momenta_of_inertia[1]}\n\nMoment of inertia from joint 2:\n\n{momenta_of_inertia[2]}\n\nMoment of inertia from joint 4:\n{momenta_of_inertia[4]}\n")

angular_accelerations = calculate_the_maximum_angular_acceleration(alpha_rad, -beta_rad, peak=False)
# print(f"The maximum acceleration with continuous motor torque in this position are:\nJoint 1: {angular_accelerations[1]}\nJoint 2: {angular_accelerations[2]}\nJoint 4: {angular_accelerations[4]}\n")

peak_angular_accelerations = calculate_the_maximum_angular_acceleration(alpha_rad, -beta_rad, peak=True)
# print(f"The maximum acceleration with peak motor torque in this position are:\nJoint 1: {peak_angular_accelerations[1]}\nJoint 2: {peak_angular_accelerations[2]}\nJoint 4: {peak_angular_accelerations[4]}\n")

forces = calculate_force_on_linkage_in_global_coordinate_system(alpha_rad, -beta_rad)
# print(f"Forces at the different cuts in the global coordinate system:\nJoint 1 acceleration:\n  C1: {forces[1]['C1']} N\n  C2: {forces[1]['C2']} N\n  C3: {forces[1]['C3']} N\n  C4: {forces[1]['C4']} N\n\n"
#       f"Joint 2 acceleration:\n  C1: {forces[2]['C1']} N\n  C2: {forces[2]['C2']} N\n  C3: {forces[2]['C3']} N\n  C4: {forces[2]['C4']} N\n\n"
#       f"Joint 4 acceleration:\n  C1: {forces[4]['C1']} N\n  C2: {forces[4]['C2']} N\n  C3: {forces[4]['C3']} N\n  C4: {forces[4]['C4']} N\n\n")

linkage_forces = calculate_force_on_linkage_in_linkage_coordinate_systems(alpha_rad, -beta_rad)
print(f"Forces at the different cuts in the respecitive linkage coordinate system:\nJoint 1 acceleration:\n  C1: {linkage_forces[1]['C1']} N\n  C2: {linkage_forces[1]['C2']} N\n  C3: {linkage_forces[1]['C3']} N\n  C4: {linkage_forces[1]['C4']} N\n\n"
      f"Joint 2 acceleration:\n  C1: {linkage_forces[2]['C1']} N\n  C2: {linkage_forces[2]['C2']} N\n  C3: {linkage_forces[2]['C3']} N\n  C4: {linkage_forces[2]['C4']} N\n\n"
      f"Joint 4 acceleration:\n  C1: {linkage_forces[4]['C1']} N\n  C2: {linkage_forces[4]['C2']} N\n  C3: {linkage_forces[4]['C3']} N\n  C4: {linkage_forces[4]['C4']} N\n\n")

print(f"Moments at the different cuts in the respecitive linkage coordinate system:\nJoint 1 acceleration:\n  C1-momentum: {linkage_forces[1]['C1-momentum']} Nm\n  C2-momentum: {linkage_forces[1]['C2-momentum']} Nm\n  C3-momentum: {linkage_forces[1]['C3-momentum']} Nm\n  C4-momentum: {linkage_forces[1]['C4-momentum']} Nm\n\n"
      f"Joint 2 acceleration:\n  C1-momentum: {linkage_forces[2]['C1-momentum']} Nm\n  C2-momentum: {linkage_forces[2]['C2-momentum']} Nm\n  C3-momentum: {linkage_forces[2]['C3-momentum']} Nm\n  C4-momentum: {linkage_forces[2]['C4-momentum']} Nm\n\n"
      f"Joint 4 acceleration:\n  C1-momentum: {linkage_forces[4]['C1-momentum']} Nm\n  C2-momentum: {linkage_forces[4]['C2-momentum']} Nm\n  C3-momentum: {linkage_forces[4]['C3-momentum']} Nm\n  C4-momentum: {linkage_forces[4]['C4-momentum']} Nm\n\n")


angular_velocity = np.array([0, 0, 2])
print(f"Energy for joint 1 with angular velocity: {angular_velocity} rad/s\n       {calculate_max_necessary_emergency_brake_energy(angular_velocity)} J\n\n")


state_input = f"This is the current state of the simulated robot:\n\n" \
                  f" link 1      - the mass: {link_1.mass} kg; the diameter: {Torque_calculation_inputs.tube_diameter_1} m; the tube thickness: {Torque_calculation_inputs.tube_wall_thickness_1} m;\n       the length: {link_1.length} m\n" \
                  f" link 2      - the mass: {link_2.mass} kg; the diameter: {Torque_calculation_inputs.tube_diameter_2} m; the tube thickness: {Torque_calculation_inputs.tube_wall_thickness_2} m;\n       the length: {link_2.length} m\n\n" \
                  f"The payload  - mass: {Torque_calculation_inputs.payload} kg\n\n" \
                  f"The Countwerweight - mass: {Torque_calculation_inputs.counterweight} kg,  the distance: {Torque_calculation_inputs.counterweight_distance} m\n\n" \
                  f"The length radius 1: {Torque_calculation_inputs.radius_1} m;  The length radius 2: {Torque_calculation_inputs.radius_2} m;  The length radius 3: {Torque_calculation_inputs.radius_3}\n\n" \
                  f"Joint 23     - mass: {Torque_calculation_inputs.mass_joint_2_3} kg;  length: {Torque_calculation_inputs.joint_length_2} m\n" \
                  f"Joint 45     - mass: {Torque_calculation_inputs.mass_joint_4_5} kg\n" \
                  f"End-effector - Mass is: {Torque_calculation_inputs.end_effector_mass} kg;  the length of the end-effector is: {Torque_calculation_inputs.end_effector_distance} m\n\n" \
                  f"The holding torques:\n" \
                  f"At angle alpha: {Torque_calculation_inputs.alpha}°; beta: {Torque_calculation_inputs.beta}°\n\n" \
                  f"Holding Torque Joint 1: {static_torques[1]} Nm\nHolding Torque Joint 2: {static_torques[2]} Nm\nHolding Torque Joint 4: {static_torques[4]} Nm\n\n" \
                  f"Full stretched length: {Torque_calculation_inputs.full_length}"


positions = calculate_the_position_of_the_center_of_mass(angle_alpha=math.radians(Torque_calculation_inputs.alpha),angle_beta=math.radians(-Torque_calculation_inputs.beta))
print(positions[-2])

save_data = False

number = 7
filename = f"Static_load_case_{number}"
if(save_data):
    with open(f"data/{filename}.txt", 'a') as f:
        f.write(state_input)
print(state_input)

# print(f"alpha: {alpha_rad/math.pi * 180}; beta: {beta_rad/math.pi * 180}")
# create_display_depending_on_angles(alpha_rad, -beta_rad)
