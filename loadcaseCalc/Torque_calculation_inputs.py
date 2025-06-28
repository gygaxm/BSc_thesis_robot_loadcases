import math


length = 0.3

length_link_1 = length
length_link_2 = length

# number_of_segments_in_link_1 = 4
# number_of_segments_in_link_2 = 5




# mass_link_1 = 10
# mass_link_2 = 5.5

# angles in °
alpha = 98
beta = 70


# payload in kg
payload = 12

# counterweight in kg
counterweight = 0



# counterweight distance
counterweight_distance = 0.3

# Constants
# length_per_segment_link_1 = 0.14
# length_per_segment_link_2 = 0.12

radius_1 = 0.1
radius_2 = 0.1
joint_length_2 = 0.1
radius_3 = 0.1

# mass of the joint connection which encompasses the joint 2 & 3 in [kg]
mass_joint_2_3 = 16

# mass of the joint connection which encompasses the joint 2 & 3 in [kg]
mass_joint_4_5 = 8

# the safety factor for the holding torque without overheating the motors
safety_factor = 0.9


tube_diameter_1 = 0.18
tube_wall_thickness_1 = 0.005

tube_diameter_2 = 0.15
tube_wall_thickness_2 = 0.004



end_effector_distance = 0.2
end_effector_mass = 1.5



full_length = length_link_1 + length_link_2 + radius_1 + radius_2 + radius_3 + joint_length_2 + end_effector_distance


Max_angular_acceleration_joint_0 = 5 * 100
Max_angular_acceleration_joint_1 = 5 * 100
Max_angular_acceleration_joint_2 = 6 * 100
Max_angular_acceleration_joint_3 = 7 * 100
Max_angular_acceleration_joint_4 = 8 * 100
Max_angular_acceleration_joint_5 = 8 * 100


# the continuous and peak torques
Motor_0 = {'continuous': 140, 'peak': 400}
Motor_1 = {'continuous': 350, 'peak': 800}
Motor_2 = {'continuous': 140, 'peak': 400}
Motor_3 = {'continuous': 40, 'peak': 90}
Motor_4 = {'continuous': 40, 'peak': 90}
Motor_5 = {'continuous': 8, 'peak': 15}


print(f"Total Length: {length_link_1 + length_link_2 + radius_1 + radius_2 + radius_3 + joint_length_2 + end_effector_distance}\n\n")
