import math
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import axes3d
import Torque_calculation
import numpy as np
import Torque_calculation_inputs
import datetime
import pickle

filename_input_infos = input("""additional part for filename:\n"YYYY-MM-DD_simulation-data-": """)

alphas = []
betas = []
forces = []
# torques = []

number_of_alphas = 20 * 5
number_of_betas = 18 * 10

range_of_alphas = 100
range_of_betas = 180

for i in range(number_of_alphas):
    alpha = range_of_alphas*i/number_of_alphas
    print(f"\n{i+1}", end='')
    for e in range(number_of_betas):
        print('.', end='')
        beta = (-number_of_betas/2+e)/number_of_betas * range_of_betas

        # this condition exists so that the two linkages to not "overlap"
        if(alpha + beta < (180)):
            alpha_rad = alpha/180 * math.pi
            beta_rad = beta/180 * math.pi
            force = Torque_calculation.calculate_force_on_linkage_in_linkage_coordinate_systems(angle_alpha=alpha_rad, angle_beta=-beta_rad)
            # torque = Torque_calculation.calculate_force_on_linkage_in_global_coordinate_system(alpha_rad, -beta_rad)

            # torques.append(torque)
            alphas.append(alpha)
            betas.append(beta)
            forces.append(force)
print()
print()
print()
forces_better_organised = [[[i[1]["C1"] for i in forces], [i[1]["C2"] for i in forces], [i[1]["C3"] for i in forces], [i[1]["C4"] for i in forces]],
                           [[i[2]["C1"] for i in forces], [i[2]["C2"] for i in forces], [i[2]["C3"] for i in forces], [i[2]["C4"] for i in forces]],
                           [[i[4]["C1"] for i in forces], [i[4]["C2"] for i in forces], [i[4]["C3"] for i in forces], [i[4]["C4"] for i in forces]]]

# print(forces_better_organised)

joints = [1, 2, 4]

"""
for joint in range(3):
    for c_point in range(4):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        X = np.array(alphas)
        Y = np.array(betas)

        Z1 = np.array([i[0] for i in forces_better_organised[joint][c_point]])
        Z2 = np.array([i[1] for i in forces_better_organised[joint][c_point]])
        Z3 = np.array([i[2] for i in forces_better_organised[joint][c_point]])
        # print(len(Z1))
        ax.plot3D(X, Y, Z1, label='1')
        # ax.plot_surface(X, Y, Z1, edgecolor='royalblue', lw=0.5, rstride=8, cstride=8,
        #         alpha=0.3)
        ax.plot3D(X, Y, Z2, label='2')
        ax.plot3D(X, Y, Z3, label='3')
        plt.legend()
        # print(f"Joint {joints[joint]}, cutting point: {c_point + 1}")
        # plt.show()

"""


forces_organised_for_argmax = [[i[1]["C1"] for i in forces] + [i[2]["C1"] for i in forces] + [i[4]["C1"] for i in forces],
                                [i[1]["C2"] for i in forces] + [i[2]["C2"] for i in forces] + [i[4]["C2"] for i in forces],
                                [i[1]["C3"] for i in forces] + [i[2]["C3"] for i in forces] + [i[4]["C3"] for i in forces],
                                [i[1]["C4"] for i in forces] + [i[2]["C4"] for i in forces] + [i[4]["C4"] for i in forces]]

torques_organised_for_argmax = [[i[1]["C1-momentum"] for i in forces] + [i[2]["C1-momentum"] for i in forces] + [i[4]["C1-momentum"] for i in forces],
                                [i[1]["C2-momentum"] for i in forces] + [i[2]["C2-momentum"] for i in forces] + [i[4]["C2-momentum"] for i in forces],
                                [i[1]["C3-momentum"] for i in forces] + [i[2]["C3-momentum"] for i in forces] + [i[4]["C3-momentum"] for i in forces],
                                [i[1]["C4-momentum"] for i in forces] + [i[2]["C4-momentum"] for i in forces] + [i[4]["C4-momentum"] for i in forces]]

compression_forces_organised_for_argmax = [-forces[i][1]["C1"][0] + forces[i][1]["C1-momentum"][2]/2 * 12.5 for i in range(len(forces))] \
                                           + [-forces[i][2]["C1"][0] + forces[i][2]["C1-momentum"][2]/2 * 12.5 for i in range(len(forces))] \
                                           + [-forces[i][4]["C1"][0] + forces[i][4]["C1-momentum"][2]/2 * 12.5 for i in range(len(forces))]

# print(forces_organised_for_argmax)
# print(torques_organised_for_argmax)

date = datetime.date.today()
filename = f"{date.isoformat()}_simulation-data-{filename_input_infos}.txt"

with open(f"data/{filename}", 'a') as f:
    # f.write('Hello\n')
    state_input = f"This is the current state of the simulated robot:\n\n" \
                  f" link 1      - the mass: {Torque_calculation.link_1.mass} kg;  the length: {Torque_calculation.link_1.length} m\n" \
                  f" link 2      - the mass: {Torque_calculation.link_2.mass} kg;  the length: {Torque_calculation.link_2.length} m\n\n" \
                  f"The payload  - mass: {Torque_calculation_inputs.payload} kg\n\n" \
                  f"The Countwerweight - mass: {Torque_calculation_inputs.counterweight} kg,  the distance: {Torque_calculation_inputs.counterweight_distance} m\n\n" \
                  f"The length radius 1: {Torque_calculation_inputs.radius_1} m;  The length radius 2: {Torque_calculation_inputs.radius_2} m;  The length radius 3: {Torque_calculation_inputs.radius_3}\n\n" \
                  f"Joint 23     - mass: {Torque_calculation_inputs.mass_joint_2_3} kg;  length: {Torque_calculation_inputs.joint_length_2} m\n" \
                  f"Joint 45     - mass: {Torque_calculation_inputs.mass_joint_4_5} kg\n" \
                  f"End-effector - Mass is: {Torque_calculation_inputs.end_effector_mass} kg;  the length of the end-effector is: {0.3} m\n\n" \
                  f"The motor torques:\n" \
                  f"  Motor 0     - continuous: {Torque_calculation.joint_0.continuous} Nm,  peak: {Torque_calculation.joint_0.peak} Nm\n" \
                  f"  Motor 1     - continuous: {Torque_calculation.joint_1.continuous} Nm,  peak: {Torque_calculation.joint_1.peak} Nm\n" \
                  f"  Motor 2     - continuous: {Torque_calculation.joint_2.continuous} Nm,  peak: {Torque_calculation.joint_2.peak} Nm\n" \
                  f"  Motor 3     - continuous: {Torque_calculation.joint_3.continuous} Nm,  peak: {Torque_calculation.joint_3.peak} Nm\n" \
                  f"  Motor 4     - continuous: {Torque_calculation.joint_4.continuous} Nm,  peak: {Torque_calculation.joint_4.peak} Nm\n" \
                  f"  Motor 5     - continuous: {Torque_calculation.joint_5.continuous} Nm,  peak: {Torque_calculation.joint_5.peak} Nm\n\n" \
                  f"The maximum angular accelerations:\n" \
                  f"  The maximum accelerations - Joint 0: {Torque_calculation_inputs.Max_angular_acceleration_joint_0} rad/s^2\n" \
                  f"  The maximum accelerations - Joint 1: {Torque_calculation_inputs.Max_angular_acceleration_joint_1} rad/s^2\n" \
                  f"  The maximum accelerations - Joint 2: {Torque_calculation_inputs.Max_angular_acceleration_joint_2} rad/s^2\n" \
                  f"  The maximum accelerations - Joint 3: {Torque_calculation_inputs.Max_angular_acceleration_joint_3} rad/s^2\n" \
                  f"  The maximum accelerations - Joint 4: {Torque_calculation_inputs.Max_angular_acceleration_joint_4} rad/s^2\n" \
                  f"  The maximum accelerations - Joint 5: {Torque_calculation_inputs.Max_angular_acceleration_joint_5} rad/s^2\n" \
                  f"\n\n\n------------------------------------------------------------------------------------------\n\n"

    f.write(state_input)

    for i in range(4):
        # figure, axis = plt.subplots(4, 2)
        figure = plt.figure()

        plot1 = plt.subplot2grid((4, 4), (1, 0))

        plot2 = plt.subplot2grid((4, 4), (1, 1))
        plot3 = plt.subplot2grid((4, 4), (1, 2))
        plot4 = plt.subplot2grid((4, 4), (1, 3))
        plot5 = plt.subplot2grid((4, 4), (2, 2))
        plot6 = plt.subplot2grid((4, 4), (3, 1))
        plot7 = plt.subplot2grid((4, 4), (3, 2))

        plt.tight_layout()


        print(f"Cut {i+1}")
        plt.suptitle(f"Cut {i+1}")

        f.write(f"Cut {i+1}\n")
        max_x_force_pos = np.argmax(np.array([i[0] for i in forces_organised_for_argmax[i]]))
        min_x_force_pos = np.argmin(np.array([i[0] for i in forces_organised_for_argmax[i]]))
        infos = f"Max x-force: {forces_organised_for_argmax[i][max_x_force_pos]} N\n" \
                f"The corresponding torque: {torques_organised_for_argmax[i][max_x_force_pos]} Nm\n" \
                f"acceleration of the joint: {(max_x_force_pos//len(forces)) + 1}\n" \
                f"at angle: alpha = {alphas[max_x_force_pos % len(forces)]}°, beta = {betas[max_x_force_pos % len(forces)]}°\n\n"
        f.write(infos)
        print(infos)
        Torque_calculation.create_display_depending_on_angles(alphas[max_x_force_pos % len(forces)] * math.pi / 180, -betas[max_x_force_pos % len(forces)] * math.pi / 180, plot1, "Max x-force")

        infos = f"Min x-force: {forces_organised_for_argmax[i][min_x_force_pos]} N\n" \
                f"the corresponding torque: {torques_organised_for_argmax[i][min_x_force_pos]} Nm\n" \
                f"acceleration of the joint: {(min_x_force_pos//len(forces)) + 1}\n" \
                f"at angle: alpha = {alphas[min_x_force_pos % len(forces)]}°, beta = {betas[min_x_force_pos % len(forces)]}°\n\n"
        f.write(infos)
        print(infos)
        Torque_calculation.create_display_depending_on_angles(alphas[min_x_force_pos % len(forces)] * math.pi / 180, -betas[min_x_force_pos % len(forces)] * math.pi / 180, plot2, "Min x-force")

        max_y_force_pos = np.argmax(np.array([i[1] for i in forces_organised_for_argmax[i]]))
        min_y_force_pos = np.argmin(np.array([i[1] for i in forces_organised_for_argmax[i]]))
        infos = f"Max y-force: {forces_organised_for_argmax[i][max_y_force_pos]} N\n" \
                f"the corresponding torque: {torques_organised_for_argmax[i][max_y_force_pos]} Nm\n" \
                f"acceleration of the joint: {(max_y_force_pos//len(forces)) + 1}\n" \
                f"at angle: alpha = {alphas[max_y_force_pos % len(forces)]}°, beta = {betas[max_y_force_pos % len(forces)]}°\n\n"
        f.write(infos)
        print(infos)
        Torque_calculation.create_display_depending_on_angles(alphas[max_y_force_pos % len(forces)] * math.pi / 180, -betas[max_y_force_pos % len(forces)] * math.pi / 180, plot3, "Max y-force")

        infos = f"Min y-force: {forces_organised_for_argmax[i][min_y_force_pos]} N\n" \
                f"the corresponding torque: {torques_organised_for_argmax[i][min_y_force_pos]} Nm\n" \
                f"acceleration of the joint: {(min_y_force_pos//len(forces)) + 1}\n" \
                f"at angle: alpha = {alphas[min_y_force_pos % len(forces)]}°, beta = {betas[min_y_force_pos % len(forces)]}°\n\n"
        f.write(infos)
        print(infos)
        Torque_calculation.create_display_depending_on_angles(alphas[min_y_force_pos % len(forces)] * math.pi / 180, -betas[min_y_force_pos % len(forces)] * math.pi / 180, plot4, "Min y-force")

        max_total_force_pos = np.argmax(np.array([np.linalg.norm(i) for i in forces_organised_for_argmax[i]]))
        infos = f"Max total-force: {forces_organised_for_argmax[i][max_total_force_pos]} N\n" \
                f"the corresponding torque: {torques_organised_for_argmax[i][max_total_force_pos]} Nm\n" \
                f"acceleration of the joint: {(max_total_force_pos//len(forces)) + 1}\n" \
                f"at angle: alpha = {alphas[max_total_force_pos % len(forces)]}°, beta = {betas[max_total_force_pos % len(forces)]}°\n\n"
        f.write(infos)
        print(infos)
        Torque_calculation.create_display_depending_on_angles(alphas[max_total_force_pos % len(forces)] * math.pi / 180, -betas[max_total_force_pos % len(forces)] * math.pi / 180, plot5, "Max total force")

        max_z_torque_pos = np.argmax(np.array([i[2] for i in torques_organised_for_argmax[i]]))
        min_z_torque_pos = np.argmin(np.array([i[2] for i in torques_organised_for_argmax[i]]))
        infos = f"Max z-torque: {torques_organised_for_argmax[i][max_z_torque_pos]} Nm\n" \
                f"the corresponding force: {forces_organised_for_argmax[i][max_z_torque_pos]} N\n" \
                f"acceleration of the joint: {(max_z_torque_pos//len(forces)) + 1}\n" \
                f"at angle: alpha = {alphas[max_z_torque_pos % len(forces)]}°, beta = {betas[max_z_torque_pos % len(forces)]}°\n\n"
        f.write(infos)
        print(infos)
        Torque_calculation.create_display_depending_on_angles(alphas[max_z_torque_pos % len(forces)] * math.pi / 180, -betas[max_z_torque_pos % len(forces)] * math.pi / 180, plot6, "Max z-torque")

        infos = f"Min z-torque: {torques_organised_for_argmax[i][min_z_torque_pos]} Nm\n" \
                f"the corresponding force: {forces_organised_for_argmax[i][min_z_torque_pos]} N\n" \
                f"acceleration of the joint: {(min_z_torque_pos//len(forces)) + 1}\n" \
                f"at angle: alpha = {alphas[min_z_torque_pos % len(forces)]}°, beta = {betas[min_z_torque_pos % len(forces)]}°\n\n"
        f.write(infos)
        print(infos)
        Torque_calculation.create_display_depending_on_angles(alphas[min_z_torque_pos % len(forces)] * math.pi / 180, -betas[min_z_torque_pos % len(forces)] * math.pi / 180, plot7, "Min z-torque")
        f.write("\n\n-----------------------------\n\n")
        # plt.show()
        image_name = f"data/{date.isoformat()}_simulation-data-{filename_input_infos}-Cut-{i+1}.png"
        plt.savefig(image_name)
        print()
        print()
        print("-------")
        print()



    print("\n\n\nCut 1:")
    f.write("Cut 1:\n")
    max_compression_force_pos = np.argmax(np.array(compression_forces_organised_for_argmax))
    infos = f"Max compression force on a point 10 cm away from the center: {compression_forces_organised_for_argmax[max_compression_force_pos]} N\n" \
            f"the corresponding force: {forces_organised_for_argmax[0][max_compression_force_pos]} N\n" \
            f"The corresponding torque: {torques_organised_for_argmax[0][max_compression_force_pos]} Nm\n" \
            f"acceleration of the joint: {(max_compression_force_pos//len(forces)) + 1}\n" \
            f"at angle: alpha = {alphas[max_compression_force_pos % len(forces)]}°, beta = {betas[max_compression_force_pos % len(forces)]}°\n\n"
    print(infos)
    f.write(infos)
    f.write("----------------------------------------------------------------------------------------\n\n")
    print()
    print()



"""fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# Grab some test data.
X, Y, Z = axes3d.get_test_data(0.05)
print(X)
print(Z)
# Plot a basic wireframe.
ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)

plt.show()"""


alpha = Torque_calculation_inputs.alpha
beta = Torque_calculation_inputs.beta





alpha_rad = alpha/180 * math.pi
beta_rad = -beta/180 * math.pi

