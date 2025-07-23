import roboticstoolbox as rtb
import numpy as np
import itertools
import pandas as pd
import datetime
import math
import matplotlib.pyplot as plt
from tqdm import tqdm
from spatialmath.base.argcheck import getvector, isscalar, verifymatrix, getmatrix # Argument checking
from spatialmath.base import tr2jac, tr2eul, tr2rpy, t2r, trlog, rotvelxform # Transformations
import spatialmath.base.symbolic as sym # For symbolic models (used internally, keep for full compatibility)




def _cross(a, b):
    return np.r_[
        a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]
    ]

# MODIFIED rne_python function (copied from RoboticToolbox source)
def python_rne_with_wrenches(
    robot_obj,
    Q,
    QD=None,
    QDD=None,
    gravity=None,
    fext=None,
    debug=False,
    base_wrench=False,
):
    """
    Compute inverse dynamics via recursive Newton-Euler formulation.
    Modified to return internal link wrenches (forces and moments) as well.
    """
    # Alias for consistency with original code's 'self'
    self = robot_obj

    # Handle base transformation
    if np.array_equal(self.base.A, np.eye(4)):
        base = None
    else:
        base = self.base.A

    # Helper for debug printing (from original source)
    def removesmall(x):
        return x

    n = self.n # Number of robot joints

    if self.symbolic:
        dtype = "O" # Object type for symbolic variables
    else:
        dtype = None # Default numerical type

    z0 = np.array([0, 0, 1], dtype=dtype) # Z-axis unit vector

    if gravity is None:
        gravity = self.gravity  # Use default gravity from the robot object
    else:
        gravity = getvector(gravity, 3) # Ensure gravity is a 3-element vector

    if fext is None:
        fext = np.zeros((6,), dtype=dtype) # Default to zero wrench
    else:
        fext = getvector(fext, 6) # Ensure fext is a 6-element vector

    if debug:
        print("grav", gravity)
        print("fext", fext)

    # Unpack joint coordinates and derivatives (handling single configuration or trajectory inputs)
    if Q is not None and QD is None and QDD is None:
        # Case: single argument P = [q, qd, qdd]
        Q = getmatrix(Q, (None, self.n * 3))
        q = Q[:, 0:n]
        qd = Q[:, n : 2 * n]
        qdd = Q[:, 2 * n :]
    else:
        # Case: three arguments q, qd, qdd
        q = getmatrix(Q, (None, self.n))
        qd = getmatrix(QD, (None, self.n))
        qdd = getmatrix(QDD, (None, self.n))

    nk = q.shape[0] # Number of trajectory points (usually 1 for a single configuration)

    tau = np.zeros((nk, n), dtype=dtype) # Array to store joint torques/forces output

    if base_wrench:
        wbase = np.zeros((nk, 6), dtype=dtype) # Dimension changed to 6 for full wrench

    # START OF MODIFICATION
    # Initialize a list to store internal wrenches
    # Each element will be a list of wrenches (force and moment vectors) for each joint in the order of the backward recursion (from end-effector towards base).
    internal_wrenches_per_timestep = [None] * nk
    # END OF MODIFICATION

    for k in range(nk):
        q_k = q[k, :]
        qd_k = qd[k, :]
        qdd_k = qdd[k, :]

        if debug:
            print("q_k", q_k)
            print("qd_k", qd_k)
            print("qdd_k", qdd_k)
            print()

        Fm = np.zeros((3, n), dtype=dtype)
        Nm = np.zeros((3, n), dtype=dtype)
        pstarm = np.zeros((3, n), dtype=dtype)
        Rm = []

        w = np.zeros((3,), dtype=dtype)
        wd = np.zeros((3,), dtype=dtype)
        vd = -gravity

        if base is not None:
            Rb = t2r(base).T
            w = Rb @ w
            wd = Rb @ wd
            vd = Rb @ gravity

        for j in range(n):
            link = self.links[j]

            if link.sigma == 0:
                Tj = link.A(q_k[j]).A
                d = link.d
            else:
                Tj = link.A(link.theta).A
                d = q_k[j]

            alpha = link.alpha
            if self.mdh:
                pstar = np.r_[link.a, -d * sym.sin(alpha), d * sym.cos(alpha)]
                if j == 0:
                    if base:
                        Tj = base @ Tj
                        # START OF MODIFICATION
                        pstar = t2r(base) @ pstar
                        # END OF MODIFICATION
            else:
                pstar = np.r_[link.a, d * sym.sin(alpha), d * sym.cos(alpha)]

            Rm.append(t2r(Tj))
            pstarm[:, j] = pstar

        for j, link in enumerate(self.links):
            Rt = Rm[j].T
            pstar = pstarm[:, j]
            r = link.r

            if self.mdh:
                if link.isrevolute:
                    w_ = Rt @ w + z0 * qd_k[j]
                    wd_ = Rt @ wd + z0 * qdd_k[j] + _cross(Rt @ w, z0 * qd_k[j])
                    vd_ = Rt @ _cross(wd, pstar) + _cross(w, _cross(w, pstar)) + vd
                else:
                    w_ = Rt @ w
                    wd_ = Rt @ wd
                    vd_ = (
                        Rt @ (_cross(wd, pstar) + _cross(w, _cross(w, pstar)) + vd)
                        + 2 * _cross(Rt @ w, z0 * qd_k[j])
                        + z0 * qdd_k[j]
                    )
                w = w_
                wd = wd_
                vd = vd_
            else:
                if link.isrevolute:
                    wd = Rt @ (wd + z0 * qdd_k[j] + _cross(w, z0 * qd_k[j]))
                    w = Rt @ (w + z0 * qd_k[j])
                    vd = _cross(wd, pstar) + _cross(w, _cross(w, pstar)) + Rt @ vd
                else:
                    w = Rt @ w
                    wd = Rt @ wd
                    vd = (
                        Rt @ (z0 * qdd_k[j] + vd)
                        + _cross(wd, pstar)
                        + 2 * _cross(w, Rt @ z0 * qd_k[j])
                        + _cross(w, _cross(w, pstar))
                    )

            vhat = _cross(wd, r) + _cross(w, _cross(w, r)) + vd
            Fm[:, j] = link.m * vhat
            Nm[:, j] = link.I @ wd + _cross(w, link.I @ w)

            if debug:
                print("w:     ", removesmall(w))
                print("wd:    ", removesmall(wd))
                print("vd:    ", removesmall(vd))
                print("vdbar: ", removesmall(vhat))
                print()

        if debug:
            print("Fm\n", Fm)
            print("Nm\n", Nm)

        f = fext[:3]
        nn = fext[3:]

        # TART OF MODIFICATION
        current_k_wrenches_temp = [None] * n
        # END OF MODIFICATION

        for j in range(n - 1, -1, -1):
            link = self.links[j]
            r = link.r

            if self.mdh:
                if j == (n - 1):
                    R = np.eye(3, dtype=dtype)
                    pstar = np.zeros((3,), dtype=dtype)
                else:
                    R = Rm[j + 1]
                    pstar = pstarm[:, j + 1]

                f_ = R @ f + Fm[:, j]
                nn_ = (
                    R @ nn
                    + _cross(pstar, R @ f)
                    + _cross(pstar, Fm[:, j])
                    + Nm[:, j]
                )
                f = f_
                nn = nn_

            else:
                pstar = pstarm[:, j]
                if j == (n - 1):
                    R = np.eye(3, dtype=dtype)
                else:
                    R = Rm[j + 1]

                nn = (
                    R @ (nn + _cross(R.T @ pstar, f))
                    + _cross(pstar + r, Fm[:, j])
                    + Nm[:, j]
                )
                f = R @ f + Fm[:, j]

            if debug:
                print("f: ", removesmall(f))
                print("n: ", removesmall(nn))

            # START OF MODIFICATION
            # Store the internal wrench (combined force 'f' and moment 'nn')
            # This wrench represents the force/moment exerted BY link (j-1) ON link (j),
            # expressed in the coordinate frame of link (j).
            # The format is [fx, fy, fz, nx, ny, nz].
            current_k_wrenches_temp[j] = np.r_[f, nn]
            # END OF MODIFICATION

            R = Rm[j]
            if self.mdh:
                if link.isrevolute:
                    t = nn @ z0
                else:
                    t = f @ z0
            else:
                if link.isrevolute:
                    t = nn @ (R.T @ z0)
                else:
                    t = f @ (R.T @ z0)

            tau[k, j] = (
                t
                + link.G**2 * link.Jm * qdd_k[j]
                - link.friction(qd_k[j], coulomb=not self.symbolic)
            )
            if debug:
                print(
                    f"j={j:}, G={link.G:}, Jm={link.Jm:},"
                    f" friction={link.friction(qd_k[j], coulomb=False):}"
                )
                print()

        # START OF MODIFICATION
        internal_wrenches_per_timestep[k] = list(current_k_wrenches_temp)
        # END OF MODIFICATION

        if base_wrench:
            R = Rm[0]
            nn = R @ nn
            f = R @ f
            wbase[k, :] = np.r_[f, nn]

    # START OF MODIFICATION
    internal_wrenches_np = np.array(internal_wrenches_per_timestep)
    # END OF MODIFICATION

    if base_wrench:
        if tau.shape[0] == 1:
            # START OF MODIFICATION
            return tau.flatten(), wbase.flatten(), internal_wrenches_np.squeeze()
            # END OF MODIFICATION

        else:
            # START OF MODIFICATION
            # Return tau, base wrench, and internal wrenches.
            return tau, wbase, internal_wrenches_np
            # END OF MODIFICATION

    else:
        if tau.shape[0] == 1:
            # START OF MODIFICATION
            # Return tau (flattened) and internal wrenches (squeezed for single timestep).
            return tau.flatten(), internal_wrenches_np.squeeze()
            # END OF MODIFICATION

        else:
            # START OF MODIFICATION
            # Return tau and internal wrenches.
            return tau, internal_wrenches_np
            # END OF MODIFICATION



# obot Kinematic Model ad Parameter Definition
deg = np.deg2rad

L0 = rtb.RevoluteDH(d=0.324, a=0, alpha=deg(-90), offset=deg(-90), qlim=[-np.pi/2, np.pi/2])
L1 = rtb.RevoluteDH(d=0, a=0.493, alpha=deg(0), offset=deg(90), qlim=[-np.pi/3, np.pi/2])   #minus 20cm for shortest link (we want 35cm optimized links) original value: a=0.693
L2 = rtb.RevoluteDH(d=0.165, a=0, alpha=deg(90), offset=deg(-90), qlim=[-np.pi, np.pi])
L3 = rtb.RevoluteDH(d=0.465, a=0, alpha=deg(90), offset=deg(90), qlim=[-np.pi, np.pi])      #minus 20cm for shortest link (we want 35cm optimized links) original value: d=0.665
L4 = rtb.RevoluteDH(d=0.1271, a=0, alpha=deg(-90), offset=deg(0), qlim=[-np.pi, np.pi])
L5 = rtb.RevoluteDH(d=0.255, a=0, alpha=deg(-90), offset=deg(0), qlim=[-np.pi, np.pi])

robot = rtb.DHRobot([L0, L1, L2, L3, L4, L5], name="ASIRA")

q_plot = [0, 0, 0, 0, 0, 0]
max_reach = robot.reach + 0.5
plot_limits = [-max_reach, max_reach, -max_reach, max_reach, -0.1, max_reach]
print(f"Robot estimated reach: {robot.reach:.2f} meters")
#print(f"Plot limits set to: {plot_limits}")
#robot.plot(q_plot, limits=plot_limits)
#plt.show()

# -rbt proprts
robot.links[0].m = 140
robot.links[0].r = np.array([0, 0, L0.d / 2])
robot.links[0].I = np.array([
                            [4.4174e+06, -15.3919698, 3806.3707],
                            [-15.3919698, 2.49667e+06, -3155.6593],
                            [3806.3707, -3155.6593, 3.38967e+06]
                            ]) * 1e-6 # kg mm^2 to kg m^2

robot.links[1].m = 18
robot.links[1].r = np.array([0, -(L1.a/3), 0])
robot.links[1].I = np.array([
                            [1.18959e+06, 8224.5820, -187.7660],
                            [8224.5820, 1.21665e+06, 27513.2121],
                            [-187.7660, 27513.2121, 1.17085e+06]
                            ]) * 1e-6 # kg mm^2 to kg m^2

robot.links[2].m = 15
robot.links[2].r = np.array([0, 0, L2.d/3])
robot.links[2].I = np.array([
                            [131003, -898.6573, 40256.7832],
                            [-898.6573, 164123, -1094.7972],
                            [40256.7832, -1094.7972, 67745.7324]
                            ]) * 1e-6 # kg mm^2 to kg m^2

robot.links[3].m = 6.5
robot.links[3].r = np.array([0, 0, L3.d / 2])
robot.links[3].I = np.array([
                            [290592, 9.9085, -35.3614],
                            [9.9085, 289672, -3282.0726],
                            [-35.3614, -3282.0726, 32998.5850]
                            ]) * 1e-6 # kg mm^2 to kg m^2

robot.links[4].m = 5
robot.links[4].r = np.array([0, 0, L4.d/3])
robot.links[4].I = np.array([
                            [40734.1776, -288.5800, 748.6597],
                            [-288.5800, 31112.2640, -9046.1309],
                            [748.6597, -9046.1309, 15143.5161]
                            ]) * 1e-6 # kg mm^2 to kg m^2

payload_mass = 4.0
end_effector_mass = 1.0
end_effector_length = 0.2
robot.links[5].m = 1.0 + payload_mass + end_effector_mass
robot.links[5].r = np.array([0, 0, L5.d/2])
robot.links[5].I = np.array([
                            [20579.0323, 12.7478, -0.9320],
                            [12.7478, 4715.3192, 6.9208],
                            [-0.9320, 6.9208, 17475.0969]
                            ]) * 1e-6 # kg mm^2 to kg m^2

robot.gravity = np.array([0, 0, 0])


# sweep parameters
num_q_steps = 1
num_qd_steps = 2
num_qdd_steps = 2

total_overall_combinations = num_q_steps**6 * num_qd_steps**6 * num_qdd_steps**6
print("Calculating", total_overall_combinations ,"loadcases.")

q_ranges = []
qd_ranges = []
qdd_ranges = []

#torques:
# DT85L --> 480nm
# DT50M --> 300nm
# DT38M --> 200nm


#speeds:
# DT85L --> 3.1 rad/s (Joint 0/1/2)
# DT50M --> 4.66 rad/s (Joint 3/4)
# DT38M --> 5.2 rad/s (Joint 5)

#accelerations:
# J0 = 5 rad/s**2
# J1 = 5 rad/s**2
# J2 = 6 rad/s**2
# J3 = 7 rad/s**2
# J4 = 8 rad/s**2
# J5 = 8 rad/s**2


real_qd_max_values = np.array([3.1, 3.1, 3.1, 4.66, 4.66, 5.2])
real_qdd_max_values = np.array([5.0, 5.0, 6.0, 7.0, 8.0, 8.0])


for i, link in enumerate(robot.links):
    q_ranges.append(np.linspace(link.qlim[0], link.qlim[1], num_q_steps))
    qd_ranges.append(np.linspace(-real_qd_max_values[i], real_qd_max_values[i], num_qd_steps))
    qdd_ranges.append(np.linspace(-real_qdd_max_values[i], real_qdd_max_values[i], num_qdd_steps))


all_q_combinations_outer = itertools.product(*q_ranges)

#min max tracker init
max_min_tracker = {}
component_names_list = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']

initial_q_state = [0.0] * robot.n
initial_qd_state = [0.0] * robot.n
initial_qdd_state = [0.0] * robot.n
initial_wrench_data = [0.0] * 6

for j in range(robot.n):
    for comp_name in component_names_list: # Fx, Fy, Fz, Mx, My, Mz
        max_min_tracker[(comp_name, j, 'max')] = [-float('inf'), initial_wrench_data, initial_q_state, initial_qd_state, initial_qdd_state]
        max_min_tracker[(comp_name, j, 'min')] = [float('inf'), initial_wrench_data, initial_q_state, initial_qd_state, initial_qdd_state]
    max_min_tracker[('F_resultant', j, 'max')] = [-float('inf'), initial_wrench_data, initial_q_state, initial_qd_state, initial_qdd_state]
    max_min_tracker[('M_resultant', j, 'max')] = [-float('inf'), initial_wrench_data, initial_q_state, initial_qd_state, initial_qdd_state]

current_combination_count = 0

# sim loop
# This loop itrates through all possible combnations of joint position, velocities, and accelerations.
with tqdm(total=total_overall_combinations, desc="Simulating Combinations") as pbar:
    for q_tuple in all_q_combinations_outer:
        all_qd_combinations_inner = itertools.product(*qd_ranges)

        for qd_tuple in all_qd_combinations_inner:
            all_qdd_combinations_innermost = itertools.product(*qdd_ranges)

            for qdd_tuple in all_qdd_combinations_innermost:
                q = np.array(q_tuple)
                qd = np.array(qd_tuple)
                qdd = np.array(qdd_tuple)

                pbar.update(1)

                try:
                    tau, internal_wrenches_np = python_rne_with_wrenches(
                        robot, q, qd, qdd, gravity=robot.gravity
                    )

                    #update maxmins
                    for j in range(robot.n):
                        wrench_data_for_this_joint = internal_wrenches_np[j]

                        for comp_idx, comp_value in enumerate(wrench_data_for_this_joint):
                            comp_name = component_names_list[comp_idx]

                            max_key = (comp_name, j, 'max')
                            if comp_value > max_min_tracker[max_key][0]:
                                max_min_tracker[max_key] = [comp_value, wrench_data_for_this_joint.tolist(), q.tolist(), qd.tolist(), qdd.tolist()]

                            min_key = (comp_name, j, 'min')
                            if comp_value < max_min_tracker[min_key][0]:
                                max_min_tracker[min_key] = [comp_value, wrench_data_for_this_joint.tolist(), q.tolist(), qd.tolist(), qdd.tolist()]

                        # Update Max Resltant Force Magnitude
                        f_resultant_magnitude = np.linalg.norm(wrench_data_for_this_joint[:3])
                        max_f_key = ('F_resultant', j, 'max')
                        if f_resultant_magnitude > max_min_tracker[max_f_key][0]:
                            max_min_tracker[max_f_key] = [f_resultant_magnitude, wrench_data_for_this_joint.tolist(), q.tolist(), qd.tolist(), qdd.tolist()]

                        #Update Max Resltant Moment Magnitude
                        m_resultant_magnitude = np.linalg.norm(wrench_data_for_this_joint[3:])
                        max_m_key = ('M_resultant', j, 'max')
                        if m_resultant_magnitude > max_min_tracker[max_m_key][0]:
                            max_min_tracker[max_m_key] = [m_resultant_magnitude, wrench_data_for_this_joint.tolist(), q.tolist(), qd.tolist(), qdd.tolist()]

                    current_combination_count += 1

                except Exception as e:
                    pass



print(f"\nFinished processing {current_combination_count} total combinations.")

# construct Final Critical LoadCases Data
critical_load_cases_coupled = []

for j in range(robot.n):
    for comp_name in component_names_list:
        max_key = (comp_name, j, 'max')
        max_value, max_wrench, max_q, max_qd, max_qdd = max_min_tracker[max_key]
        if max_value != -float('inf'):
            row_data_max = {
                'ID': f'MAX_{comp_name}_Joint{j}',
                'Fx': max_wrench[0], 'Fy': max_wrench[1], 'Fz': max_wrench[2],
                'Mx': max_wrench[3], 'My': max_wrench[4], 'Mz': max_wrench[5],
                'q': max_q, 'qd': max_qd, 'qdd': max_qdd
            }
            critical_load_cases_coupled.append(row_data_max)

        min_key = (comp_name, j, 'min')
        min_value, min_wrench, min_q, min_qd, min_qdd = max_min_tracker[min_key]
        if min_value != float('inf'):
            row_data_min = {
                'ID': f'MIN_{comp_name}_Joint{j}',
                'Fx': min_wrench[0], 'Fy': min_wrench[1], 'Fz': min_wrench[2],
                'Mx': min_wrench[3], 'My': min_wrench[4], 'Mz': min_wrench[5],
                'q': min_q, 'qd': min_qd, 'qdd': min_qdd
            }
            critical_load_cases_coupled.append(row_data_min)

    # Max Resultant Force
    max_f_key = ('F_resultant', j, 'max')
    max_f_mag, max_f_wrench, max_f_q, max_f_qd, max_f_qdd = max_min_tracker[max_f_key]
    if max_f_mag != -float('inf'):
        row_data_max_f = {
            'ID': f'MAX_F_RESULTANT_Joint{j}', # e.g., MAX_F_RESULTANT_Joint0
            'Fx': max_f_wrench[0], 'Fy': max_f_wrench[1], 'Fz': max_f_wrench[2],
            'Mx': max_f_wrench[3], 'My': max_f_wrench[4], 'Mz': max_f_wrench[5],
            'q': max_f_q, 'qd': max_f_qd, 'qdd': max_f_qdd
        }
        critical_load_cases_coupled.append(row_data_max_f)

    # Max Resultant Moment
    max_m_key = ('M_resultant', j, 'max')
    max_m_mag, max_m_wrench, max_m_q, max_m_qd, max_m_qdd = max_min_tracker[max_m_key]
    if max_m_mag != -float('inf'):
        row_data_max_m = {
            'ID': f'MAX_M_RESULTANT_Joint{j}', # e.g., MAX_M_RESULTANT_Joint0
            'Fx': max_m_wrench[0], 'Fy': max_m_wrench[1], 'Fz': max_m_wrench[2],
            'Mx': max_m_wrench[3], 'My': max_m_wrench[4], 'Mz': max_m_wrench[5],
            'q': max_m_q, 'qd': max_m_qd, 'qdd': max_m_qdd
        }
        critical_load_cases_coupled.append(row_data_max_m)


df_critical_loads_coupled = pd.DataFrame(critical_load_cases_coupled)

# export
now = datetime.datetime.now()
time_str = now.strftime("%H%M%S")
output_filename = f"coupled_extreme_joint_wrenches_{now.date().isoformat()}_{time_str}.csv"

output_columns = ['ID', 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz', 'q', 'qd', 'qdd']

if not df_critical_loads_coupled.empty:
    df_critical_loads_coupled[output_columns].to_csv(output_filename, index=False)
    print(f"\nCoupled extreme joint wrenches for all DH frames exported to {output_filename}")

    #print("\n--- Important Notes for Generative Design (Coupled Extreme Load Cases) ---")
    #print("1. Each row in the CSV represents a fully coupled load case that produced an extreme value for one specific force or moment component (or resultant magnitude).")
    #print("2. The columns Fx, Fy, Fz, Mx, My, Mz provide the complete wrench (all 6 components) acting at the DH origin frame specified in the 'ID' column (e.g., 'Joint0').")
    #print("3. The wrench values are expressed in the LOCAL coordinate system of that specific DH origin frame.")
    #print("4. The (q, qd, qdd) values in that row represent the robot configuration that generated this specific coupled wrench.")
    #print("For each critical load row, apply ALL 6 force/moment components (Fx, Fy, Fz, Mx, My, Mz) at the origin of the corresponding DH frame (e.g., DH Frame 2 for 'Joint2').")

else:
    print("\nNo coupled extreme load cases to export because no valid combinations were processed.")