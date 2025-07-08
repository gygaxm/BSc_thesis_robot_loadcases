import roboticstoolbox as rtb
import numpy as np
import itertools
import pandas as pd
import datetime
import math

# --- Imports for the copied _modified_rne_python function ---
# These are internal imports of roboticstoolbox.
# They are required because _modified_rne_python uses them directly.
from spatialmath import SE3 # Used for base/tool transformations
from spatialmath.base.argcheck import getvector, isscalar, verifymatrix, getmatrix # Argument checking
from spatialmath.base import tr2jac, tr2eul, tr2rpy, t2r, trlog, rotvelxform # Transformations
import spatialmath.base.symbolic as sym # For symbolic models (used internally, keep for full compatibility)
# --- End of Imports for _modified_rne_python ---


# --- Helper function from original rne_python source ---
def _cross(a, b):
    # This is a helper function used internally by rne_python for cross product
    return np.r_[
        a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]
    ]

# --- MODIFIED rne_python function (copied from RoboticToolbox source) ---
# This function is now local to your script.
# Original source code for rne_python is here for reference:
# https://github.com/petercorke/roboticstoolbox/blob/master/roboticstoolbox/robot/DHRobot.py#L716
# Modifications: Added return values for all_f_wrenches and all_n_wrenches
# IMPORTANT: This function expects the robot object as its first argument (robot_obj)
# to access link properties and other robot attributes.
def _modified_rne_python(
    robot_obj, # Pass the robot object as the first argument (corresponds to 'self' in original)
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

    # Determine data type for arrays based on symbolic mode
    if self.symbolic:
        dtype = "O" # Object type for symbolic variables
    else:
        dtype = None # Default numerical type

    z0 = np.array([0, 0, 1], dtype=dtype) # Z-axis unit vector

    # Handle gravity input
    if gravity is None:
        gravity = self.gravity  # Use default gravity from the robot object
    else:
        gravity = getvector(gravity, 3) # Ensure gravity is a 3-element vector

    # Handle external force input (on end-effector)
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

    # Optional: Array to store base wrench (if base_wrench=True)
    if base_wrench:
        wbase = np.zeros((nk, n), dtype=dtype)

    # Lists to store internal forces (f) and moments (nn) for each link for each trajectory point.
    # all_f_wrenches[k][j] will store the force part of the wrench on link j at trajectory point k.
    # all_n_wrenches[k][j] will store the moment part of the wrench on link j at trajectory point k.
    all_f_wrenches = [[None for _ in range(n)] for _ in range(nk)]
    all_n_wrenches = [[None for _ in range(n)] for _ in range(nk)]

    # --- Main Loop for Trajectory Points (k) ---
    # This loop processes each set of (q, qd, qdd) from a trajectory.
    # For a single combination, nk will be 1, and this loop runs once.
    for k in range(nk):
        # Extract current joint states for this specific trajectory point
        q_k = q[k, :]
        qd_k = qd[k, :]
        qdd_k = qdd[k, :]

        if debug:
            print("q_k", q_k)
            print("qd_k", qd_k)
            print("qdd_k", qdd_k)
            print()

        # Initialize internal variables for forward recursion pass
        Fm = np.zeros((3, n), dtype=dtype) # Inertial force for each link's CoM
        Nm = np.zeros((3, n), dtype=dtype) # Inertial moment for each link's CoM
        pstarm = np.zeros((3, n), dtype=dtype) # Position vector from O_{j-1} to O_j in {j} frame
        Rm = [] # List of rotation matrices from link j-1 to link j

        # Initialize base angular and linear velocities/accelerations
        w = np.zeros((3,), dtype=dtype) # Angular velocity of link frame
        wd = np.zeros((3,), dtype=dtype) # Angular acceleration of link frame
        vd = -gravity  # Linear acceleration (initial due to gravity, transformed if base is set)

        # Apply base transformation to initial base velocities/accelerations
        if base is not None:
            Rb = t2r(base).T # Transpose of base rotation matrix (converts from world to base frame)
            w = Rb @ w
            wd = Rb @ wd
            vd = Rb @ gravity # Transform gravity vector into base frame if base transform exists

        # --- Initialize Link-specific variables for Forward Recursion ---
        for j in range(n):
            link = self.links[j]

            # Compute the link transformation matrix and determine 'd' value (prismatic vs revolute)
            if link.sigma == 0: # Revolute joint
                Tj = link.A(q_k[j]).A # A(q) computes the transformation for the link's joint angle
                d = link.d
            else: # Prismatic joint
                Tj = link.A(link.theta).A # Prismatic joint's rotation is fixed
                d = q_k[j] # For prismatic, 'd' is the variable joint coordinate

            # Compute pstar: vector from O_{j-1} to O_j in {j} frame
            alpha = link.alpha
            if self.mdh: # Modified DH convention
                pstar = np.r_[link.a, -d * sym.sin(alpha), d * sym.cos(alpha)]
                if j == 0 and base is not None: # Apply base transform for first link if base is present
                    Tj = base @ Tj
                    # pstar = base @ pstar # Original code has this, but pstar is a vector, not a pose
            else: # Standard DH convention
                pstar = np.r_[link.a, d * sym.sin(alpha), d * sym.cos(alpha)]

            # Store rotation matrix (from j-1 to j) and pstar vector for backward pass
            Rm.append(t2r(Tj)) # Rotation matrix from link j-1 to link j frame
            pstarm[:, j] = pstar # pstar vector for link j

        # -----------------  Forward Recursion (from base to end-effector) -------------------- #
        # Calculates angular/linear velocities/accelerations, and inertial forces/moments for each link
        for j, link in enumerate(self.links):
            Rt = Rm[j].T # Rotation matrix from link j to link j-1 (transpose of Rm[j])
            pstar = pstarm[:, j] # Position vector from O_{j-1} to O_j in {j} frame
            r = link.r # Center of Mass vector in link j frame

            # Update angular (w, wd) and linear (vd) velocities/accelerations for link j
            if self.mdh: # Modified DH
                if link.isrevolute:
                    w_ = Rt @ w + z0 * qd_k[j] # Angular velocity
                    wd_ = Rt @ wd + z0 * qdd_k[j] + _cross(Rt @ w, z0 * qd_k[j]) # Angular acceleration
                    vd_ = Rt @ (_cross(wd, pstar) + _cross(w, _cross(w, pstar)) + vd) # Linear acceleration
                else: # Prismatic
                    w_ = Rt @ w
                    wd_ = Rt @ wd
                    vd_ = (
                        Rt @ (_cross(wd, pstar) + _cross(w, _cross(w, pstar)) + vd)
                        + 2 * _cross(Rt @ w, z0 * qd_k[j]) # Coriolis term for prismatic
                        + z0 * qdd_k[j] # Linear acceleration due to joint motion
                    )
                w = w_
                wd = wd_
                vd = vd_
            else: # Standard DH
                if link.isrevolute:
                    wd = Rt @ (wd + z0 * qdd_k[j] + _cross(w, z0 * qd_k[j])) # Angular acceleration
                    w = Rt @ (w + z0 * qd_k[j]) # Angular velocity
                    vd = _cross(wd, pstar) + _cross(w, _cross(w, pstar)) + Rt @ vd # Linear acceleration
                else: # Prismatic
                    w = Rt @ w
                    wd = Rt @ wd
                    vd = (
                        Rt @ (z0 * qdd_k[j] + vd) # Propagated linear acceleration + joint acceleration
                        + _cross(wd, pstar) # Tangential acceleration due to angular acceleration
                        + 2 * _cross(w, Rt @ z0 * qd_k[j]) # Coriolis term for prismatic
                        + _cross(w, _cross(w, pstar)) # Centripetal acceleration
                    )

            # Calculate inertial forces (Fm) and moments (Nm) for link j's CoM
            vhat = _cross(wd, r) + _cross(w, _cross(w, r)) + vd # Linear acceleration of link's CoM
            Fm[:, j] = link.m * vhat # Inertial force acting at CoM (F = ma)
            Nm[:, j] = link.I @ wd + _cross(w, link.I @ w) # Inertial moment about CoM (N = I*alpha + omega x I*omega)

            if debug:
                print("w:     ", removesmall(w))
                print("wd:    ", removesmall(wd))
                print("vd:    ", removesmall(vd))
                print("vdbar: ", removesmall(vhat))
                print()

        if debug:
            print("Fm\n", Fm)
            print("Nm\n", Nm)

        # -----------------  Backward Recursion (from end-effector to base) -------------------- #
        # Propagates forces and moments, calculates joint torques.
        # f and nn are the net forces and moments acting on the distal end of the current link.
        f = fext[:3]  # Initial force component of end-effector wrench (in EE frame)
        nn = fext[3:] # Initial moment component of end-effector wrench (in EE frame)

        for j in range(n - 1, -1, -1): # Iterate from last link (n-1) down to first link (0)
            link = self.links[j]
            r = link.r # CoM vector in link j frame

            # Propagate forces (f) and moments (nn) from link j+1 to link j
            if self.mdh: # Modified DH convention
                if j == (n - 1): # If processing the last link
                    R = np.eye(3, dtype=dtype) # Identity rotation for end-effector to link
                    pstar = np.zeros((3,), dtype=dtype) # Zero vector from O_n to EE if no tool
                else: # For other links (not the last)
                    R = Rm[j + 1] # Rotation from link j to link j+1
                    pstar = pstarm[:, j + 1] # Vector from O_j to O_{j+1} in {j+1} frame

                f_ = R @ f + Fm[:, j] # Force: propagated from j+1 + current link's inertial force
                nn_ = (
                    R @ nn # Propagated moment from j+1
                    + _cross(pstar, R @ f) # Moment due to propagated force about new origin (O_j)
                    + _cross(pstar, Fm[:, j]) # Moment due to current link's inertial force about O_j
                    + Nm[:, j] # Current link's inertial moment
                )
                f = f_
                nn = nn_

            else: # Standard DH convention
                pstar = pstarm[:, j] # Vector from O_{j-1} to O_j in {j} frame
                if j == (n - 1): # If processing the last link
                    R = np.eye(3, dtype=dtype)
                else: # For other links
                    R = Rm[j + 1] # Rotation from link j to link j+1

                nn = (
                    R @ (nn + _cross(R.T @ pstar, f)) # Propagated moment + moment due to propagated force about O_j
                    + _cross(pstar + r, Fm[:, j]) # Moment due to current link's inertial force about O_j
                    + Nm[:, j] # Current link's inertial moment
                )
                f = R @ f + Fm[:, j] # Propagated force + current link's inertial force

            if debug:
                print("f: ", removesmall(f))
                print("n: ", removesmall(nn))

            # --- Store the extracted wrenches for each link ---
            # These are the forces (f) and moments (nn) acting ON THE CURRENT LINK (link j)
            # at its connection point to the *previous* link (joint j).
            # This is precisely the wrench transmitted through joint j from all subsequent links and current link's inertia.
            all_f_wrenches[k][j] = f.flatten() # Store force vector
            all_n_wrenches[k][j] = nn.flatten() # Store moment vector

            # Calculate joint torque (t) from the net moment (nn) or force (f)
            R = Rm[j] # Rotation matrix from j-1 to j
            if self.mdh:
                if link.isrevolute:
                    t = nn @ z0 # Moment about Z-axis for revolute (in link j frame)
                else: # Prismatic
                    t = f @ z0 # Force along Z-axis for prismatic (in link j frame)
            else: # Standard DH
                if link.isrevolute:
                    t = nn @ (R.T @ z0) # Moment about joint axis (transformed to j-1 frame for dot product)
                else: # Prismatic
                    t = f @ (R.T @ z0) # Force along joint axis (transformed to j-1 frame for dot product)

            # Add joint inertia and friction contributions to torque
            tau[k, j] = (
                t
                + link.G**2 * link.Jm * qdd_k[j] # Armature inertia contribution
                - link.friction(qd_k[j], coulomb=not self.symbolic) # Joint friction
            )
            if debug:
                print(
                    f"j={j:}, G={link.G:}, Jm={link.Jm:},"
                    f" friction={link.friction(qd_k[j], coulomb=False):}"
                )
                print()

        # Compute the base wrench if requested (wrench on link 0 from the base/ground)
        if base_wrench:
            R = Rm[0] # Rotation from base to link 0
            nn = R @ nn
            f = R @ f
            wbase[k, :] = np.r_[f, nn]

    # --- Return values based on requested output ---
    # The [0] index is used if only a single trajectory point (nk=1) was processed
    # to return a 1D array instead of a 2D array.
    if base_wrench:
        if tau.shape[0] == 1:
            return tau.flatten(), wbase.flatten(), all_f_wrenches[0], all_n_wrenches[0]
        else:
            return tau, wbase, all_f_wrenches, all_n_wrenches
    else:
        if tau.shape[0] == 1:
            return tau.flatten(), all_f_wrenches[0], all_n_wrenches[0]
        else:
            return tau, all_f_wrenches, all_n_wrenches


# --- 1. Robot Kinematic Model & Parameter Definition ---
# Replace with actual DH parameters for Stäubli TX2-40 or your 6-DOF robot
# Format: RevoluteDH(d, a, alpha, qlim=[min_angle, max_angle])
# d: link offset (m)
# a: link length (m)
# alpha: link twist (rad)
# qlim: joint angle limits (rad)

# These are placeholders! You MUST find the correct DH parameters for your robot.
# Example DH values (simplified for a generic 6-DOF, NOT accurate for Stäubli TX2-40)
# Refer to Peter Corke's book or Stäubli documentation for actual values.
# Joint 1 (base rotation)
L1 = rtb.RevoluteDH(d=0.3, a=0, alpha=math.pi / 2, qlim=[-math.pi, math.pi])
# Joint 2 (shoulder lift)
L2 = rtb.RevoluteDH(d=0, a=0.5, alpha=0, qlim=[-math.pi / 2, math.pi / 2])
# Joint 3 (elbow bend)
L3 = rtb.RevoluteDH(d=0, a=0.4, alpha=0, qlim=[-math.pi * 3 / 4, math.pi * 3 / 4])
# Joint 4 (wrist rotation)
L4 = rtb.RevoluteDH(d=0.4, a=0, alpha=math.pi / 2, qlim=[-math.pi, math.pi])
# Joint 5 (wrist bend)
L5 = rtb.RevoluteDH(d=0, a=0, alpha=-math.pi / 2, qlim=[-math.pi / 2, math.pi / 2])
# Joint 6 (tool flange rotation)
L6 = rtb.RevoluteDH(d=0.1, a=0, alpha=0, qlim=[-math.pi, math.pi])

robot = rtb.DHRobot([L1, L2, L3, L4, L5, L6], name="Generic6DOF")

# --- 2. Robot Inertial Properties (Simplified Placeholders) ---
# You MUST replace these with realistic estimates for your robot links.
# Assign mass (m), center of mass (r, a 3-element vector), and inertia tensor (I, a 3x3 matrix)
# Inertia tensor is typically about the link's center of mass, aligned with its local frame.

# Example properties for links (mass in kg, r in m, I in kg*m^2)
# For Link 1
robot.links[0].m = 5.0  # kg
robot.links[0].r = np.array([0, 0, L1.d / 2])  # CoM for a base link
robot.links[0].I = np.diag([0.1, 0.1, 0.05])  # Inertia tensor (simplified diagonal)

# For Link 2
robot.links[1].m = 10.0
robot.links[1].r = np.array([L2.a / 2, 0, 0])  # CoM in link frame
robot.links[1].I = np.diag([0.2, 0.8, 0.8])

# For Link 3 (the one you want to redesign - between J3 and J4)
# Make these values as accurate as possible for the target link
robot.links[2].m = 8.0  # Example mass for Link 3
robot.links[2].r = np.array([L3.a / 2, 0, 0])  # CoM for Link 3
robot.links[2].I = np.diag([0.1, 0.5, 0.5])  # Example inertia for Link 3

# For Link 4
robot.links[3].m = 3.0
robot.links[3].r = np.array([0, 0, L4.d / 2])
robot.links[3].I = np.diag([0.05, 0.05, 0.02])

# For Link 5
robot.links[4].m = 2.0
robot.links[4].r = np.array([0, 0, 0])  # Assuming CoM at joint for small wrist link
robot.links[4].I = np.diag([0.02, 0.02, 0.01])

# For Link 6 (with end-effector/payload)
# Add end-effector and payload mass/inertia to the last link
payload_mass = 20.0  # kg
end_effector_mass = 1.0  # kg
end_effector_length = 0.1  # m (from last joint to payload CoM)
robot.links[5].m = 1.0 + payload_mass + end_effector_mass  # Combined mass
robot.links[5].r = np.array([end_effector_length / 2, 0, 0])  # CoM adjusted for end-effector
# Simplified inertia for Link 6 + end-effector + payload
robot.links[5].I = np.diag([0.01, 0.01, 0.005]) + np.diag(
    [0, payload_mass * (end_effector_length / 2) ** 2, payload_mass * (end_effector_length / 2) ** 2])

# --- Set Microgravity Environment ---
robot.gravity = np.array([0, 0, 0])  # No gravity vector


# --- 3. Define Sweep Parameters ---
# Number of steps for each joint's q, qd, qdd
# WARNING: The total number of combinations is N_q^6 * N_qd^6 * N_qdd^6
# Start with small numbers for testing, then increase strategically.
num_q_steps = 10
num_qd_steps = 2  # e.g., max_neg_vel, 0, max_pos_vel
num_qdd_steps = 2  # e.g., max_neg_accel, max_pos_accel

q_ranges = []
qd_ranges = []
qdd_ranges = []

for i, link in enumerate(robot.links):
    # Joint angle range (q)
    q_ranges.append(np.linspace(link.qlim[0], link.qlim[1], num_q_steps))

    # Joint velocity range (qd) - estimate reasonable max values if not in DH
    # Adjust these max velocity/acceleration values to your robot's actual limits
    qd_max = math.pi*2  # rad/s, placeholder for max angular velocity
    qdd_max = math.pi * 4  # rad/s^2, placeholder for max angular acceleration

    qd_ranges.append(np.linspace(-qd_max, qd_max, num_qd_steps))
    qdd_ranges.append(np.linspace(-qdd_max, qdd_max, num_qdd_steps))

# Combine all ranges for iteration
all_q_combinations = itertools.product(*q_ranges)
all_qd_combinations = itertools.product(*qd_ranges)
all_qdd_combinations = itertools.product(*qdd_ranges)

# --- 4. Dynamic Calculation and Load Case Extraction ---
all_extracted_wrenches = []
link_to_redesign_index = 2  # Link 3 (robot.links[2]) is between Joint 3 and Joint 4

total_combinations = num_q_steps ** 6 * num_qd_steps ** 6 * num_qdd_steps ** 6
print(f"Total combinations to simulate: {total_combinations}")
if total_combinations > 1e9:
    print("WARNING: This many combinations will take a very long time. Consider reducing num_steps or sampling.")

current_combination_count = 0
# Use a more efficient way to combine the iterators
for q_tuple in all_q_combinations:
    for qd_tuple in all_qd_combinations:
        for qdd_tuple in all_qdd_combinations:
            q = np.array(q_tuple)
            qd = np.array(qd_tuple)
            qdd = np.array(qdd_tuple)

            try:
                # Calculate wrenches (forces and moments) for all links using our modified function
                # The returned all_f_wrenches_per_k and all_n_wrenches_per_k are lists of arrays
                # where each element corresponds to a link's wrench in its own frame.
                # Since nk is 1 in this loop (single q, qd, qdd combo), we get the first element [0]
                _, all_f_wrenches_per_k, all_n_wrenches_per_k = _modified_rne_python(
                    robot, q, qd, qdd, gravity=robot.gravity
                )

                # --- 1. Wrench on Joint 3 side of Link 3 (from Joint 3) ---
                # This is the wrench applied *to Link 3* by Joint 3.
                # It's already expressed in Link 3's local coordinate frame.
                # These are the direct loads Link 3 experiences from the previous part of the arm.
                wrench_J3_on_L3_force = all_f_wrenches_per_k[link_to_redesign_index]
                wrench_J3_on_L3_moment = all_n_wrenches_per_k[link_to_redesign_index]

                # --- 2. Wrench on Joint 4 side of Link 3 (from Joint 4) ---
                # This is the wrench applied *to Link 3* by Joint 4.
                # It's the reaction to the wrench applied *to Link 4 by Joint 4*.
                # All values will be in Link 3's local coordinate frame.

                # Get the wrench applied *to Link 4 by Joint 4* (expressed in Link 4's frame)
                # This is all_f_wrenches_per_k[link_to_redesign_index + 1] and all_n_wrenches_per_k[link_to_redesign_index + 1]
                wrench_J4_on_L4_force_L4_frame = all_f_wrenches_per_k[link_to_redesign_index + 1]
                wrench_J4_on_L4_M_L4_frame = all_n_wrenches_per_k[link_to_redesign_index + 1]

                # Calculate the reaction wrench (negative of the load Link 4 experiences from J4)
                # This reaction is still in Link 4's coordinate frame for now.
                reaction_J4_on_L3_force_L4_frame = -wrench_J4_on_L4_force_L4_frame
                reaction_J4_on_L3_M_L4_frame = -wrench_J4_on_L4_M_L4_frame # Corrected to use M_L4_frame

                # --- Transform this reaction wrench from Link 4's frame to Link 3's frame ---
                # Get transformation from base to Link 3, and base to Link 4
                T_base_to_L3 = robot.fkine_all(q)[link_to_redesign_index] # Pose of L3's frame relative to base
                T_base_to_L4 = robot.fkine_all(q)[link_to_redesign_index + 1] # Pose of L4's frame relative to base

                # Get transformation from Link 3 to Link 4 (T_L3_to_L4 = T_base_to_L3_inv * T_base_to_L4)
                T_L3_to_L4 = T_base_to_L3.inv() * T_base_to_L4

                # Get rotation matrix from Link 4's frame to Link 3's frame (R_L4_to_L3 = inverse of R_L3_to_L4)
                # The inverse of a rotation matrix is its transpose.
                R_L4_to_L3 = T_L3_to_L4.R.T

                # Rotate the force and moment vectors into Link 3's frame
                wrench_J4_on_L3_force = R_L4_to_L3 @ reaction_J4_on_L3_force_L4_frame
                wrench_J4_on_L3_moment = R_L4_to_L3 @ reaction_J4_on_L3_M_L4_frame # Corrected to use M_L4_frame

                # --- Store both wrenches for this combination ---
                # These are the actual loads you would apply in generative design.
                # They are all expressed in Link 3's local coordinate frame.
                all_extracted_wrenches.append({
                    'q': q.tolist(),
                    'qd': qd.tolist(),
                    'qdd': qdd.tolist(),
                    # Loads at Joint 3 side of Link 3 (proximal end)
                    'Fx_J3_on_L3': wrench_J3_on_L3_force[0],
                    'Fy_J3_on_L3': wrench_J3_on_L3_force[1],
                    'Fz_J3_on_L3': wrench_J3_on_L3_force[2],
                    'Mx_J3_on_L3': wrench_J3_on_L3_moment[0],
                    'My_J3_on_L3': wrench_J3_on_L3_moment[1],
                    'Mz_J3_on_L3': wrench_J3_on_L3_moment[2],
                    # Loads at Joint 4 side of Link 3 (distal end)
                    'Fx_J4_on_L3': wrench_J4_on_L3_force[0],
                    'Fy_J4_on_L3': wrench_J4_on_L3_force[1],
                    'Fz_J4_on_L3': wrench_J4_on_L3_force[2],
                    'Mx_J4_on_L3': wrench_J4_on_L3_moment[0],
                    'My_J4_on_L3': wrench_J4_on_L3_moment[1],
                    'Mz_J4_on_L3': wrench_J4_on_L3_moment[2]
                })
            except Exception as e:
                print(f"Error calculating wrench for state {q}, {qd}, {qdd}: {e}")

            current_combination_count += 1
            if current_combination_count % 1000 == 0:
                print(f"Processed {current_combination_count} combinations...", end='\r')

print(f"\nFinished processing {len(all_extracted_wrenches)} valid combinations.")

# Convert to DataFrame for easier analysis
df_wrenches = pd.DataFrame(all_extracted_wrenches)

# --- 5. Load Case Filtering ---
# Find extreme load cases for generative design.
# We'll find max/min for individual components on BOTH sides, and max resultants.

critical_load_cases = []

# List of columns for Joint 3 side and Joint 4 side loads
cols_J3 = ['Fx_J3_on_L3', 'Fy_J3_on_L3', 'Fz_J3_on_L3', 'Mx_J3_on_L3', 'My_J3_on_L3', 'Mz_J3_on_L3']
cols_J4 = ['Fx_J4_on_L3', 'Fy_J4_on_L3', 'Fz_J4_on_L3', 'Mx_J4_on_L3', 'My_J4_on_L3', 'Mz_J4_on_L3']

# Check if DataFrame is empty before proceeding
if df_wrenches.empty:
    print("\nNo critical load cases to extract because no valid combinations were processed.")
else:
    # Max/Min for each component on Joint 3 side
    for col in cols_J3:
        max_row = df_wrenches.loc[df_wrenches[col].idxmax()].to_dict()
        min_row = df_wrenches.loc[df_wrenches[col].idxmin()].to_dict()

        critical_load_cases.append({**{'ID': f'MAX_{col}'}, **{c: max_row[c] for c in cols_J3 + cols_J4},
                                    'q': max_row['q'], 'qd': max_row['qd'], 'qdd': max_row['qdd']})
        critical_load_cases.append({**{'ID': f'MIN_{col}'}, **{c: min_row[c] for c in cols_J3 + cols_J4},
                                    'q': min_row['q'], 'qd': min_row['qd'], 'qdd': min_row['qdd']})

    # Max/Min for each component on Joint 4 side
    for col in cols_J4:
        max_row = df_wrenches.loc[df_wrenches[col].idxmax()].to_dict()
        min_row = df_wrenches.loc[df_wrenches[col].idxmin()].to_dict()

        critical_load_cases.append({**{'ID': f'MAX_{col}'}, **{c: max_row[c] for c in cols_J3 + cols_J4},
                                    'q': max_row['q'], 'qd': max_row['qd'], 'qdd': max_row['qdd']})
        critical_load_cases.append({**{'ID': f'MIN_{col}'}, **{c: min_row[c] for c in cols_J3 + cols_J4},
                                    'q': min_row['q'], 'qd': min_row['qd'], 'qdd': min_row['qdd']})

    # Max Resultant Force Magnitude (considering sum of forces from both ends for overall stress)
    # Note: For generative design, you typically apply loads at specific interfaces.
    # The 'resultant' here means the highest magnitude across all the Fx, Fy, Fz from both J3 and J4 sides.
    # It doesn't mean the vector sum of F_J3_on_L3 and F_J4_on_L3, as those are applied at different points.
    df_wrenches['F_resultant_J3'] = np.linalg.norm(df_wrenches[['Fx_J3_on_L3', 'Fy_J3_on_L3', 'Fz_J3_on_L3']].values, axis=1)
    df_wrenches['F_resultant_J4'] = np.linalg.norm(df_wrenches[['Fx_J4_on_L3', 'Fy_J4_on_L3', 'Fz_J4_on_L3']].values, axis=1)
    df_wrenches['Max_F_Overall'] = df_wrenches[['F_resultant_J3', 'F_resultant_J4']].max(axis=1)
    max_F_overall_row = df_wrenches.loc[df_wrenches['Max_F_Overall'].idxmax()].to_dict()
    critical_load_cases.append({**{'ID': 'MAX_F_OVERALL_MAGNITUDE'}, **{c: max_F_overall_row[c] for c in cols_J3 + cols_J4},
                                'q': max_F_overall_row['q'], 'qd': max_F_overall_row['qd'], 'qdd': max_F_overall_row['qdd']})

    # Max Resultant Moment Magnitude (overall)
    df_wrenches['M_resultant_J3'] = np.linalg.norm(df_wrenches[['Mx_J3_on_L3', 'My_J3_on_L3', 'Mz_J3_on_L3']].values, axis=1)
    df_wrenches['M_resultant_J4'] = np.linalg.norm(df_wrenches[['Mx_J4_on_L3', 'My_J4_on_L3', 'Mz_J4_on_L3']].values, axis=1)
    df_wrenches['Max_M_Overall'] = df_wrenches[['M_resultant_J3', 'M_resultant_J4']].max(axis=1)
    max_M_overall_row = df_wrenches.loc[df_wrenches['Max_M_Overall'].idxmax()].to_dict()
    critical_load_cases.append({**{'ID': 'MAX_M_OVERALL_MAGNITUDE'}, **{c: max_M_overall_row[c] for c in cols_J3 + cols_J4},
                                'q': max_M_overall_row['q'], 'qd': max_M_overall_row['qd'], 'qdd': max_M_overall_row['qdd']})


# Convert critical load cases to DataFrame for export
df_critical_loads = pd.DataFrame(critical_load_cases)

# --- 6. Export for nTop ---
output_filename = f"load_cases_link{link_to_redesign_index + 1}_{datetime.date.today().isoformat()}.csv"

# Define the full set of output columns explicitly for clarity
output_columns = [
    'ID',
    'Fx_J3_on_L3', 'Fy_J3_on_L3', 'Fz_J3_on_L3',
    'Mx_J3_on_L3', 'My_J3_on_L3', 'Mz_J3_on_L3',
    'Point_of_Application_J3_X', 'Point_of_Application_J3_Y', 'Point_of_Application_J3_Z', # New point for J3
    'Fx_J4_on_L3', 'Fy_J4_on_L3', 'Fz_J4_on_L3',
    'Mx_J4_on_L3', 'My_J4_on_L3', 'Mz_J4_on_L3',
    'Point_of_Application_J4_X', 'Point_of_Application_J4_Y', 'Point_of_Application_J4_Z', # New point for J4
    'q', 'qd', 'qdd'
]

# Add placeholder for Point_of_Application for nTop for BOTH ends
# Point of Application J3 side: origin of Link 3's coordinate frame
df_critical_loads['Point_of_Application_J3_X'] = 0.0
df_critical_loads['Point_of_Application_J3_Y'] = 0.0
df_critical_loads['Point_of_Application_J3_Z'] = 0.0

# Point of Application J4 side: distal end of Link 3 (where Joint 4 connects)
# This would be at (L3.a, 0, 0) if Link 3's frame is aligned with its length.
# Adjust this based on your actual Link 3 geometry relative to its frame.
df_critical_loads['Point_of_Application_J4_X'] = robot.links[link_to_redesign_index].a
df_critical_loads['Point_of_Application_J4_Y'] = 0.0
df_critical_loads['Point_of_Application_J4_Z'] = 0.0


if not df_critical_loads.empty:
    df_critical_loads[output_columns].to_csv(output_filename, index=False)
    print(f"\nCritical load cases for Link {link_to_redesign_index + 1} exported to {output_filename}")
    print("\nSample of Critical Load Cases:")
    print(df_critical_loads[output_columns].head())
    print("\n--- Important Notes for Generative Design ---")
    print("1. All Force (F) and Moment (M) values are in Link 3's LOCAL coordinate frame.")
    print("2. For each load case (each row in the CSV):")
    print("   - Apply 'Fx_J3_on_L3' to 'Mz_J3_on_L3' at Point_of_Application_J3 (origin of Link 3's frame).")
    print("   - Apply 'Fx_J4_on_L3' to 'Mz_J4_on_L3' at Point_of_Application_J4 (distal end of Link 3).")
    print("3. You MUST apply boundary conditions (constraints) at ONE of these ends in your generative design software to prevent rigid body motion. The other end will receive loads.")
    print("   - For example: Apply a 'Fixed Support' or 'Hinge Constraint' at the Joint 3 interface (Point_of_Application_J3), and then apply the loads (Fx_J4_on_L3 to Mz_J4_on_L3) at the Joint 4 interface (Point_of_Application_J4).")
    print("   - Alternatively, fix the Joint 4 interface and apply loads at the Joint 3 interface. The choice depends on your specific analysis setup.")
    print("4. These loads already include the inertial effects of Link 3 itself and all downstream components.")

else:
    print("\nNo critical load cases to export because no valid combinations were processed.")

    # --- Understanding Link 3's Local Coordinate Frame (DH Convention) ---
    # The forces and moments (wrenches) calculated for Link 3 are expressed in its own local coordinate frame.
    # This frame's orientation follows the Denavit-Hartenberg (DH) convention:
    #
    # 1. Z-axis:
    #    [cite_start]- Aligned with the axis of rotation of Joint 3 (the joint preceding Link 3)[cite: 1, 2].
    #    [cite_start]- For typical industrial robot elbow joints, this axis is generally horizontal[cite: 2].
    #
    # 2. X-axis:
    #    [cite_start]- Points along the length of Link 3, extending from Joint 3 towards Joint 4[cite: 1, 2].
    #    [cite_start]- The 'a' (link length) DH parameter for Link 3 is defined along this X-axis[cite: 1].
    #
    # 3. Y-axis:
    #    [cite_start]- Completes a right-handed coordinate system with the X and Z axes[cite: 1, 2].
    #    - Its direction is determined by the cross-product of Z and X (Y = Z x X).
    #
    # Note on "Left/Right":
    # Whether the Z-axis (or Y-axis) points "left" or "right" in a global sense depends on the specific
    # DH parameters chosen for your robot (which define the positive direction of each joint's axis)
    # [cite_start]and the current configuration of the preceding joints (Joint 1 and Joint 2)[cite: 2].
    #