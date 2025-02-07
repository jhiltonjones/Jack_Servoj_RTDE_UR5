import math
import numpy as np
mu_0 = 4* math.pi * 1e-7
length = 0.0192 #m permeability of free space
length_m = 0.007 #m length of marker 2
length_c = 10 #m legnth of magnetic catheter
s_d = 0.0025 #m diamter of marker 1
s_c = 0.001 #m diameter of catheter
EI  = 2.14e-6 #N*m**2 youngs modulus
B_ra = 1.35 #T Residual flux density of external magnet
B_rc = 1.44 #T Residual flux density fo magnetic catheter
s_a = 0.025
h_a = 0.025 #m diameter and height of of the driving magent
d_off = 0.099 #m offset distance between {E} and the flange
v_r = 0.006 #m/s insertion speed
epsilon = 0.0012 #m detection threshold
I = np.eye(3)
theta_l = np.pi / 4
theta_f = np.pi / 6
# p = [0,0.05,0.05]
# p_norm = np.linalg.norm(p)
m_c_hat = np.array([1,0,0])
m_a_hat = np.array([0,0,1])

d, h = 16, -75
y_distance = 100
x_distance = 100