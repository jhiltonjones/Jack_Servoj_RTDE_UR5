import math 
import numpy as np
from constants import *
from magnetic import *
import scipy.integrate as spi

def compute_bending_moment_l(l_c, theta_F, theta_l, F_m, T_m):
    return (l_c * np.linalg.norm(F_m) * np.sin(theta_F - theta_l) + 2 * np.linalg.norm(T_m)) / 2

def integrand(theta):
    term = (M_l**2 / (2 * EI)) + Fy * np.sin(theta) + Fx * np.cos(theta)
    return term**(-0.5)

def compute_deflection(theta_max):
    integral_result, _ = spi.quad(integrand, 0, theta_max)
    x = np.sqrt(EI / 2) * integral_result * np.cos(theta_max)
    y = np.sqrt(EI / 2) * integral_result * np.sin(theta_max)
    return x, y
def compute_curvature(theta_l, length_c):
    """Computes curvature kappa based on bending angle and catheter length."""
    return np.abs(theta_l) / length_c  # Ensure positive curvature


# Volume_catheter = volume_calculator_cyclinder(s_d, length_c)
# M_catheter, M_catheter_scaler = magneitc_moment(mu_0, B_ra, Volume_catheter, m_c_hat)
# Volume_magnet = volume_calculator_cyclinder(s_a, h_a)
# M_mag, M_mag_scaler = magneitc_moment(mu_0, B_ra, Volume_magnet, m_a_hat)
# B = magnetic_field_external_magnet(mu_0, p_hat, M_mag_scaler)
# B_norm = np.linalg.norm(B)*1e3
# T_m = magnetic_torque(M_catheter, B)
# F_m = magnetic_force(mu_0, p, M_mag, M_catheter)
# M_l = compute_bending_moment_l(length_c, theta_f, theta_l, F_m, T_m)
# Fx = np.linalg.norm(F_m) * np.sin(theta_l)
# Fy = np.linalg.norm(F_m) * np.cos(theta_l)
# C = C = (M_l**2) / (2 * EI) + Fy * np.sin(theta_l) + Fx * np.cos(theta_l)
# print("Volume of Catheter: ", Volume_catheter)
# print("Magnetic moment of Catheter:", M_catheter_scaler)
# print("Magnetic moment of Magnet:", M_mag_scaler)
# print("Magnetic field: ", B)
# print("Magnetic field norm (mT): ", B_norm)
# print("Torque: ", T_m)
# print("Force: ", F_m)


