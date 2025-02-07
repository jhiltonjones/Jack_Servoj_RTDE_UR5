import numpy as np
from scipy.optimize import minimize
from ik_modelling import compute_unit_position_vector, components, moment_cath1, compute_center_of_catheter, compute_T_m
from constants import *
from kinematic import compute_curvature
def objective_function(alpha_o, EI, theta_l, length_c, x_p, y_p, x_distance, y_distance, p_norm):

    moment_cath_unit, moment_mag_unit = moment_cath1(x_distance, y_distance, alpha_o)
    n0 = (mu_0 * moment_mag_unit * moment_cath_unit) / (4 * np.pi * p_norm**3)
    print("N0: ", n0)
    n1, n2, n3, n4 = components(x_p, y_p)

    T_m= compute_T_m(n0, theta_l, alpha_o, n1, n2, n3, n4)
    print("Torque: ", T_m)
    # T_m = T_m_vector[2]
    T_e = (EI * theta_l) / length_c  
    print("Bending Moment: ", T_e)

    residual = T_e - T_m  
    alpha_o_scalar = alpha_o.item()
    print(f"Evaluating α_o = {np.degrees(alpha_o_scalar):.3f}° → Residual Norm = {np.linalg.norm(residual):.6f}")




    return np.linalg.norm(residual)

def find_optimal_alpha_o(EI, theta_l, length_c, x_p, y_p, x_distance, y_distance, p_norm, initial_guess=np.radians(5)):

    result = minimize(objective_function, initial_guess, args=(EI, theta_l, length_c, x_p, y_p, x_distance, y_distance, p_norm), method='BFGS')
    
    return result.x[0]

  
# if __name__ == "__main__":
kappa = compute_curvature(theta_l, length_c)
print(f"Curvature kappa: {kappa}")
x_c, y_c = compute_center_of_catheter(length_c, kappa, theta_l)
print("Center of catheter: ", x_c, y_c)
p_hat = compute_unit_position_vector(x_c, y_c, d, h)
x_p, y_p, _ = p_hat  
print("Positions of catheter: ", x_p, y_p)
p_norm = np.linalg.norm(p_hat)
alpha_o_star = find_optimal_alpha_o(EI, theta_l, length_c, x_p, y_p, x_distance, y_distance, p_norm)
print(f"Optimal α_o* (degrees): {np.degrees(alpha_o_star):.3f}°")
alpha_o_degrees = np.degrees(alpha_o_star)
rotation_angle_radians = alpha_o_star