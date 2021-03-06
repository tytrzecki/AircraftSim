function [T_phi_delta_a,T_chi_phi,T_theta_delta_e,T_h_theta,T_h_Va,T_Va_delta_t,T_Va_theta,T_v_delta_r]...
    = compute_tf_model(x_trim,u_trim,P)
% x_trim is the trimmed state,
% u_trim is the trimmed input

% add stuff here
Gamma = P.Jx*P.Jz-P.Jxz^2;
Gamma3 = P.Jz/Gamma;
Gamma4 = P.Jxz/Gamma;

C_P_p = Gamma3*P.C_ell_p + Gamma4*P.C_n_p;
C_P_delta_a = Gamma3*P.C_ell_delta_a + Gamma4*P.C_n_delta_a;
a_phi1 = -.5*P.rho*P.Va^2*P.S_wing*P.b*C_P_p*P.b/2/P.Va;
a_phi2 = .5*P.rho*P.Va^2*P.S_wing*P.b*C_P_delta_a;

Va_trim = sqrt(P.u0^2 + P.v0^2 + P.w0^2);
a_theta1 = -P.rho*P.Va^2*P.c*P.S_wing/2/P.Jy*P.C_M_q*P.c/2/P.Va;
a_theta2 = -P.rho*P.Va^2*P.c*P.S_wing/2/P.Jy*P.C_M_alpha;
a_theta3 = P.rho*P.Va^2*P.c*P.S_wing/2/P.Jy*P.C_M_delta_e;

alpha_trim = atan(P.w0/P.u0);
theta_trim = P.theta0;
C_D_alpha = (P.C_L_0 + P.C_L_alpha*alpha_trim)^2/pi/P.e/(P.b^2/P.S_wing);
a_V1 = P.rho*Va_trim*P.S_wing/P.mass*(P.C_D_0 + C_D_alpha*alpha_trim + P.C_D_delta_e*u_trim(1)) ...
    + P.rho*P.S_prop/P.mass*P.C_prop*Va_trim;
a_V2 = P.rho*P.S_prop/P.mass*P.C_prop*P.k_motor^2*u_trim(4);
a_V3 = P.gravity*cos(theta_trim - alpha_trim); 

a_beta1 = -P.rho*P.Va*P.S_wing/2/P.mass*P.C_Y_beta;
a_beta2 = P.rho*P.Va*P.S_wing/2/P.mass*P.C_Y_delta_r;
    
% define transfer functions
T_phi_delta_a   = tf([a_phi2],[1,a_phi1,0]);
T_chi_phi       = tf([P.gravity/Va_trim],[1,0]);
T_theta_delta_e = tf(a_theta3,[1,a_theta1,a_theta2]);
T_h_theta       = tf([Va_trim],[1,0]);
T_h_Va          = tf([theta_trim],[1,0]);
T_Va_delta_t    = tf([a_V2],[1,a_V1]);
T_Va_theta      = tf([-a_V3],[1,a_V1]);
T_v_delta_r     = tf([a_beta2],[1,-a_beta1]);
