% forces_moments.m
%   Computes the foreces and moments acting on the airframe. 
%
%   Output is
%       F     - forces
%       M     - moments
%       Va    - airspeed
%       alpha - angle of attack
%       beta  - sideslip angle
%       wind  - wind vector in the inertial frame
%
%  Revised:
%   2/2/2010 - RB 
%   5/14/2010 - RB

function out = forces_moments(x, delta, wind, P)

    % relabel the inputs
    pn      = x(1);
    pe      = x(2);
    pd      = x(3);
    u       = x(4);
    v       = x(5);
    w       = x(6);
    phi     = x(7);
    theta   = x(8);
    psi     = x(9);
    p       = x(10);
    q       = x(11);
    r       = x(12);
    delta_e = delta(1);
    delta_a = delta(2);
    delta_r = delta(3);
    delta_t = delta(4);
    w_ns    = wind(1); % steady wind - North
    w_es    = wind(2); % steady wind - East
    w_ds    = wind(3); % steady wind - Down
    u_wg    = wind(4); % gust along body x-axis
    v_wg    = wind(5); % gust along body y-axis    
    w_wg    = wind(6); % gust along body z-axis
    
    %Rotation matracies
    Rb_v2 = [...
          1, 0, 0;...
          0, cos(phi), sin(phi);...
          0, -sin(phi), cos(phi)];
    Rv2_v1 = [...
          cos(theta), 0, -sin(theta);...
          0, 1, 0;...
          sin(theta), 0, cos(theta)];
    Rv1_v = [...
          cos(psi), sin(psi), 0;...
          -sin(psi), cos(psi), 0;...
          0, 0, 1];
    Rb_v = Rb_v2*Rv2_v1*Rv1_v;
    
    % compute wind vector in the inertial frame
    Wg_ned = Rb_v'*[u_wg; v_wg; w_wg];
    w_n = w_ns + Wg_ned(1);
    w_e = w_es + Wg_ned(2);
    w_d = w_ds + Wg_ned(3);
    
    Vb_a = [u; v; w] - Rb_v*[w_ns; w_es; w_ds] + [u_wg; v_wg; w_wg];
    
    % compute airspeed Va, angle-of-attack alpha, side-slip beta
    Va    = sqrt(Vb_a(1)^2 + Vb_a(2)^2 + Vb_a(3)^2);
    alpha = atan(Vb_a(3)/Vb_a(1));
    beta  = asin(Vb_a(2)/Va);

    %Definitions
    sigma = (1+exp(-P.M*(alpha-P.alpha0))+exp(P.M*(alpha-P.alpha0)))/((1+exp(-P.M*(alpha-P.alpha0)))*(1+exp(P.M*(alpha+P.alpha0))));
    C_L = (1-sigma)*(P.C_L_0+P.C_L_alpha*alpha) + sigma*(2*sign(alpha)*(sin(alpha)^2)*cos(alpha));
    C_D = P.C_D_0 + (P.C_L_0 + P.C_L_alpha*alpha)^2/(pi*P.e*(P.b^2)/P.S_wing);
    C_X = -C_D*cos(alpha) + C_L*sin(alpha);
    C_X_q = -P.C_D_q*cos(alpha) + P.C_L_q*sin(alpha);
    C_X_delta_e = -P.C_D_delta_e*cos(alpha) + P.C_L_delta_e*sin(alpha);
    C_Z = -C_D*sin(alpha) - C_L*cos(alpha);
    C_Z_q = -P.C_D_q*sin(alpha) - P.C_L_q*cos(alpha);
    C_Z_delta_e = -P.C_D_delta_e*sin(alpha) - P.C_L_delta_e*cos(alpha);
    
    mg = P.mass*P.gravity;
    Force = [-mg*sin(theta);...
             mg*cos(theta)*sin(phi);...
             mg*cos(theta)*cos(phi)] + ...
             .5*P.rho*Va^2*P.S_wing*[C_X + C_X_q*P.c*q*.5/Va + C_X_delta_e*delta_e;...
                 P.C_Y_0+P.C_Y_beta*beta+P.C_Y_p*P.b*p*.5/Va+P.C_Y_r*P.b*r*.5/Va+P.C_Y_delta_a*delta_a+P.C_Y_delta_r*delta_r;...
                                     C_Z + C_Z_q*P.c*q*.5/Va + C_Z_delta_e*delta_e] + ...
             .5*P.rho*P.S_prop*P.C_prop*[(P.k_motor*delta_t)^2 - Va^2; 0; 0];
    
    Torque = .5*P.rho*Va^2*P.S_wing*...
             [P.b*(P.C_ell_0+P.C_ell_beta*beta+P.C_ell_p*P.b*p*.5/Va+P.C_ell_r*P.b*r*.5/Va+P.C_ell_delta_a*delta_a+P.C_ell_delta_r*delta_r);...
              P.c*(P.C_M_0+P.C_M_alpha*alpha+P.C_M_q*P.c*q*.5/Va+P.C_M_delta_e*delta_e);...
              P.b*(P.C_n_0+P.C_n_beta*beta+P.C_n_p*P.b*p*.5/Va+P.C_n_r*P.b*r*.5/Va+P.C_n_delta_a*delta_a+P.C_n_delta_r*delta_r)] +...
              [-P.k_T_p*(P.k_omega*delta_t)^2; 0; 0];

    out = [Force; Torque; Va; alpha; beta; w_n; w_e; w_d];

end



