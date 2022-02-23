% estimate_states
%   - estimate the MAV states using gyros, accels, pressure sensors, and
%   GPS.
%
% Outputs are:
%   pnhat    - estimated North position, 
%   pehat    - estimated East position, 
%   hhat     - estimated altitude, 
%   Vahat    - estimated airspeed, 
%   alphahat - estimated angle of attack
%   betahat  - estimated sideslip angle
%   phihat   - estimated roll angle, 
%   thetahat - estimated pitch angel, 
%   chihat   - estimated course, 
%   phat     - estimated roll rate, 
%   qhat     - estimated pitch rate, 
%   rhat     - estimated yaw rate,
%   Vghat    - estimated ground speed, 
%   wnhat    - estimate of North wind, 
%   wehat    - estimate of East wind
%   psihat   - estimate of heading angle
% 
% 
% Modified:  3/15/2010 - RB
%            5/18/2010 - RB
%

function xhat = estimate_states(uu, P)

   % rename inputs
   y_gyro_x      = uu(1);
   y_gyro_y      = uu(2);
   y_gyro_z      = uu(3);
   y_accel_x     = uu(4);
   y_accel_y     = uu(5);
   y_accel_z     = uu(6);
   y_static_pres = uu(7);
   y_diff_pres   = uu(8);
   y_gps_n       = uu(9);
   y_gps_e       = uu(10);
   y_gps_h       = uu(11);
   y_gps_Vg      = uu(12);
   y_gps_course  = uu(13);
   t             = uu(14);
   
    % define persistent variables
    persistent alpha  % constant for low pass filter - only compute once
    persistent alpha1  % constant for low pass filter - only compute once
    persistent lpf_gyro_x   % low pass filter of x-gyro
    persistent lpf_gyro_y   % low pass filter of y-gyro
    persistent lpf_gyro_z   % low pass filter of z-gyro
    persistent lpf_static   % low pass filter of static pressure sensor
    persistent lpf_diff     % low pass filter of diff pressure sensor
    persistent lpf_accel_x  % low pass filter of x-accelerometer
    persistent lpf_accel_y  % low pass filter of y-accelerometer
    persistent lpf_accel_z  % low pass filter of z-accelerometer
    persistent xhat_a       % estimate of roll and pitch
    persistent P_a          % error covariance for roll and pitch angles
    persistent xhat_p       % estimate of pn, pe, Vg, chi, wn, we, psi
    persistent P_p          % error covariance for pn, pe, Vg, chi, wn, we, psi
    persistent y_gps_n_old  % last measurement of gps_n - used to detect new GPS signal
    persistent y_gps_e_old  % last measurement of gps_e - used to detect new GPS signal
    persistent y_gps_Vg_old % last measurement of gps_Vg - used to detect new GPS signal
    persistent y_gps_course_old  % last measurement of gps_course - used to detect new GPS signal
    
    
    
    % initialize persistent variables
    lpf_a = 50;
    lpf_a1 = 50;
    if t==0,
        alpha = exp(-lpf_a*P.Ts);
        alpha1 = exp(-lpf_a1*P.Ts);
        lpf_gyro_x   = 0;
        lpf_gyro_y   = 0;
        lpf_gyro_z   = 0;
        lpf_static   = P.rho*P.gravity*(-P.pd0);
        lpf_diff     = 1/2*P.rho*P.Va^2;
        lpf_accel_x  = 0;
        lpf_accel_y  = 0;
        lpf_accel_z  = 0;        
        xhat_a       = [0; 0];
        P_a          = (20*pi/180)^2*[1, 0; 0, 1];
        xhat_p       = [P.pn0; P.pe0; P.Va0; P.psi0; 0; 0; P.psi0];
        P_p          = diag([.03, .03, .1^2, (5*pi/180), .2^2, .22^2, (5*pi/180)]);
        y_gps_n_old  = -9999;
        y_gps_e_old  = -9999;
        y_gps_Vg_old = -9999;
        y_gps_course_old  = -9999;
    end
    
    %------------------------------------------------------------------
    % low pass filter gyros to estimate angular rates
    lpf_gyro_x = alpha*lpf_gyro_x + (1-alpha)*y_gyro_x;
    lpf_gyro_y = alpha*lpf_gyro_y + (1-alpha)*y_gyro_y;
    lpf_gyro_z = alpha*lpf_gyro_z + (1-alpha)*y_gyro_z;
    phat = lpf_gyro_x;
    qhat = lpf_gyro_y;
    rhat = lpf_gyro_z;
    
    %------------------------------------------------------------------
    % low pass filter static pressure sensor and invert to estimate
    % altitude
    lpf_static = alpha1*lpf_static + (1-alpha1)*y_static_pres;
    hhat = lpf_static/P.rho/P.gravity;
    
    % low pass filter diff pressure sensor and invert to estimate Va
    lpf_diff = alpha1*lpf_diff + (1-alpha1)*y_diff_pres;
    Vahat = sqrt(2/P.rho*lpf_diff);
    
    %------------------------------------------------------------------
    % low pass filter accelerometers
    lpf_accel_x = alpha*lpf_accel_x + (1-alpha)*y_accel_x;
    lpf_accel_y = alpha*lpf_accel_y + (1-alpha)*y_accel_y;
    lpf_accel_z = alpha*lpf_accel_z + (1-alpha)*y_accel_z;
    
    % invert accels to estimate phi and theta
    if abs(lpf_accel_z)<.001,
        phihat_accel = 0;
    else
        phihat_accel = atan(lpf_accel_y/lpf_accel_z);
    end
    thetahat_accel = asin(lpf_accel_x/P.gravity);
        
    %-------------------------------------------------------------------
    % implement continous-discrete EKF to estimate roll and pitch angles
    Q_a = 1*[0.000001, 0; 0, 0.000000001];
    Qinput = 1.0*diag([...
        P.sigma_gyro^2,...   % p
        P.sigma_gyro^2,...   % q
        P.sigma_gyro^2,...   % r
        ]);
    R_accel = P.sigma_accel^2;
    % R_accel = 0.0025*9.8;
    
    N = 10;
    % prediction step
    for i=1:N,
        cp = cos(xhat_a(1));  % cos(phi)
        sp = sin(xhat_a(1));  % sin(phi)
        tt = tan(xhat_a(2));  % tan(theta)
        ct = cos(xhat_a(2));  % cos(theta)
        f_a = [...
            phat + (qhat*sp+rhat*cp)*tt;...
            qhat*cp-rhat*sp;...
            ];
        A_a = [...
            (qhat*cp-rhat*sp)*tt, (qhat*sp+rhat*cp)/ct/ct;...
            -qhat*sp-rhat*cp, 0;...
            ];
        G_a = [ 1, sp*tt, cp*tt; ...
                0, cp, -sp];

        xhat_a = xhat_a + (P.Ts/N)*f_a;
        % P_a = P_a + (P.Ts/N)*(A_a*P_a + P_a*A_a' + Q_a +
        % G_a*Qinput*G_a'); % Tried this but it didn't work so well.
        P_a = P_a + (P.Ts/N)*(A_a*P_a + P_a*A_a' + Q_a);
    end
    % measurement updates
    cp = cos(xhat_a(1));  % cos(phi)
    sp = sin(xhat_a(1));  % sin(phi)
    ct = cos(xhat_a(2));  % cos(theta)
    st = sin(xhat_a(2));  % sin(theta)
    % x-axis accelerometer
    h_a = qhat*Vahat*st+P.gravity*st;
    C_a = [0, qhat*Vahat*ct+P.gravity*ct];
    L_a = P_a*C_a'/(R_accel+C_a*P_a*C_a');
    P_a = (eye(2)-L_a*C_a)*P_a;
    xhat_a = xhat_a + L_a*(y_accel_x - h_a);
    % y-axis accelerometer
    h_a = rhat*Vahat*ct-phat*Vahat*st-P.gravity*ct*sp;
    C_a = [-P.gravity*cp*ct, -rhat*Vahat*st-phat*Vahat*ct+P.gravity*st*sp];
    L_a = P_a*C_a'/(R_accel+C_a*P_a*C_a');
    P_a = (eye(2)-L_a*C_a)*P_a;
    xhat_a = xhat_a + L_a*(y_accel_y - h_a);
    % z-axis accelerometer
    h_a = -qhat*Vahat*ct-P.gravity*ct*cp;
    C_a = [P.gravity*sp*ct, (qhat*Vahat+P.gravity*cp)*st];
    L_a = P_a*C_a'/(R_accel+C_a*P_a*C_a');
    P_a = (eye(2)-L_a*C_a)*P_a;
    xhat_a = xhat_a + L_a*(y_accel_z - h_a);
    
     
    phihat   = xhat_a(1);
    thetahat = xhat_a(2);
    
        %-------------------------------------------------------------------
    % implement continous-discrete EKF to estimate pn, pe, chi, Vg
    Q_p = diag([...
        .0001,...  % pn
        .0001,...  % pe
        .0001,...  % Vg
        .000001,...  % chi
        .0001,...  % wn
        .0001,...  % we
        .0001,...  % psi
        ]);
    R_p = diag([...
        P.sigma_n_gps^2,...      % y_gps_n
        P.sigma_e_gps^2,...      % y_gps_e
        P.sigma_Vg_gps^2,...     % y_gps_Vg
        P.sigma_course_gps^2,... % y_gps_course
        0.001,...              % pseudo measurement #1
        0.001,...              % pseudo measurement #2
        ]);
    
    N = 10;
    % prediction step
    for i=1:N,
        psidot = (qhat*sin(phihat)+rhat*cos(phihat))/cos(thetahat);
        tmp = -psidot*Vahat*(xhat_p(5)*cos(xhat_p(7))+xhat_p(6)*sin(xhat_p(7)))/xhat_p(3);
        Vgdot = ((Vahat*cos(xhat_p(7))+xhat_p(5))*(-psidot*Vahat*sin(xhat_p(7)))...
            + (Vahat*sin(xhat_p(7))+xhat_p(6))*(psidot*Vahat*cos(xhat_p(7))))/xhat_p(3);
        f_p = [...
            xhat_p(3)*cos(xhat_p(4));...
            xhat_p(3)*sin(xhat_p(4));...
            Vgdot;...
            P.gravity/xhat_p(3)*tan(phihat)*cos(xhat_p(4)-xhat_p(7));...
            0;...
            0;...
            psidot;...
            ];
        A_p = [...
            0, 0,  cos(xhat_p(4)), -xhat_p(3)*sin(xhat_p(4)), 0, 0, 0;...
            0, 0,  sin(xhat_p(4)),  xhat_p(3)*cos(xhat_p(4)), 0, 0, 0;...
            0, 0,  -Vgdot/xhat_p(3), 0, -psidot*Vahat*sin(xhat_p(7))/xhat_p(3), psidot*Vahat*cos(xhat_p(7))/xhat_p(3), tmp;...
            0, 0,  -P.gravity/(xhat_p(3)^2)*tan(phihat)*cos(xhat_p(4)-xhat_p(7)), -P.gravity/xhat_p(3)*tan(phihat)*sin(xhat_p(4)-xhat_p(7)), 0, 0, P.gravity/xhat_p(3)*tan(phihat)*sin(xhat_p(4)-xhat_p(7));...
            0, 0, 0, 0, 0, 0, 0;...
            0, 0, 0, 0, 0, 0, 0;...
            0, 0, 0, 0, 0, 0, 0;...
            ];
        xhat_p = xhat_p + (P.Ts/N)*f_p;
        P_p = P_p + (P.Ts/N)*(A_p*P_p + P_p*A_p' + Q_p);
    end
    
    % measurement updates
    if   (y_gps_n~=y_gps_n_old)...
        |(y_gps_e~=y_gps_e_old)...
        |(y_gps_Vg~=y_gps_Vg_old)...
        |(y_gps_course~=y_gps_course_old),
    
        % gps North position
        h_p = xhat_p(1);
        C_p = [1, 0, 0, 0, 0, 0, 0];
        L_p = P_p*C_p'/(R_p(1,1)+C_p*P_p*C_p');
        P_p = (eye(7)-L_p*C_p)*P_p;
        xhat_p = xhat_p + L_p*(y_gps_n - h_p);
        % gps East position
        h_p = xhat_p(2);
        C_p = [0, 1, 0, 0, 0, 0, 0];
        L_p = P_p*C_p'/(R_p(2,2)+C_p*P_p*C_p');
        P_p = (eye(7)-L_p*C_p)*P_p;
        xhat_p = xhat_p + L_p*(y_gps_e - h_p);  
        % gps ground speed
        h_p = xhat_p(3);
        C_p = [0, 0, 1, 0, 0, 0, 0];
        L_p = P_p*C_p'/(R_p(3,3)+C_p*P_p*C_p');
        P_p = (eye(7)-L_p*C_p)*P_p;
        xhat_p = xhat_p + L_p*(y_gps_Vg - h_p);  
        % gps course
        % wrap course measurement
        while (y_gps_course - xhat_p(4))>pi, y_gps_course = y_gps_course - 2*pi; end
        while (y_gps_course - xhat_p(4))<-pi, y_gps_course = y_gps_course + 2*pi; end
        h_p = xhat_p(4);
        C_p = [0, 0, 0, 1, 0, 0, 0];
        L_p = P_p*C_p'/(R_p(4,4)+C_p*P_p*C_p');
        P_p = (eye(7)-L_p*C_p)*P_p;
        xhat_p = xhat_p + L_p*(y_gps_course - h_p);  
        % pseudo measurement #1 y_1 = Va*cos(psi)+wn-Vg*cos(chi)
        h_p = Vahat*cos(xhat_p(7))+xhat_p(5)-xhat_p(3)*cos(xhat_p(4));  % pseudo measurement
        C_p = [...
            0,... 
            0,... 
            -cos(xhat_p(4)),...
            xhat_p(3)*sin(xhat_p(4)),...
            1,...
            0,...
            -Vahat*sin(xhat_p(7)),...
            ];
        L_p = P_p*C_p'/(R_p(5,5)+C_p*P_p*C_p');
        P_p = (eye(7)-L_p*C_p)*P_p;
        xhat_p = xhat_p + L_p*(0 - h_p);
        % pseudo measurement #2 y_2 = Va*sin(psi) + we - Vg*sin(chi)
        h_p = Vahat*sin(xhat_p(7))+xhat_p(6)-xhat_p(3)*sin(xhat_p(4));  % pseudo measurement
        C_p = [...
            0,... 
            0,... 
            -sin(xhat_p(4)),...
            -xhat_p(3)*cos(xhat_p(4)),...
            0,...
            1,...
            Vahat*cos(xhat_p(7)),...
            ];
        L_p = P_p*C_p'/(R_p(6,6)+C_p*P_p*C_p');
        P_p = (eye(7)-L_p*C_p)*P_p;
        xhat_p = xhat_p + L_p*(0 - h_p);

        % update stored GPS signals
        y_gps_n_old      = y_gps_n;
        y_gps_e_old      = y_gps_e;
        y_gps_Vg_old     = y_gps_Vg;
        y_gps_course_old = y_gps_course;
    end
     
    pnhat    = xhat_p(1);
    pehat    = xhat_p(2);
    Vghat    = xhat_p(3);
    chihat   = xhat_p(4); 
    wnhat    = xhat_p(5);
    wehat    = xhat_p(6);
    psihat   = xhat_p(7);
  
    % not estimating these states 
    alphahat = 0;
    betahat  = 0;
    bxhat    = 0;
    byhat    = 0;
    bzhat    = 0;
    
      xhat = [...
        pnhat;...
        pehat;...
        hhat;...
        Vahat;...
        alphahat;...
        betahat;...
        phihat;...
        thetahat;...
        chihat;...
        phat;...
        qhat;...
        rhat;...
        Vghat;...
        wnhat;...
        wehat;...
        psihat;...
        ];
end