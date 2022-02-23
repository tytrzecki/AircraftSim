function mav_dynamics_L2(block)
% Level-2 MATLAB file S-Function for Chapter 3 HW solution -- mavsim_chap3.mdl.

  setup(block);
  
%endfunction

function setup(block)
  
  %% Register number of dialog parameters   
  block.NumDialogPrms = 17;

  %% Register number of input and output ports
  block.NumInputPorts  = 1;
  block.NumOutputPorts = 1;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
 
  block.InputPort(1).Dimensions        = 6;
  block.InputPort(1).DirectFeedthrough = false;
  
  block.OutputPort(1).Dimensions       = 12;
  
  %% Set block sample time to continuous
  block.SampleTimes = [0 0];
  
  %% Setup Dwork
  block.NumContStates = 12;

  %% Set the block simStateCompliance to default (i.e., same as a built-in block)
  block.SimStateCompliance = 'DefaultSimState';

  %% Register methods
  block.RegBlockMethod('InitializeConditions',    @InitConditions);  
  block.RegBlockMethod('Outputs',                 @Output);  
  block.RegBlockMethod('Derivatives',             @Derivative);  
  
%endfunction

function InitConditions(block)

  %% Initialize Dwork
  %% Initial conditions on states
  block.ContStates.Data(1)  = block.DialogPrm(6).Data;  % pn0
  block.ContStates.Data(2)  = block.DialogPrm(7).Data;  % pe0
  block.ContStates.Data(3)  = block.DialogPrm(8).Data;  % pd0
  block.ContStates.Data(4)  = block.DialogPrm(9).Data;  % u0
  block.ContStates.Data(5)  = block.DialogPrm(10).Data; % v0
  block.ContStates.Data(6)  = block.DialogPrm(11).Data; % w0
  block.ContStates.Data(7)  = block.DialogPrm(12).Data; % phi0
  block.ContStates.Data(8)  = block.DialogPrm(13).Data; % theta0
  block.ContStates.Data(9)  = block.DialogPrm(14).Data; % psi0
  block.ContStates.Data(10) = block.DialogPrm(15).Data; % p0
  block.ContStates.Data(11) = block.DialogPrm(16).Data; % q0
  block.ContStates.Data(12) = block.DialogPrm(17).Data; % r0 

%endfunction

function Output(block)

  block.OutputPort(1).Data = block.ContStates.Data;
  
%endfunction

function Derivative(block)

    mass = block.DialogPrm(1).Data;
    Jx = block.DialogPrm(2).Data;
    Jy = block.DialogPrm(3).Data;
    Jz = block.DialogPrm(4).Data;
    Jxz = block.DialogPrm(5).Data;
    
    x = block.ContStates.Data;
    uu = block.InputPort(1).Data;

    pn    = x(1);
    pe    = x(2);
    pe    = x(3);
    u     = x(4);
    v     = x(5);
    w     = x(6);
    phi   = x(7);
    theta = x(8);
    psi   = x(9);
    p     = x(10);
    q     = x(11);
    r     = x(12);
    fx    = uu(1);
    fy    = uu(2);
    fz    = uu(3);
    ell   = uu(4);
    m     = uu(5);
    n     = uu(6);
    
    pdot = [cos(theta)*cos(psi),  sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi), cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);...
            cos(theta)*sin(psi),  sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi), cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);...
            -sin(theta),          sin(phi)*cos(theta),                            cos(phi)*cos(theta)]    *     [u; v; w];
    pndot = pdot(1);
    pedot = pdot(2);
    pddot = pdot(3);
    
    veldot = [r*v - q*w;...
              p*w - r*u;...
              q*u - p*v] + (1/mass)*[fx; fy; fz];
    udot = veldot(1);
    vdot = veldot(2); 
    wdot = veldot(3);
    
    angledot = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);...
                0, cos(phi),            -sin(phi);...
                0, sin(phi)/cos(theta), cos(phi)/cos(theta)]   *    [p; q; r];
    phidot = angledot(1);
    thetadot = angledot(2);
    psidot = angledot(3);
    
    Gamma = Jx*Jz-Jxz^2;
    Gamma1 = (Jxz*(Jx-Jy+Jz))/Gamma;
    Gamma2 = (Jz*(Jz-Jy)+Jxz^2)/Gamma;
    Gamma3 = Jz/Gamma;
    Gamma4 = Jxz/Gamma;
    Gamma5 = (Jz-Jx)/Jy;
    Gamma6 = Jxz/Jy;
    Gamma7 = ((Jx-Jy)*Jx+Jxz^2)/Gamma;
    Gamma8 = Jx/Gamma;
    
    spindot = [Gamma1*p*q-Gamma2*q*r;... 
               Gamma5*p*r-Gamma6*(p^2-r^2);...
               Gamma7*p*q-Gamma1*q*r]...
               + ...
              [Gamma3*ell+Gamma4*n; 
               m/Jy; 
               Gamma4*ell+Gamma8*n];
    pdot = spindot(1);
    qdot = spindot(2);
    rdot = spindot(3);    

    block.Derivatives.Data = [pndot; pedot; pddot; udot; vdot; wdot; phidot; thetadot; psidot; pdot; qdot; rdot];
     
%endfunction

