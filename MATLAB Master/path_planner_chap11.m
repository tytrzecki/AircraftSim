% path planner
%
% Modified:  
%   - 4/06/2010 - RWB
%
% output is a vector containing P.num_waypoints waypoints
%
% input is the map of the environment
function out = path_planner_chap11(in,P,map)

  NN = 0;
  % pn        = in(1+NN);
  % pe        = in(2+NN);
  % h         = in(3+NN);
  % Va      = in(4+NN);
  % alpha   = in(5+NN);
  % beta    = in(6+NN);
  % phi     = in(7+NN);
  % theta   = in(8+NN);
  % chi     = in(9+NN);
  % p       = in(10+NN);
  % q       = in(11+NN);
  % r       = in(12+NN);
  % Vg      = in(13+NN);
  % wn      = in(14+NN);
  % we      = in(15+NN);
  % psi     = in(16+NN);
  % flag_new_waypoints =  in(17+NN);
  NN = NN + 17;
  % t         =  in(1+NN);


  num_waypoints = 4;
  % format for each point is [pn, pe, pd, chi, Va^d] where the position
  % of the waypoint is (pn, pe, pd), the desired course at the waypoint
  % is chi, and the desired airspeed between waypoints is Va
  % if chi!=-9999, then Dubins paths will be used between waypoints.
  if 0,%CHANGEABLE: 1 means Fillet, 0 means Dubins
 chix=-9999;
 chix2=chix;
 chix3=chix;
 sgn=0;
  else  % Dubins
      chix=45*pi/180;
      chix2=0;
      chix3=chix*-3;
      sgn=-50;

  end
    wpp= [...
             0,   0,   -50, chix2, P.Va0;...
             440, 20,   -15+sgn, chix, P.Va0;...
             60,   350, -15+sgn, chix, P.Va0;...
            420, 400, -50, chix3, P.Va0;...
%             50, 50, -250, chix, P.Va0;...
%             0,   0,   -250, chix, P.Va0;...
           ];
%  end      
  for i=5:P.size_waypoint_array,
      wpp = [...
          wpp;...
          -9999, -9999, -9999, -9999, -9999;...
          ];
  end
  
  out = [num_waypoints; reshape(wpp', 5*P.size_waypoint_array, 1)]; 

end