%-----------------------------------------------------------------------
%|*********************************************************************|
%|*         User Defined Data for 20 Story 5 Bay Building Data        *|
%|*                                                                   *|
%|*                  University of Notre Dame                         *|
%|*                       November, 1999                              *|
%|*                                                                   *|
%|*               Coded by      Y.Ohtori                              *|
%|*                             R.E.Christenson                       *|
%|*               Supervised by B.F.Spencer, Jr.  	              *|
%|*********************************************************************|
%-----------------------------------------------------------------------

% ****************KEY*****************
% *Dir. = 1: Horizontal              *
% *       2: Vertical                *
% *       3: Rotation                *
% *                                  *
% *Resp. = 1: Absolute acceleration  *
% *        2: Relative velocity      *
% *        3: Relative displacement  *
% ************************************

% ------------------------------------------------
% --- Define File Name of Structural Data File ---
% ---      for the SIMULINK S-function         ---
% ------------------------------------------------ 
  fid10 = fopen('User_Def.data','w');

% ------------------------------------------------
% --- Observation Points for Evaluations: (ye) ---
% ------------------------------------------------
%  obs(i,1): No.
%  obs(i,2): Node No. (from Fig 1, 2 or 3 in paper)
%  obs(i,3): Direction (1, 2, or 3)
%  obs(i,4): Response (1, 2, or 3)
% ----------------------------------
   obs = 0;
   no =0;
   for No_node=19:6:133
       for comp = 1:3
           no = no + 1;
           obs(no,1) = no;
           obs(no,2) = No_node;
           obs(no,3) = 1;       % Horizontal Component
           obs(no,4) = comp;    % Required Response 
       end
    end

   Num_obs = size(obs,1);

% ---------------------------------------------
% --- Sensor Positions for Aquisition: (ym) ---
% ---------------------------------------------
%  snr(i,1): No.
%  snr(i,2): Node No. (from Fig 1, 2 or 3 in paper)
%  snr(i,3): Direction (1, 2, or 3)
%  snr(i,4): Response (1, 2, or 3)
% ----------------------------------
   snr = 0;
   no =0;
   for No_node=[37 61 85 109 133];	% on floors 4, 8, 12, 16 and 20
       no = no + 1;
       snr(no,1) = no;
       snr(no,2) = No_node;
       snr(no,3) = 1;       % Horizontal Component
       snr(no,4) = 1;       % Required Response 
    end

    Num_snr = size(snr,1);

% -----------------------------------------
% --- Connection Points of Devices (yc) ---
% -----------------------------------------
%  cps(i,1): No.
%  cps(i,2): Node No. (from Fig 1, 2 or 3 in paper)
%  cps(i,3): Direction (1, 2, or 3)
%  cps(i,4): Response (1, 2, or 3)
% ----------------------------------
   cps = 0;
   no =0;
   for No_node=19:6:133;
       no = no + 1;
       cps(no,1) = no;
       cps(no,2) = No_node;
       cps(no,3) = 1;       % Horizontal Direction
       cps(no,4) = 3;       % Required Response 
    end
   for No_node=19:6:133;
       no = no + 1;
       cps(no,1) = no;
       cps(no,2) = No_node;
       cps(no,3) = 1;       % Horizontal Direction
       cps(no,4) = 2;       % Required Response 
    end

    Num_cps = size(cps,1);

% ---------------------------------------
% --- Location of Control Forces: (f) ---
% ---------------------------------------
%  cf(i,1): No.
%  cf(i,2): Node No. (from Fig 1, 2 or 3 in paper)
%  cf(i,3): Direction (1, 2, or 3)
% ----------------------------------
   cf = zeros(39,3);
   no =0;
   for No_node=19:6:127;
       no = no + 1;
       cf(no,1) = no;
       cf(no,2) = No_node;
       cf(no,3) = 1;       % Horizontal Direction
       no = no + 1;
       cf(no,1) = no;
       cf(no,2) = No_node;
       cf(no,3) = 1;       % Horizontal Direction
    end
       cf(39,:) = [39 133 1];

    Num_cf = size(cf,1);














