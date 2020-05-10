%-----------------------------------------------------------------------
%|*********************************************************************|
%|*        Define Controller for 20 Story 5 Bay Building Model        *|
%|*                Example is Active Tendon Control                   *|
%|*                                                                   *|
%|*                  University of Notre Dame                         *|
%|*                       November, 1999                              *|
%|*                                                                   *|
%|*               Coded by      Y.Ohtori                              *|
%|*                             R.E.Christenson                       *|
%|*               Supervised by Prof. B.F.Spencer, Jr.                *|
%|*********************************************************************|
%-----------------------------------------------------------------------

% -------------------------------------------
% --- Load Reduced Order Structural Model --- 
% -------------------------------------------
load red_mod			% for the design model, unlike the evaluation
				% model, the responses are given as all 
				% displacements, all velocities, and all abs.
				% accelerations.

% ------------------------------------
% --- Define Controller Parameters ---
% ------------------------------------
nms 	= 20; 				% number of floors 
nout 	= 3*nms;			% total number of model outputs
nst	= length(Ared(:,1));		% number of design model states
flr 	= [1:20];			% identify floors with actuators
ndev_flr= [4 2 2 ones(1,17)];	 	% number of actuators on each floor
nflr	= length(flr);			% number of floors with actuators
ndev	= sum(ndev_flr);		% total number of actuators
q1 	= 3e9;				% weighting on accel responses
r1	= 1/16;				% weighting on control forces

% -------------------------
% --- Design Controller ---
% -------------------------
% MODIFY STATE SPACE B & D MATRICES SUCH THAT INPUT IS DEVICE FORCE
K1      = diag(ndev_flr);	% multiple actuators per floor
K2	= eye(20);		% each actuator is connected between 2 floors
K2(1:19,2:20) = K2(1:19,2:20)-eye(19);	% force is = and opposite

Bmod    = Bred(:,2:nflr+1)*K2*K1;
Dmod    = Dred(:,2:nflr+1)*K2*K1;

% DEFINE WEIGHTING MATRICES
zvec	= [41:60];
Q 	= q1*eye(length(zvec));
R 	= diag([r1*ones(1,nflr)]);

% STATE FEEDBACK GAINS
Cout	= [Cred(zvec,:)];
Dout	= [Dmod(zvec,:)];
Kgain	= lqry(Ared,Bmod,Cout,Dout,Q,R);

% STATE ESTIMATOR
fbvec	= [4 8 12 16 20]+40;	% hor accel at floors 5, 9, 13, 17 and roof
numfb	= length(fbvec);
SW	= 25;
SV	= eye(numfb);
Cfb	= [Cred(fbvec,:)];
Dfb 	= [Dred(fbvec,:)];
Lgain   = lqe2(Ared,Bred(:,1),Cfb,SW,SV+Dfb(:,1)*SW*[Dfb(:,1)]',SW*[Dfb(:,1)]');

% FORM CONTINUOUS STATE-SPACE CONTROLLER
Ac	= Ared-Bmod*Kgain-Lgain*Cfb+Lgain*Dfb(:,2:nflr+1)*Kgain;
Bc	= Lgain;
Cc	= -Kgain;
Dc	= zeros(nflr,numfb);

% -----------------------
% --- Check Stability ---
% -----------------------
[Acl,Bcl,Ccl,Dcl]=feedback(Ared,[Bred(:,1) Bmod],Cred,[Dred(:,1) Dmod],...
                           Ac,Bc,Cc,Dc,[2:nflr+1],[fbvec]);

% CHECK STABILITY OF CONTROLLER
if (max(real(eig(Ac))) > 0) 
	disp('UNSTABLE CONTROLLER !')
	return
end

% CHECK STABILITY OF CLOSED LOOP SYSTEM (w/CONTROL DESIGN MODEL)
if (max(real(eig(Acl))) > 0) 
	disp('CLOSED LOOP UNSTABLE !')
	return
end

% --------------------------------
% --- Convert to Discrete Form ---
% --------------------------------
Tcd=0.01;
[Acd,Bcd,Ccd,Dcd] = c2dm(Ac,Bc,Cc,Dc,Tcd,'zoh');

% CHECK STABILITY OF DISCRETE CONTROLLER
EIGVAL = eig(Acd);
if (max(abs(EIGVAL)) > 1) 
	disp('UNSTABLE CONTROLLER (DISCRETE)!')
	return
end

 Max_frc = 1000e3;		% Maximum Force of Devices
 Max_voltd = 10;		% Maximum Voltage of Control Signal
 gain_ctr = Max_frc/Max_voltd;	% Gain

% -----------------------------------------------------
% --- Data for Sensor (for Structure: acceleration) ---
% -----------------------------------------------------
 Num_decimate = round( dt_out / dt_cal);
 Max_acc   = 9.81;              % Maximum Acceleration of Sensors
 Max_volts  = 10;                % Maximum Voltage of Sensors
 gain_snr = Max_volts/Max_acc;  % Gain                              
 Ds        = gain_snr*eye(5);  % Sensor Model (Gain Matrix)

 NP_snr   = 0.03/1000*(10);    % Noise Power for Measurement Noise 
 ST_snr   = 0.01;              % Sample Time for Measurement Noise 

% find sufficient white noise seeds
 numv	  = numfb;
rand('seed',123);
Seed_snr  = round(100000*rand(numv,1));

% ---------------------------------------
% --- Data for A/D and D/A Converters ---
% ---------------------------------------

 SL_ctr   = Max_voltd;          % Saturation Level                  
 QI_ctr   = 2*Max_voltd/(2^16);  % Quantization Interval             

 SL_snr   = Max_volts;          % Saturation Level                  
 QI_snr   = 2*Max_volts/(2^16);  % Quantization Interval             

 Bcd 	= Bcd/gain_snr;		% Adjust controller matrices accordingly
 Ccd 	= Ccd/gain_ctr;
 Dcd 	= Dcd/gain_snr/gain_ctr;

% -------------------------------
% --- Data for Ideal Actuator ---
% -------------------------------

K_load(1:2:39,:) = eye(20);		% load force vector on the building
K_load(2:2:39,2:20) = -eye(19);
K_mult = diag(ndev_flr);		% multiple actuators per level

Kf      = K_load*K_mult*gain_ctr;	% output force gain matrix
