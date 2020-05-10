%-----------------------------------------------------------------------
%|*********************************************************************|
%|*      Calculate Evaluation Criteria for the 20 Story Building      *|
%|*                                                                   *|
%|*                  University of Notre Dame                         *|
%|*                       November, 1999                              *|
%|*                                                                   *|
%|*               Coded by      Y.Ohtori                              *|
%|*                             R.E.Christenson                       *|
%|*               Supervised by B.F.Spencer, Jr.  	              *|
%|*********************************************************************|
%-----------------------------------------------------------------------

clear

% --------------------------------------------------------------
% --- Pre-Define some Building Data Specific to the 20 Story ---
% ---  Building Controlled in the Sample Design              ---
% --------------------------------------------------------------
hi = [ 5.4864 3.9624 3.9624 3.9624 3.9624 3.9624 3.9624 3.9624 ...
       3.9624 3.9624 3.9624 3.9624 3.9624 3.9624 3.9624 3.9624 ...
       3.9624 3.9624 3.9624 3.9624 ];
mi = [ 563320 551640 551640 551640 551640 551640 551640 551640 551640 ...
       551640 551640 551640 551640 551640 551640 551640 551640 551640 ...
       551640 583760];
W = sum(mi*9.81);
No_bld  = 20;
ndev_flr= [4 2 2 ones(1,17)];
dt	= 0.01;
load unctrl20

% --------------------------------------------------------
% --- Find the Evaluation Criteria for Each Earthquake ---
% --------------------------------------------------------

for isim=1:10
 load(['BLD' num2str(No_bld) '_' num2str(isim) '.mat'])

idx_acc = [1:3:size(ye,2)];
idx_vel = [2:3:size(ye,2)];
idx_disp= [3:3:size(ye,2)];

% -----------------------------------------------------
% --- Building Response Evaluation Criteria (J1-J6) ---
% -----------------------------------------------------

di   = [ye(:,idx_disp(1)) diff(ye(:,idx_disp)')'];
xddi = ye(:,idx_acc);

% Peak Interstory Drift [J1]

J(1,isim) = max(max(abs(di))./hi) / d_max(isim);

% Peak Floor Acceleration [J2]

J(2,isim) = max(max(abs(xddi))) / xdd_max(isim);

% Peak Base Shear [J3]

J(3,isim) = max(abs(xddi*mi')) / F_max(isim);

% Normed Based Interstory Drift [J4]

J(4,isim) = max(sqrt((1/tf)*sum(di.*di)*dt)./hi) / d_max_norm(isim);

% Normed Based Floor Acceleration [J5]

J(5,isim) = max(sqrt((1/tf)*sum(xddi.*xddi)*dt)) / xdd_max_norm(isim);

% Normed Based Base Shear [J6]

J(6,isim) = sqrt((1/tf)*sum((xddi*mi').*(xddi*mi'))*dt) / F_max_norm(isim);


% ----------------------------------------------------
% --- Building Damage Evaluation Criteria (J7-J10) ---
% ----------------------------------------------------

phii      = Mem_damage(:,5);
phii_norm = Mem_damage(:,6);
dEi       = Mem_damage(:,7);

% Ductility Factor [J7]

J(7,isim) = max(abs(phii)) / phi_max(isim);

% Dissipated Energy of the Curvatures at the End of Members [J8]

if N_d(isim)>0
 J(8,isim) = max(abs(dEi)) / E_max(isim);
else
 J(8,isim) = 999999;
end

% Ratio of Plastic Connections

if N_d(isim)>0
 J(9,isim) = size(find(Mem_damage(:,5)>1.0),1)/N_d(isim);
else
 J(9,isim) = 999999;
end

% Normed Basis Ductility Factor

J(10,isim) = max(abs(phii_norm)) / phi_max_norm(isim);


% ----------------------------------------------------
% --- Control Device Evaluation Criteria (J11-J14) ---
% ----------------------------------------------------

f  = yf(:,1:No_bld);
yi = [yf(:,No_bld+1) diff(yf(:,No_bld+1:2*No_bld),[],2)];
ydi= [yf(:,2*No_bld+1) diff(yf(:,2*No_bld+1:3*No_bld),[],2)];
P  = abs(ydi.*yf(:,1:No_bld));

% Peak Control Force

J(11,isim) = max(max(abs(f))) / W;

% Peak Control Device Stroke

J(12,isim) = max(max(abs(yi))) / x_max(isim);

% Peak Control Power

J(13,isim) = max(P*ndev_flr') / (xd_max(isim)*W);

% Total Control Power

J(14,isim) = (1/tf)*sum(P)*dt*ndev_flr' / (xd_max(isim)*W);

% Maximum Control Device Values

Max_y(isim) = max(max(abs(yi)));

Max_f(isim) = max(max(abs(f)));

% ------------------------------------------------------
% --- Control Strategy Evaluation Criteria (J15-J17) ---
% ------------------------------------------------------

% Number of Control Devices

J(15,isim) = sum(ndev_flr);

% Number of Required Sensors

J(16,isim) = 5;		% Obtained previously from design

% Computational Resourse

J(17,isim) = 20;			% Obtained previously from design

end
