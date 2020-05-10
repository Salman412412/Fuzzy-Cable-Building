%-----------------------------------------------------------------------
%|*********************************************************************|
%|*      Simulation Control Program for Nonlinear Benchmark Study     *|
%|*                                                                   *|
%|*                  University of Notre Dame                         *|
%|*                       November, 1999                              *|
%|*                                                                   *|
%|*               Coded by      Y.Ohtori                              *|
%|*                             R.E.Christenson                       *|
%|*               Edited by     Salman Masoumi(salman412412@gmail.com)*|
%|*********************************************************************|
%-----------------------------------------------------------------------

clear

%---------------------------------------
% --- Define Building for Simulation ---
%---------------------------------------
No_bld     = 20;			% 20 Story Building

%---------------------------------------------
% --- Build Benchmark Model for Simulation ---
%---------------------------------------------
sprintf(' ---> BUILDING MODEL FOR SIMULATION')
Bld_NLBM

%----------------------------------------------
% --- Simulate El Centro Earthquake Records ---
%----------------------------------------------
OPTIONS = simset('solver','ode5','FixedStep',dt_cal);
isim = 0;
tf       = 100;            % Duration (sec)

% Load the best Genes which was obtained by TrainFuzzy.m
load('TheBest.mat');
global x1;
x1 = zeros(tf*100+1, 96);
x1(:,1) = (0:0.01:tf)';
for i = 1:(tf*100+1)
    x1(i,2:96) = (bestGenes);
end


for intensity=0.5:0.5:1.5
    isim = isim+1;
    
    sprintf([' ---> SIMULATION: ElCentro ' num2str(intensity,2)])
    sim('Sim_NLBM_2016',[0 tf],OPTIONS,[])
    load Mem_damage.out -ascii
    
    save(['BLD' num2str(No_bld) '_' num2str(isim) '.mat'],...
        'No_bld','yf','ye','Mem_damage','tf')
end

%----------------------------------------------
% --- Simulate Hachinohe Earthquake Records ---
%----------------------------------------------
delete_line('Sim_NLBM_2016','Ground Accel/1','intensity/1')
add_line   ('Sim_NLBM_2016','Ground Accel/2','intensity/1')
tf       = 100;            % Duration (sec)
for intensity=0.5:0.5:1.5
    isim = isim+1;
    
    sprintf([' ---> SIMULATION: Hachinohe ' num2str(intensity,2)])
    sim('Sim_NLBM_2016',[0 tf],OPTIONS,[])
    load Mem_damage.out -ascii
    
    save(['BLD' num2str(No_bld) '_' num2str(isim) '.mat'],...
        'No_bld','yf','ye','Mem_damage','tf')
end

%-----------------------------------------------
% --- Simulate Northridge Earthquake Records ---
%-----------------------------------------------
delete_line('Sim_NLBM_2016','Ground Accel/2','intensity/1')
add_line   ('Sim_NLBM_2016','Ground Accel/3','intensity/1')
tf       = 100;            % Duration (sec)
for intensity=0.5:0.5:1.0
    isim = isim+1;
    
    sprintf([' ---> SIMULATION: Northridge ' num2str(intensity,2)])
    sim('Sim_NLBM_2016',[0 tf],OPTIONS,[])
    load Mem_damage.out -ascii
    
    save(['BLD' num2str(No_bld) '_' num2str(isim) '.mat'],...
        'No_bld','yf','ye','Mem_damage','tf')
end

%-----------------------------------------
% --- Simulate Kobe Earthquake Records ---
%-----------------------------------------
delete_line('Sim_NLBM_2016','Ground Accel/3','intensity/1')
add_line   ('Sim_NLBM_2016','Ground Accel/4','intensity/1')
tf       = 180;            % Duration (sec)
for intensity=0.5:0.5:1.0
    isim = isim+1;
    
    sprintf([' ---> SIMULATION: Kobe ' num2str(intensity,2)])
    sim('Sim_NLBM_2016',[0 tf],OPTIONS,[])
    load Mem_damage.out -ascii
    
    save(['BLD' num2str(No_bld) '_' num2str(isim) '.mat'],...
        'No_bld','yf','ye','Mem_damage','tf')
end

%------------------------------------------------------
%--- RETURN SIMULATION TO ORIGINAL CONFIGURATION ------
%------------------------------------------------------
delete_line('Sim_NLBM_2016','Ground Accel/4','intensity/1')
add_line   ('Sim_NLBM_2016','Ground Accel/1','intensity/1')

EvalNLBM;