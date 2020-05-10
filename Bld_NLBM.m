%-----------------------------------------------------------------------
%|*********************************************************************|
%|*     Build Nonlinear Benchmark Models (3-, 9- and 20- stories)     *|
%|*                                                                   *|
%|*                  University of Notre Dame                         *|
%|*                       November, 1999                              *|
%|*                                                                   *|
%|*               Coded by      Y.Ohtori                              *|
%|*                             R.E.Christenson                       *|
%|*               Supervised by B.F.Spencer, Jr.  	              *|
%|*********************************************************************|
%-----------------------------------------------------------------------

%----------------------------------------
% --- Define Building Characteristics ---
%----------------------------------------
 Idx_linear = 0;			% 1:Linear 0:Nonlinear
 Damp_type  = 3;			% 1:Mass, 2:Stiffness, 3:Rayleigh
 Idx_Band   = 1;			% Solver  1:Banded, 0: Symmetric

% ------------------------
% --- Load input files ---
% ------------------------

%--- LA 3 Story Building ---
  if No_bld == 3
     Str_file = 'struct3';
     IO_file  = 'inout_3';		% Participant must define
     Cnt_file= 'ctrlr_3';		% Participant must define
  end
%

%--- LA 9 Story Building ---
  if No_bld == 9
     Str_file = 'struct9';
     IO_file  = 'inout_9';		% Participant must define
     Cnt_file= 'ctrlr_9';		% Participant must define
  end
%
%--- LA 20 Story Building ---
 if No_bld == 20
    Str_file = 'struct20';
    IO_file  = 'inout_20';		% Sample Control Design
    Cnt_file= 'ctrlr_20';		% Sample Control Design
 end
%
%*********************************
     eval(Str_file);
     eval(IO_file);
% ------------------------------
% Calculate the nodal coordinates
% ------------------------------
     EPS = 1.0E-14;

     y0 = 0.0;
     i  = 0;
     for no_story=1:Num_story+1
        x0 = 0.0;
        for no_bay=1:Num_bay+1
           i  = i + 1;
           x(i) =  x0;
           y(i) = -y0;
           if no_bay<=Num_bay 
              x0 = x0 + width(no_bay);
           end
        end
        if no_story<=Num_story
           y0 = y0 + height(no_story);
        end
     end

% ----------------------------
% Set Node Numbers of Elements
% -------------------------------------------------------------------
%   element_tbl(:,1) = Node_a
%   element_tbl(:,2) = Node_b
%   element_tbl(:,3) = Material Table Number
%   element_tbl(:,4) = Type of Element (1: Column, 2: Beam, 3:Brace)
% -------------------------------------------------------------------
%  Judge Linear or Nonlinear Analysis
% ------------------------------------
  if Idx_linear == 1
     mat_tbl(:,7) = 10*mat_tbl(:,7);
     mat_tbl(:,8) = 10*mat_tbl(:,8);
  end
% -----------------------------------
% ----------------------
% Set Element Properties
% ----------------------

   for i=1:Num_elem
     na = element_tbl(i,1);
     nb = element_tbl(i,2);
     elem_prop(i,1) = sqrt((x(na)-x(nb))^2 + (y(na)-y(nb))^2);    % Length
     elem_prop(i,2) = mat_tbl(element_tbl(i,3),2);                % EI1
     elem_prop(i,3) = mat_tbl(element_tbl(i,3),3);                % EI2
     elem_prop(i,4) = mat_tbl(element_tbl(i,3),4);                % EI3
     elem_prop(i,5) = mat_tbl(element_tbl(i,3),5);                % EA
     elem_prop(i,6) = mat_tbl(element_tbl(i,3),6);                % GA
     elem_prop(i,7) = mat_tbl(element_tbl(i,3),7);                % d1
     elem_prop(i,8) = mat_tbl(element_tbl(i,3),8);                % d2
     elem_prop(i,9) = mat_tbl(element_tbl(i,3),9);                % Type
   end

% ------------------------------------
% Assemble Global Matricies [M] & [K]
% ------------------------------------

   MM = zeros(Num_DOF);
   KK = zeros(Num_DOF);
   CC = zeros(Num_DOF);

% --- Stiffness Matrix Fixed Fixed [K] ---

   for i=1:Num_elem
    na = element_tbl(i,1);
    nb = element_tbl(i,2);
    Ia = 3*(na-1)+[1 2 3];
    Ib = 3*(nb-1)+[1 2 3];
    EI = elem_prop(i,2);
    EA = elem_prop(i,5);
    invL=1.0/elem_prop(i,1);
    T=diag(ones(6,1));
    c(i)=(x(nb)-x(na))*invL;
    s(i)=(y(nb)-y(na))*invL;
    a=[1;2];
    b=a+3;
    T(a,a)=[c(i)  s(i);-s(i) c(i)];
    T(b,b)=[c(i)  s(i);-s(i) c(i)];

  if element_tbl(i,4)<=2			%      |-----|
    g=EA/EI;
    k1=[g;0;0;-g;0;0];
    k2=[0;12*invL*invL;6*invL;0;-12*invL*invL;6*invL];
    k3=[0;6*invL;4;0;-6*invL;2];
    k6=[0;6*invL;2;0;-6*invL;4];
    khat=(EI*invL)*[k1 k2 k3 -k1 -k2 k6];
    km=T'*khat*T;

  elseif element_tbl(i,4)==3			%      o-----o
    khat=zeros(6);
    k1=(EA*invL)*[1 -1;-1 1];
    khat([1 4],[1 4])=k1;
    km=T'*khat*T;

  elseif element_tbl(i,4)==4			%      o-----|
    khat=zeros(6);
    k1=(EA*invL)*[1 -1;-1 1];
    k2=(3*EI*invL)*[1*invL*invL -1*invL*invL 1*invL;-1*invL*invL 1*invL*invL -1*invL;1*invL -1*invL 1];
    khat([1 4],[1 4])=k1;
    khat([2 5 6],[2 5 6])=k2;
    km=T'*khat*T;

  elseif element_tbl(i,4)==5			%      |-----o
    khat=zeros(6);
    k1=(EA*invL)*[1 -1;-1 1];
    k2=(3*EI*invL)*[1*invL*invL 1*invL -1*invL*invL;1*invL 1 -1*invL;-1*invL*invL -1*invL 1*invL*invL];
    khat([1 4],[1 4])=k1;
    khat([2 3 5],[2 3 5])=k2;
    km=T'*khat*T;
  end

    KK([Ia Ib],[Ia Ib]) = KK([Ia Ib],[Ia Ib])+km;

  end

% --- Lumped Mass Matrix [M] ---
% only the beam elements have mass associated with them

    alpha = 1e-06;

   for i=1:Num_elem
    na = element_tbl(i,1);
    nb = element_tbl(i,2);
    Ia = 3*(na-1)+[1 2 3];
    Ib = 3*(nb-1)+[1 2 3];
    mi = Element_mass(i)/2*diag([1 1 alpha*(elem_prop(i,1)^2)/210 ...
                               1 1 alpha*(elem_prop(i,1)^2)/210]);
    MM([Ia Ib],[Ia Ib]) = MM([Ia Ib],[Ia Ib])+mi;
   end

% --- Place an arbitrary rotational mass on the basement pinned locations ---

   pinned_base=(Fix_node(:,1)<=Num_bay+1).*(1-Fix_node(:,4));
  
   for i = 1:length(pinned_base)
    if pinned_base(i)==1
     MM(Fix_node(i,1)*3,Fix_node(i,1)*3) = ...
     MM(Fix_node(i,1)*3+(Num_bay+1)*3,Fix_node(i,1)*3+(Num_bay+1)*3);
    end
   end

% ------------------------------------
% Eliminate Boundary Condition DOFs
% ------------------------------------

    Ref_BND = zeros(1,Num_DOF);

  for i=1:Num_BND
    Ref_BND(3*(Fix_node(i,1)-1) + [1 2 3]) = Fix_node(i,1+[1 2 3]);
  end

    Order = 1:Num_DOF;

    free_vec = Order(Ref_BND==0);

% --- Eliminate the Fixed Boundary DOF ---

    M_free  = MM(free_vec,free_vec);
    K_free  = KK(free_vec,free_vec);

    Ord_vec = [Order(Ref_BND==0) Order(Ref_BND==1)];
    Out_idx(Ord_vec) = Order;
    Num_free = length(free_vec);

% ------------------------------------
% Eliminate Slave Nodes
% ------------------------------------

   T_rigid = eye(Num_DOF);

   for i=1:Num_mst
     i_mst = 3*(slv_tbl(i,1)-1) + slv_tbl(i,2);
     for j=1:slv_tbl(i,3)
         i_slv = 3*(slv_tbl(i,3+j)-1) + slv_tbl(i,2);
         T_rigid(i_slv,i_slv) = 0;
         T_rigid(i_slv,i_mst) = 1;
     end
   end

   Num_slv = max(slv_tbl(:,3));			    % Max Number of Slave Nodes

   Tr_rigid = T_rigid(free_vec,free_vec);

   act_vec = find(diag(Tr_rigid));

   Tr2_rigid = Tr_rigid(:,act_vec);                   % Num_free * Num_act

   Num_act = length(act_vec);

   diagT = cumsum(diag(Tr_rigid));
   for i=1:Num_mst
     i_mst = 3*(slv_tbl(i,1)-1) + slv_tbl(i,2);
     for j=1:slv_tbl(i,3)
         i_slv = 3*(slv_tbl(i,3+j)-1) + slv_tbl(i,2);
         diagT(Out_idx(i_slv)) = diagT(Out_idx(i_mst));
     end
   end

% --- Remove Rigid Body -----------------

    M = Tr2_rigid' * M_free * Tr2_rigid;
    K = Tr2_rigid' * K_free * Tr2_rigid;

% ------------------------------------
% Assemble Damping Matrix [C]
% ------------------------------------

 invM = inv(M);

 [eig_vec,eig_val] = eig(invM*K);
 [omeg,w_order]    = sort(sqrt(diag(eig_val)));
 mode_vec = eig_vec(:,w_order);

 C_bar = zeros(Num_DOF);
 zeta1 = zeta_cr;
 
if Damp_type==1 
    zeta_vec = zeta1*omeg(1)./omeg;
    C_bar    = diag(2.0*zeta_vec.*omeg);

 elseif Damp_type==2
    zeta_vec = zeta1/omeg(1)*omeg;
    zeta_vec = min(zeta_vec,h_max);
    C_bar    = diag(2.0*zeta_vec.*omeg);

 elseif Damp_type==3
     for i=1:size(omeg,1)
         zeta_vec(i) = zeta1*(omeg(1)*omeg(nCutOff) + omeg(i)^2) ...
                       / (omeg(i)*(omeg(1)+omeg(nCutOff)));
     end

     C_bar    = diag(2*zeta_vec'.*omeg);

 else
 end

 C = M*mode_vec*C_bar*inv(mode_vec);

% ------------------------------------
% Determine Global Damping Matrix [C]
% ------------------------------------

 invTr2 = pinv(Tr2_rigid);
 invTr2t= pinv(Tr2_rigid');
 C_free = invTr2t * C * invTr2;
 for i=1:Num_DOF
    for j=1:Num_DOF
       if i<=Num_free & j<=Num_free
          CC(Ord_vec(i),Ord_vec(j)) = C_free(i,j);
       else
          CC(Ord_vec(i),Ord_vec(j)) = 0.0;
       end
    end
 end

% ------------------------------
% Output for SIMLINK S-function 
% ------------------------------

    fid1 = fopen(f1_name,'w');

% --- Search the Band width of the Stiffness Matrix ---

 Idx_nz = find(KK);           
 mat_size = size(KK);
 iv = rem(Idx_nz-1,mat_size(1))+1;
 jv = (Idx_nz -iv)/mat_size(1)+1;
 Band_width1 = max(abs(iv-jv+1));

 Idx_nz = find(K);           
 mat_size = size(K);
 iv = rem(Idx_nz-1,mat_size(1))+1;
 jv = (Idx_nz -iv)/mat_size(1)+1;
 Band_width2 = max(abs(iv-jv+1));
 
 Band_width = max(Band_width1,Band_width2);
 Band_width = min(Band_width,Num_act);

 if (Idx_Band ~= 1) 
    Band_width = Num_act;
 end

% --- Print Data to file -------------------------------------------------

 Num_plast = 0; 		% Number of plastic elements

 for i=1:Num_elem
     if element_tbl(i,4)<=2
        Num_plast = Num_plast+1;
     end
 end

   fprintf(fid1,'\n');
   fprintf(fid1,'%6.3f   %6.3f  %6.3f  %6.3f\n',T_str,T_end,dt_cal,dt_out);
   fprintf(fid1,'%6.3f   %6.3f\n',beta_val,gamma_val);
   fprintf(fid1,'\n');

   fprintf(fid1,'%4d   %4d   %4d\n',Num_story,Num_bay,Num_MRF);  
   fprintf(fid1,'%4d   %4d   %4d\n',Num_node,Num_DOF,Num_free);
   fprintf(fid1,'%4d   %4d   %4d\n',Num_act, Num_mst,Num_slv);
   fprintf(fid1,'%4d   %4d   %4d\n',Band_width,Num_elem,Num_plast);

   fprintf(fid1,'\n');

% --- Material Table ---
 for i=1:Num_elem
 if element_tbl(i,4)<=2
      fprintf(fid1,'%4d  %4d  %4d  %4d  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e', ...
                    element_tbl(i,1), element_tbl(i,2), ... 
                    element_tbl(i,3), element_tbl(i,4), ...
             elem_prop  (i,1), elem_prop  (i,6), elem_prop  (i,5), ...
             elem_prop  (i,2), elem_prop  (i,3), elem_prop  (i,4));
      fprintf(fid1,'  %12.5e  %12.5e  %12.5e  %12.5e %12.5e\n', ...
             elem_prop(i,7), elem_prop(i,8), elem_prop(i,9),c(i),s(i));
 end
 end
 fprintf(fid1,'\n');

% --- Master and Slave Data ---

 for i=1:Num_mst
     fprintf(fid1,'%4d  %4d  %4d', slv_tbl(i,1),slv_tbl(i,2),slv_tbl(i,3));
     for j=4:(3+Num_slv)
         fprintf(fid1,' %4d ', slv_tbl(i,j));
     end
     fprintf(fid1,'\n');
 end
 fprintf(fid1,'\n');

% --- Output Response Vector Index ---

      for i=0:fix(Num_DOF/10)-1
         fprintf(fid1,' %4d  %4d  %4d  %4d  %4d  %4d  %4d  %4d  %4d  %4d', ...
                 Ord_vec(10*i+1),Ord_vec(10*i+2),Ord_vec(10*i+3), ...
                 Ord_vec(10*i+4),Ord_vec(10*i+5), ...
                 Ord_vec(10*i+6),Ord_vec(10*i+7),Ord_vec(10*i+8), ...
                 Ord_vec(10*i+9),Ord_vec(10*i+10));
         fprintf(fid1,'\n');
      end
      if 10*(fix(Num_DOF/10)-1)+11 <= Num_DOF
         for i=10*(fix(Num_DOF/10)-1)+11:Num_DOF
            fprintf(fid1,' %4d ',Ord_vec(i));
         end 
      end
      fprintf(fid1,'\n\n');

% --- Reduced Matrix Index (After Removing Rigid Floor) ---

      for i=0:fix(Num_free/10)-1
         fprintf(fid1,' %4d  %4d  %4d  %4d  %4d  %4d  %4d  %4d  %4d  %4d', ...
                 diagT(10*i+1),diagT(10*i+2),diagT(10*i+3), ...
                 diagT(10*i+4),diagT(10*i+5), ...
                 diagT(10*i+6),diagT(10*i+7),diagT(10*i+8), ...
                 diagT(10*i+9),diagT(10*i+10));
         fprintf(fid1,'\n');
      end
      if 10*(fix(Num_free/10)-1)+11 <= Num_free
         for i=10*(fix(Num_free/10)-1)+11:Num_free
            fprintf(fid1,' %4d ',diagT(i));
         end 
      end
      fprintf(fid1,'\n\n');

%-----------------------------------
%   Sensor Positions for Aquitition
%-----------------------------------

    fprintf(fid10,'\n %d  %d  %d  %d\n\n',Num_snr,Num_obs,Num_cps,Num_cf);
    for i=1:Num_snr
       fprintf(fid10,'%d  %d  %d\n',snr(i,2),snr(i,3),snr(i,4));
    end
    fprintf(fid10,'\n');

%-------------------------------------
%   Observation Points for evaluation
%-------------------------------------

    for i=1:Num_obs
       fprintf(fid10,'%d  %d  %d\n',obs(i,2),obs(i,3),obs(i,4));
    end
    fprintf(fid10,'\n');

%----------------------------------------
%   Connection Points of Passive Devices
%----------------------------------------

    for i=1:Num_cps
       fprintf(fid10,'%d  %d  %d\n',cps(i,2),cps(i,3),cps(i,4));
    end
    fprintf(fid10,'\n');

%------------------------------
%   Location of Control Forces
%------------------------------

    for i=1:Num_cf
       fprintf(fid10,'%d  %d\n',cf(i,2),cf(i,3));
    end
    fprintf(fid10,'\n');

% -----------------------------------
    fclose (fid1);
    fclose (fid10);

% -----------------------------------
%  Set up SIMULINK Parameters
% -----------------------------------

    eval(Cnt_file);
