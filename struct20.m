%-----------------------------------------------------------------------
%|*********************************************************************|
%|*           Structural Data of the 20 Story 5 Bay Building          *|
%|*                                                                   *|
%|*                  University of Notre Dame                         *|
%|*                       November, 1999                              *|
%|*                                                                   *|
%|*               Coded by      Y.Ohtori                              *|
%|*                             R.E.Christenson                       *|
%|*               Supervised by B.F.Spencer, Jr.  	              *|
%|*********************************************************************|
%-----------------------------------------------------------------------

% UNITS: N - m

% ------------------------------------------------
% --- Define File Name of Structural Data File ---
% ---  for the SIMULINK S-function             ---
% ------------------------------------------------
     f1_name = 'structure.dat1';       % Name of Structural Data for SIMULINK

% --------------------------------------------
% --- Define Initial Data for Calculations ---
% --------------------------------------------
     T_str   = 0.0;                     % Start time of the Response Analysis
     T_end   = 100.0;                   % End time of the Response Analysis
     dt_cal  = 0.005;                   % Calculation Interval Delta_t
     dt_out  = 0.01;                   % Output Interval for file

     beta_val  = 1.0/4.0;              % beta  value for Newmark-beta Method
     gamma_val = 1.0/2.0;              % gamma value for Newmark-beta Method

% ----------------------
% --- Define Damping ---
% ----------------------
     zeta_cr   = 0.02;                 % Critical Damping
     h_max     = 0.02;                 % Maximum Damping for Type 2
     nCutOff   = 5;                    % Number of mode for cutting off

% -------------------------------------
% --- Set Data for Nodal Coordinate ---
% -------------------------------------
     Num_story = 20+2;                       % Number of stories + basements
     Num_bay   = 5;                          % Number of Bays
     Num_MRF   = 2;                          % Number of Moment Resisting Frame
     Num_node  = (Num_bay+1)*(Num_story+1);  % Number of Nodes
     Num_DOF   = 3*Num_node;                 % Number of Degree of Freedoms

     height = [ 3.6576 3.6576 5.4864 3.9624 3.9624 3.9624 3.9624 ...
                3.9624 3.9624 3.9624 3.9624 3.9624 3.9624 3.9624 ...
                3.9624 3.9624 3.9624 3.9624 3.9624 3.9624 3.9624 3.9624 ];  

     width  = [ 6.0960 6.0960 6.0960 6.0960 6.0960 ];

% ---------------------------------------------
% --- Define Material Table (left to right) ---
% ---------------------------------------------
     Num_beam  = Num_bay*Num_story;              % Number of Beams
     Num_col   = (Num_bay+1)*Num_story;          % Number of Coloms
     Num_elem  = Num_beam + Num_col;             % Number of Elements

     element_tbl = zeros(Num_elem,4);

%  - Column -

   for i=1:Num_col
     element_tbl(i,1) = i;
     element_tbl(i,2) = i + (Num_bay+1);
     element_tbl(i,4) = 1;
   end

  element_tbl(1:Num_col,3) = [  
                         1     2     2     2     2     1  ...   %  1st basement
                         1     2     2     2     2     1  ...   %  2nd basement
                         1     2     2     2     2     1  ...   %  1st floor
                        17     2     2     2     2    17  ...   %  2nd floor
                         3     2     2     2     2     3  ...   %  3rd floor
                         3     2     2     2     2     3  ...   %  4th floor
                        18    19    19    19    19    18  ...   %  5th floor
                         4     5     5     5     5     4  ...   %  6th floor
                         4     5     5     5     5     4  ...   %  7th floor
                         4     5     5     5     5     4  ...   %  8th floor
                         4     5     5     5     5     4  ...   %  9th floor
                         4     5     5     5     5     4  ...   %  10th floor
                         4    20    20    20    20     4  ...   %  11th floor
                         4     6     6     6     6     4  ...   %  12th floor
                         4     6     6     6     6     4  ...   %  13th floor
                        21    22    22    22    22    21  ...   %  14th floor
                         7     8     8     8     8     7  ...   %  15th floor
                         7     8     8     8     8     7  ...   %  16th floor
                         7    23    23    23    23     7  ...   %  17th floor
                         7     9     9     9     9     7  ...   %  18th floor
                        24    25    25    25    25    24  ...   %  19th floor
                        10    11    11    11    11    10  ...   %  20th floor
                                                        ]';

%  - Beam   -

   elem_no = Num_col;
   for i = 1:Num_story
     for j = 1:Num_bay
         elem_no = elem_no + 1;
         element_tbl(elem_no,1) = i*(Num_bay+1) + j;
         element_tbl(elem_no,2) = element_tbl(elem_no,1) + 1;
         element_tbl(elem_no,4) = 2;
     end
   end

  element_tbl(Num_col+1:Num_col+Num_beam,3) = [ 
                        12    12    12    12    12  ...   %  1st basement
                        12    12    12    12    12  ...   %  1st floor
                        12    12    12    12    12  ...   %  2nd floor
                        12    12    12    12    12  ...   %  3rd floor
                        12    12    12    12    12  ...   %  4th floor
                        12    12    12    12    12  ...   %  5th floor
                        13    13    13    13    13  ...   %  6th floor
                        13    13    13    13    13  ...   %  7th floor
                        13    13    13    13    13  ...   %  8th floor
                        13    13    13    13    13  ...   %  9th floor
                        13    13    13    13    13  ...   %  10th floor
                        13    13    13    13    13  ...   %  11th floor
                        12    12    12    12    12  ...   %  12th floor
                        12    12    12    12    12  ...   %  13th floor
                        12    12    12    12    12  ...   %  14th floor
                         8     8     8     8     8  ...   %  15th floor
                         8     8     8     8     8  ...   %  16th floor
                         8     8     8     8     8  ...   %  17th floor
                        14    14    14    14    14  ...   %  18th floor
                        14    14    14    14    14  ...   %  19th floor
                        15    15    15    15    15  ...   %  20th floor
                        16    16    16    16    16  ...   %  roof
                         			       ]'; 

  element_tbl(Num_col+[1 2 3 4 5],4) = 3;

% -----------------------------------
% --- Define Master - Slave Nodes ---
% -----------------------------------
  slv_tbl = [
%   Master  Dir  Num_slv  Slv_1  Slv_2  Slv_3  Slv_4  Slv_5
      9      1     5        7      8      10     11     12   %  1st basement
     21      1     5       19     20      22     23     24   %  2nd floor
     27      1     5       25     26      28     29     30   %  3rd floor
     33      1     5       31     32      34     35     36   %  4th floor
     39      1     5       37     38      40     41     42   %  5th floor
     45      1     5       43     44      46     47     48   %  6th floor
     51      1     5       49     50      52     53     54   %  7th floor
     57      1     5       55     56      58     59     60   %  8th floor
     63      1     5       61     62      64     65     66   %  9th floor
     69      1     5       67     68      70     71     72   %  10th floor
     75      1     5       73     74      76     77     78   %  11th floor
     81      1     5       79     80      82     83     84   %  12th floor
     87      1     5       85     86      88     89     90   %  13th floor
     93      1     5       91     92      94     95     96   %  14th floor
     99      1     5       97     98     100    101    102   %  15th floor
    105      1     5      103    104     106    107    108   %  16th floor
    111      1     5      109    110     112    113    114   %  17th floor
    117      1     5      115    116     118    119    120   %  18th floor
    123      1     5      121    122     124    125    126   %  19th floor
    129      1     5      127    128     130    131    132   %  20th floor
    135      1     5      133    134     136    137    138   %  roof
						          ];      
  Num_mst = length(slv_tbl(:,1));                    % Number of Master Nodes

% --------------------------
% --- Define Fixed Nodes ---
% --------------------------
%                       Fix: 1 , Free: 0 
%              Node_no  Hor.  Vert.  Rot.
  Fix_node = [    1       1     1     0 ;
                  2       1     1     0 ;
                  3       1     1     0 ;
                  4       1     1     0 ;
                  5       1     1     0 ;
                  6       1     1     0 ;
                 13       1     0     0 ;
                 14       1     0     0 ;
                 15       1     0     0 ;
                 16       1     0     0 ;
                 17       1     0     0 ;
                 18       1     0     0 ;
                  ];
  Num_BND   = length(Fix_node(:,1));         % Number of Fixed Nodes

% ---------------------------
% --- Define Element mass ---
% ---------------------------
  seismic_mass1 = 2.6561e5;
  seismic_mass2 = 2.8166e5;
  seismic_mass3 = 2.7582e5;
  seismic_mass4 = 2.9188e5;

  Element_mass = [
   zeros(1,Num_col)...
   seismic_mass1*(1/(Num_bay))*ones(1,Num_bay)...   %  1st basement
   seismic_mass1*(1/(Num_bay))*ones(1,Num_bay)...   %  1st floor
   seismic_mass2*(1/(Num_bay))*ones(1,Num_bay)...   %  2nd floor
   seismic_mass3*(1/(Num_bay))*ones(1,Num_bay)...   %  3rd floor
   seismic_mass3*(1/(Num_bay))*ones(1,Num_bay)...   %  4th floor
   seismic_mass3*(1/(Num_bay))*ones(1,Num_bay)...   %  5th floor
   seismic_mass3*(1/(Num_bay))*ones(1,Num_bay)...   %  6th floor
   seismic_mass3*(1/(Num_bay))*ones(1,Num_bay)...   %  7th floor
   seismic_mass3*(1/(Num_bay))*ones(1,Num_bay)...   %  8th floor
   seismic_mass3*(1/(Num_bay))*ones(1,Num_bay)...   %  9th floor
   seismic_mass3*(1/(Num_bay))*ones(1,Num_bay)...   %  10th floor
   seismic_mass3*(1/(Num_bay))*ones(1,Num_bay)...   %  11th floor
   seismic_mass3*(1/(Num_bay))*ones(1,Num_bay)...   %  12th floor
   seismic_mass3*(1/(Num_bay))*ones(1,Num_bay)...   %  13th floor
   seismic_mass3*(1/(Num_bay))*ones(1,Num_bay)...   %  14th floor
   seismic_mass3*(1/(Num_bay))*ones(1,Num_bay)...   %  15th floor
   seismic_mass3*(1/(Num_bay))*ones(1,Num_bay)...   %  16th floor
   seismic_mass3*(1/(Num_bay))*ones(1,Num_bay)...   %  17th floor
   seismic_mass3*(1/(Num_bay))*ones(1,Num_bay)...   %  18th floor
   seismic_mass3*(1/(Num_bay))*ones(1,Num_bay)...   %  19th floor
   seismic_mass3*(1/(Num_bay))*ones(1,Num_bay)...   %  20th floor
   seismic_mass4*(1/(Num_bay))*ones(1,Num_bay)...   %  roof
                 				  ]';

% -----------------------------
% --- Define Material Table --- 
% -----------------------------
% Type = 0: Spread Plasticity Model
%        1: Concentrated Plasticity Model
  mat_tbl = [
%Mat_no  EI1       EI2        EI3        EA         GA         d1        d2   type mod.
 1  2.4959e+08 2.2584e+08 7.4866e+06 1.3416e+10 8.8960e+15 1.2659e-02 1.3325e-02  1
 2  9.9037e+08 8.9628e+08 2.9710e+07 1.2693e+10 8.8960e+15 6.3675e-03 6.7026e-03  1
 3  1.8176e+08 1.6450e+08 5.4532e+06 8.8686e+09 8.8960e+15 1.2090e-02 1.2726e-02  1
 4  1.5305e+08 1.3847e+08 4.5904e+06 7.2239e+09 8.8960e+15 1.1899e-02 1.2525e-02  1
 5  6.3667e+08 5.7617e+08 1.9099e+07 8.6687e+09 8.8960e+15 6.5645e-03 6.9100e-03  1
 6  5.2099e+08 4.7145e+08 1.5629e+07 7.2626e+09 8.8960e+15 6.6336e-03 6.9828e-03  1
 7  1.2076e+08 1.0927e+08 3.6222e+06 5.5147e+09 8.8960e+15 1.1708e-02 1.2324e-02  1
 8  3.3456e+08 3.0277e+08 1.0036e+07 4.9664e+09 8.8960e+15 6.8374e-03 7.1972e-03  1
 9  2.9462e+08 2.6660e+08 8.8380e+06 4.4375e+09 8.8960e+15 6.8621e-03 7.2233e-03  1
 10 8.4639e+07 7.6629e+07 2.5401e+06 3.7410e+09 8.8960e+15 1.1518e-02 1.2124e-02  1
 11 1.9724e+08 1.7849e+08 5.9170e+06 3.1863e+09 8.8960e+15 7.0212e-03 7.3908e-03  1
 12 3.3207e+08 2.8596e+08 9.9615e+06 3.7539e+09 8.8960e+15 4.9618e-03 5.2362e-03  1
 13 3.7201e+08 3.5399e+08 1.1160e+07 4.0892e+09 8.8960e+15 4.9116e-03 5.1575e-03  1
 14 2.3719e+08 2.2226e+08 7.1154e+06 3.1992e+09 8.8960e+15 5.4325e-03 5.7087e-03  1
 15 1.2900e+08 1.2362e+08 3.8698e+06 2.3478e+09 8.8960e+15 6.2635e-03 6.5748e-03  1
 16 8.1893e+07 7.1509e+07 2.4567e+06 1.8963e+09 8.8960e+15 7.0934e-03 7.4803e-03  1
                ];
% Add splice elements
  mat_tbl = [mat_tbl
17         (6/13)*mat_tbl(1,2:8)+(7/13)*mat_tbl(3,2:8)                            1
18         (6/13)*mat_tbl(3,2:8)+(7/13)*mat_tbl(4,2:8)                            1
19         (6/13)*mat_tbl(2,2:8)+(7/13)*mat_tbl(5,2:8)                            1
20         (6/13)*mat_tbl(5,2:8)+(7/13)*mat_tbl(6,2:8)                            1
21         (6/13)*mat_tbl(4,2:8)+(7/13)*mat_tbl(7,2:8)                            1
22         (6/13)*mat_tbl(6,2:8)+(7/13)*mat_tbl(8,2:8)                            1
23         (6/13)*mat_tbl(8,2:8)+(7/13)*mat_tbl(9,2:8)                            1
24         (6/13)*mat_tbl(7,2:8)+(7/13)*mat_tbl(10,2:8)                           1
25         (6/13)*mat_tbl(9,2:8)+(7/13)*mat_tbl(11,2:8)                           1
                ];
