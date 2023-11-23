%==== This M.file declares the tracking problem structure ===
% here we declare a 2-D constant velocity model with coordinate
% measurements

%=== the following parameters is directly assessed by 
model.x_dim= 7;   %dimension of state vector
model.z_dim= 4;   %dimension of observation vector
% model.M_max= 20;      %max. number of targets

model.P_S= .99;  %probability of target death
model.P_D= .98;  %probability of detection in measurements

model.P_S_tempd = 0.96*model.P_S;
model.P_D_tempd = 0.96*model.P_D;

model.Q_S= 1-model.P_S; %probability of death
model.Q_D = 1-model.P_D;   %probability of misdetection

model.Q_S_tempd = 1-model.P_S_tempd;
model.Q_D_tempd = 1-model.P_D_tempd;


model.N_max = 50; %max number of terms to calculate the cardinality distribution
%WARNING: DO NOT set N_max  <= max |Z_k| (out of bounds error will occur)
%NOTE:    Advisable to set N_max >> max |Z_k| 
%         and              N_max >> max |X_k| 
%This will ensure capturing the tail of the clutter and state cardinality!!
model.lambda_c = 1; %make sure lambda_c << N_max !!!
model.log_lambda_c = log(model.lambda_c);
model.cdn_clutter = poisspdf([0:model.N_max]',model.lambda_c); %cardinality distribution of clutter - must be Poisson
model.log_cdn_clutter = log(model.cdn_clutter+eps(0));
model.clutterpdf = 1/Vc;
model.log_clutterpdf = log(model.clutterpdf);
model.run_flag = 'disp';
model.disp_flag = 'vis';

%=== the following parameters are more problem dependent
% 'siggen', 'phdfilter', etc. do not assess these parameters directly.
% These parameters are used by problem dependent functions such as
%
% gen_newstate, gen_birthstate, gen_observation, compute_likelihood
% 
% we create a structure array 'model' to store these specific parameters

%===here we set up the state spc eqn. x= Ax_old + Bv
T= 1;   %sampling period
model.A = [1 0 0 0 1 0 0; 0 1 0 0 0 1 0; 0 0 1 0 0 0 1; 0 0 0 1 0 0 0; 0 0 0 0 1 0 0; 0 0 0 0 0 1 0; 0 0 0 0 0 0 1];
model.Q = [1 0 0 0 0 0 0; 0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 1 0 0 0; 0 0 0 0 0.01 0 0; 0 0 0 0 0 0.01 0; 0 0 0 0 0 0 0.0001];
% 
% A0= [ 1 T; 0 1 ];                       
% model.A= [ A0 zeros(2,2); zeros(2,2) A0 ];
% B0= [ (T^2)/2; T ];
% model.B= [ B0 zeros(2,1); zeros(2,1) B0 ];
% model.sigma_v= 5;
% model.Q= (model.sigma_v)^2* model.B*model.B';

%=== parameters for the observation
model.C_posn = [1 0 0 0 0 0 0; 0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 1 0 0 0];
model.R = [1 0 0 0; 0 1 0 0; 0 0 10 0; 0 0 0 10];

% model.C_posn = [1 0 0 0; 0 1 0 0];
% model.R_posn = diag([10; 10]);
model.R_posn = diag([10; 10]);


L_s = 0; model.lambda_s = []; %for compatibility with phd filter code