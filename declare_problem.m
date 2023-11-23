%==== This M.file declares the tracking problem structure ===
% here we declare a 2-D constant velocity model with coordinate
% measurements

%=== the following parameters is directly assessed by 
% siggen.m and phdfilter.m
% K= 100;  %data length
x_dim= 7;   %dimension of state vector
z_dim= 4;   %dimension of observation vector
M_max= 20;      %max. number of targets

P_S= .99;  %probability of target death
% P_D= .80;  %probability of detection in measurements
P_D = .98;

P_S_tempd = 0.96*P_S;
P_D_tempd = 0.96*P_D;

Q_S= 1-P_S; %probability of death
Q_D = 1-P_D;   %probability of misdetection

Q_S_tempd = 1-P_S_tempd;
Q_D_tempd = 1-P_D_tempd;


N_max = 100; %max number of terms to calculate the cardinality distribution
%WARNING: DO NOT set N_max  <= max |Z_k| (out of bounds error will occur)
%NOTE:    Advisable to set N_max >> max |Z_k| 
%         and              N_max >> max |X_k| 
%This will ensure capturing the tail of the clutter and state cardinality!!
lambda_c = 1; %make sure lambda_c << N_max !!!
log_lambda_c = log(lambda_c);
cdn_clutter = poisspdf([0:N_max]',lambda_c); %cardinality distribution of clutter - must be Poisson
log_cdn_clutter = log(cdn_clutter+eps(0));
range_c= [ 0 576; 0 768 ];    %clutter intervals
clutterpdf = 1/prod(range_c(:,2)-range_c(:,1));
log_clutterpdf = log(clutterpdf);
run_flag = 'disp';

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
model.R_posn = diag([10; 10]);
% model.C_posn= [ 1 0 0 0 ; 0 0 1 0 ];  
% model.D= diag([ 10; 10 ]);   %std for angle and range noise
% model.R= model.D*model.D';

%===MB BIRTH --- here is the parameter for the birth target states
% L_b= 4;         %no. of MB birth terms
% 
% N_birth= zeros(L_b,1);
% model.bar_q= zeros(L_b,1);
% model.lambda_b= cell(L_b,1);
% model.bar_x= cell(L_b,1);
% model.bar_B= cell(L_b,1); 
% model.bar_Q= cell(L_b,1);
% 
% N_birth(1)= 1;          %no. of Gaussians in birth term 1
% model.bar_q(1)=0.03;                            %prob of existence for term 1
% model.lambda_b{1}(1)= 1;                        %weight of Gaussians
% model.bar_x{1}(:,1)= [ 0; 0; 0; 0 ];            %mean of Gaussians
% model.bar_B{1}(:,:,1)= diag([ 10; 10; 10; 10 ]);%std of Gaussians
% 
% N_birth(2)= 1;          %no. of Gaussians in birth term 2
% model.bar_q(2)=0.03;                            %prob of existence for term 2
% model.lambda_b{2}(1)= 1;                        %weight of Gaussians
% model.bar_x{2}(:,1)= [ 400; 0; -600; 0 ];       %mean of Gaussians
% model.bar_B{2}(:,:,1)= diag([ 10; 10; 10; 10 ]);%std of Gaussians
% 
% N_birth(3)= 1;          %no. of Gaussians in birth term 3
% model.bar_q(3)=0.03;                            %prob of existence for term 3
% model.lambda_b{3}(1)= 1;                        %weight of Gaussians
% model.bar_x{3}(:,1)= [ -800; 0; -200; 0 ];      %mean of Gaussians
% model.bar_B{3}(:,:,1)= diag([ 10; 10; 10; 10 ]);%std of Gaussians
% 
% N_birth(4)= 1;          %no. of Gaussians in birth term 4
% model.bar_q(4)=0.03;                            %prob of existence for term 4
% model.lambda_b{4}(1)= 1;                        %weight of Gaussians
% model.bar_x{4}(:,1)= [ -200; 0; 800; 0 ];       %mean of Gaussians
% model.bar_B{4}(:,:,1)= diag([ 10; 10; 10; 10 ]);%std of Gaussians
% 
% for i=1:L_b
%     for g=1:N_birth(i)
%     model.bar_Q{i}(:,:,g)= model.bar_B{i}(:,:,g)*model.bar_B{i}(:,:,g)'; %cov of Gaussians
%     end;
% end;


L_s = 0; model.lambda_s = []; %for compatibility with phd filter code