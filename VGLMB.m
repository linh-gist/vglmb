

% This is a routine for Visual GLMB tracker
% declare_problem.m is changed to setFilterParameters.m

%==== M-Best- code base from CPHD_LIN_F3

%==== Run parameters
function track = VGLMB(observation,input_frames,model,other_param)
% 
% marg_flag= 0; %0/1 on/off for MDGLMB approximation
% elim_threshold= 1e-3; %for pruning of Gaussians inside tracks
% merge_threshold= 5; %for merging of Gaussians inside tracks
% cap_threshold= 5; %for capping of Gaussians inside tracks

%KCF parameters
kernel_type = 'gaussian';
kernel.type = kernel_type;

feature_type = 'hog';
show_visualization = 1;
show_plots = 1;
interp_factor = 0.02;
		
kernel.sigma = 0.5;
kernel.poly_a = 1;
kernel.poly_b = 9;

features.gray = false;
features.hog = true;
features.hog_orientations = 9;

padding = 1.5;  %extra area surrounding the target
lambda = 1e-4;  %regularization
output_sigma_factor = 0.1;  %spatial bandwidth (proportional to target)
cell_size = 4;


% Pre-process detection data
observation = cutDetections(observation,other_param);

K = max(observation.fr);
sc_threshold = quantile(observation.r,0.15);

detection = cell(K,1);
for t=1:K
    
    idx = find(observation.fr == t);
    cur_observation.x = observation.x(idx);
    cur_observation.y = observation.y(idx);
    cur_observation.bx = observation.bx(idx);
    cur_observation.by = observation.by(idx);
    cur_observation.w = observation.w(idx);
    cur_observation.h = observation.h(idx);
    cur_observation.sc = observation.r(idx);
    cur_observation.fr = observation.fr(idx);
        
    detection{t} = [cur_observation.x cur_observation.y cur_observation.w cur_observation.h cur_observation.sc]';
    
end



Hbes= 500; %cap number of components
% Hbesreq= [30 100 200];  %min num generated internally based on cardinality std dev
Hbesreq= [100 200 500];  %min num generated internally based on cardinality std dev
stdvpts= [2 2 2]; %cardinality std dev cut off points for internal requested number of components
Cmin= 0; %min components per cardinality for pruning
chop_threshold = 1e-5; %pruning threshold for component weights

gate_flag= 1;               %0/1 on/off
P_G= 0.9999999;             %gate probability (for analog and digital gating)  
gamma= chi2inv(P_G,model.z_dim);  %for analog gating only - inv chi^2 dn gamma value
Z_orig = detection;                 %copy origina l measurements - Z{k} will be replaced with gated measurements
origclutrate= model.lambda_c; log_origclutrate= model.log_lambda_c; %original clutter rate
origclutpdf= model.clutterpdf; log_origclutpdf= model.log_clutterpdf; %original clutter pdf
Z_birt = cell(size(detection));
Z_q = cell(size(detection));

%==== Start filtering
N_max = model.N_max;
P_S = model.P_S;
P_D = model.P_D;
Q_S = model.Q_S;
Q_D = model.Q_D;
P_S_tempd = model.P_S_tempd;
P_D_tempd = model.P_D_tempd;
Q_S_tempd = model.Q_S_tempd;
Q_D_tempd = model.Q_D_tempd;

hat_N= zeros(K,1);
hat_X= cell(K,1);
hat_T= cell(K,1);

track.x = [];
track.y = [];
track.x_hat = [];
track.y_hat = [];
track.w = [];
track.h = [];
track.fr = [];
track.id = [];
% track.sc = [];
% track.isdummy = [];

cdn_update_stack  = zeros(N_max+1,K);
cdn_update_mean   = zeros(K,1);
cdn_update_mode   = zeros(K,1);
cdn_update_var    = zeros(K,1);

emm=0;
Z = cell(K,1);
for k=1:K
    emm= max(emm,size(detection{k},2));
    Z{k} = [detection{k}(1:2,:); detection{k}(1:2,:)+detection{k}(3:4,:); detection{k}(5,:)];
end
tmp_N_max= max(emm,N_max)+3;

%---precalculate constants used prediction and update calculations
nvector  = [0:N_max]';
logPSpow = [0:tmp_N_max]'*log(P_S);
logQSpow = [0:tmp_N_max]'*log(1-P_S);
logPDpow = [0:tmp_N_max]'*log(P_D);
logQDpow = [0:tmp_N_max]'*log(1-P_D);

% PULL THESE LINES OUTSIDE OF THIS SCRIPT FOR MC RUNS
% THESE ONLY DEPEND ON N_MAX AND DO NOT CHANGE FOR SIMULATION RUNS
% precalculate values of P(n,j) and C(n,j)
logCcoef = zeros(tmp_N_max+1,tmp_N_max+1);
logPcoef = zeros(tmp_N_max+1,tmp_N_max+1);
logfactorial = zeros(tmp_N_max+1,1);
logfactorial(1) = 0;
for n=1:tmp_N_max
    logfactorial(n+1) = log(n)+logfactorial(n);
end

for ell=0:tmp_N_max
    for j=0:ell
        logPcoef(ell+1,j+1) = logfactorial(ell+1)-logfactorial(ell-j+1);
        logCcoef(ell+1,j+1) = logPcoef(ell+1,j+1)-logfactorial(j+1);
    end
end
% end calculations for P(n,j) and C(n,j)

if ~exist('run_flag','var')
    run_flag = 'disp';
end

%--- i'll need to record these stuff for performance analysis
gaus_size_lg= zeros(K,1);
cpu_time_lg= zeros(K,1);

log_cdn_update = -realmax*ones(N_max+1,1); log_cdn_update(1)= 0;
cdn_update= zeros(N_max+1,1); cdn_update(1)= 1;
log_wtv_update= cell(N_max+1,1);  log_wtv_update{1}= 0;
wtv_update= cell(N_max+1,1);  wtv_update{1}= 1;
tracks_update= cell(0,1); 
hyps_update= cell(N_max+1,1);


track_count = 0;
for k=1:50
    time_start= cputime;
    
    offset = 0;
    try
        if strcmp(other_param.seq,'PETS2009') || strcmp(other_param.seq,'KITTI_train') || strcmp(other_param.seq,'KITTI_test')
            img = imread(sprintf(input_frames, k-1+offset)); %% read an image
        else
            img = imread(sprintf(input_frames, k+offset)); %% read an image
        end
    catch
        break
    end
    
    frame = img;
    iou_threshold_birth = 0.1;
    Z_gt = convert_bbox_to_z(Z{k});
%     m = size(Z_gt,2); 
    
    if gate_flag          
%         if m~=0
                
        %find measurements for measurement driven birth density
        %(operates track by track, not by predicted PHD)
            if k==1
                Z_birt{k} = Z_gt;
                Z_q{k} = .001*ones(size(Z_birt{k},2),1);
            elseif isempty(hat_X{k-1})
                Z_birt{k} = Z_gt;
            else
                valid_idx1 = [];
                m = size(Z{k-1},2);
                if m~=0
                    iou_mat = zeros(size(Z{k-1},2),size(hat_X{k-1},2));
                    for d=1:size(Z{k-1},2)
                        for t=1:size(hat_X{k-1},2)
                            iou_mat(d,t) = iou(Z{k-1}(:,d),convert_x_to_bbox(hat_X{k-1}(:,t))); % why? because we consider birth from Z{k-1}
                        end
                    end
                    [valid_idx,~] = find(iou_mat > iou_threshold_birth);
                    valid_idx1 = setdiff([1:m],valid_idx);
                    Z_birt{k} = convert_bbox_to_z(Z{k-1}(:,valid_idx1));
                    if isempty(valid_idx1)
                        
                    end
                    Z_q{k} = iou_mat(valid_idx1,:)';
    %                 Z_tr{k} = Z_gt(:,valid_idx);
                end
            end
%         end
    end
    
    
    % Construct measurement driven birth density
    Z_birth = Z_birt{k};
    L_b= size(Z_birth,2); % number of MB birth
    N_birth= zeros(L_b,1);
    model.bar_q= zeros(L_b,1);
    model.bar_q_temped = zeros(L_b,1);
    model.lambda_b= cell(L_b,1);
    model.bar_x= cell(L_b,1);
    model.bar_B= cell(L_b,1); 
    model.bar_Q= cell(L_b,1);
    
    for i=1:L_b
        N_birth(i) = 1;
        model.bar_q(i)= 0.01;
%         model.bar_q(i)= max(0.01,0.1-Z_q{k}(i));
        model.bar_q_tempd(i) = 0.05;
        model.lambda_b{i}(1)= 1;
        model.bar_x{i}(:,1) = [Z_birth(1,i);Z_birth(2,i);Z_birth(3,i);Z_birth(4,i);0;0;0];
        model.bar_B{i}(:,:,1)= diag([10; 10; 10; 10; 10000; 10000; 10000]);
    end   
    for i=1:L_b
        for g=1:N_birth(i)
            model.bar_Q{i}(:,:,g)= model.bar_B{i}; %cov of Gaussians
        end
    end
    model.N_birth= N_birth;
      
    %---stochastic component selection/allocation
    prev_cvar= ([0:N_max].^2*cdn_update(:)) - ([0:N_max]*cdn_update(:))^2;
    prev_stdv= sqrt(prev_cvar);
    Hbesuse= Hbesreq(find(cumsum(stdvpts)>prev_stdv,1)); 
    hbescell= cell(N_max+1,1);
    sampledn= resample(cdn_update,Hbesuse);
    for n=0:N_max
        nidx= n+1;
        hbescell{nidx}= zeros(length(wtv_update{nidx}),1);
        sampledc= resample(wtv_update{nidx},nnz(sampledn==nidx));
        for cidx=1:length(wtv_update{nidx})
            hbescell{nidx}(cidx)= nnz(sampledc==cidx);
        end
    end
    
    %---update
    %init params
    log_cdn_temp = -realmax*ones(N_max+1,1); log_wtv_temp= cell(N_max+1,1); tracks_predict= cell(length(model.bar_q)+length(tracks_update),1); tracks_temp= cell(length(tracks_predict)*(1+size(Z{k},2)),1); hyps_temp= cell(N_max+1,1);

    %h-best update
    
        %create birth tracks
        for tabbidx=1:length(model.bar_q)
            track_count= track_count + 1;
            tracks_predict{tabbidx}.m = model.bar_x{tabbidx};
            tracks_predict{tabbidx}.P = model.bar_Q{tabbidx};
            tracks_predict{tabbidx}.w = model.lambda_b{tabbidx}(:);
            tracks_predict{tabbidx}.l = [k;tabbidx;track_count];
            tracks_predict{tabbidx}.ah = [];
            tracks_predict{tabbidx}.avps = model.bar_q(tabbidx);
            tracks_predict{tabbidx}.avqs = 1-model.bar_q(tabbidx);
            tracks_predict{tabbidx}.avps_tempd = model.bar_q_tempd(tabbidx);
            tracks_predict{tabbidx}.avqs_tempd = 1-model.bar_q_tempd(tabbidx);
            
            pos =  tracks_predict{tabbidx}.m(1:2)';
            w = sqrt(tracks_predict{tabbidx}.m(3)*tracks_predict{tabbidx}.m(4)); h = tracks_predict{tabbidx}.m(3)./w; 
            target_sz = [h w];
            ref = kcf_tracker_init(img, pos, target_sz, ...
                        padding, kernel, lambda, output_sigma_factor, interp_factor, ...
                        cell_size, features);
            tracks_predict{tabbidx}.ref = ref;                     
        end
        
        %create predicted surviving tracks
        for tabsidx=1:length(tracks_update)
            offset= length(model.bar_q);
            [wtemp_predict,mtemp_predict,Ptemp_predict]= kalman_predict_sum_AS(1,model.A,model.Q,tracks_update{tabsidx}.w,tracks_update{tabsidx}.m,tracks_update{tabsidx}.P,tracks_update{tabsidx}.ah);
            tracks_predict{tabsidx+offset}.m = mtemp_predict;
            tracks_predict{tabsidx+offset}.P = Ptemp_predict;
            tracks_predict{tabsidx+offset}.w = wtemp_predict;
            tracks_predict{tabsidx+offset}.l = tracks_update{tabsidx}.l;
            tracks_predict{tabsidx+offset}.ah = tracks_update{tabsidx}.ah;
            tracks_predict{tabsidx+offset}.avps = P_S;
            tracks_predict{tabsidx+offset}.avqs = Q_S;
            tracks_predict{tabsidx+offset}.avps_tempd = P_S_tempd;
            tracks_predict{tabsidx+offset}.avqs_tempd = Q_S_tempd;
            tracks_predict{tabsidx+offset}.ref = tracks_update{tabsidx}.ref;
        end
        
        %track level gating
        %no of measurements
        offset_tr = 0;
        for tabidx=1:length(tracks_predict)
           if (tracks_predict{tabidx}.ah == 0)
%            if isempty(tracks_predict{tabidx}.ah)
                offset_tr = offset_tr +1; 
           end
        end
        m = size(Z_gt,2);
        iou_threshold = 0.1;
        
        %gating based on IOU
        valid_idx = []; 
        if gate_flag          
            if m~=0
                %find measurements for measurement driven birth density
                %(operates track by track, not by predicted PHD)
                iou_mat = zeros(size(Z_gt,2),track_count);
                if k==1
                    for d=1:size(Z_gt,2)
                        for t=1:length(model.bar_q)
                            iou_mat(d,t) = iou(convert_x_to_bbox(Z_gt(:,d)),convert_x_to_bbox(model.bar_x{t}));
                        end
                    end
                else
                    
                    for d=1:size(Z_gt,2)
                        for t=1:length(tracks_predict)-offset_tr
                            iou_mat(d,t) = iou(convert_x_to_bbox(Z_gt(:,d)),convert_x_to_bbox(tracks_predict{t+offset_tr}.m));
                        end
                    end
                end
%                 for t=1:length(tracks_predict)
%                     tracks_predict{t}.gate_meas= [];
%                     tracks_predict{t}.gate_meas= union(tracks_predict{t}.gate_meas,find(iou_mat(:,t) > iou_threshold ));                    
%                 end
                
                [valid_idx,~] = find(iou_mat > iou_threshold);
                Zgated{k} = Z_gt(:,valid_idx);
                Zz = Zgated{k};
            else
                Zz = [];
            end
        end
        
        cur_observation = NMS(detection{k}(:,unique(valid_idx)));
        Zt = convert_bbox_to_z(cur_observation);
        
        m = size(Zt,2);      
        tracks_temp= cell(length(tracks_predict)*(1+size(Zt,2)),1);
        %create temporary updated update tracks (legacy ones first)
        for tabidx= 1:length(tracks_predict)
            tracks_temp{tabidx}= tracks_predict{tabidx};
            tracks_temp{tabidx}.qz= NaN;
            tracks_temp{tabidx}.ah= [tracks_predict{tabidx}.ah(:); 0];
            tracks_temp{tabidx}.avpd= P_D;
            tracks_temp{tabidx}.avqd= Q_D;
            tracks_temp{tabidx}.avpd_tempd= P_D_tempd;
            tracks_temp{tabidx}.avqd_tempd= Q_D_tempd;
            tracks_temp{tabidx}.used= 0;
        end
        
        %create temporary updated update tracks (now measurement updated ones, organized in blocks of predicted tracks, one for each received measurement)
        for emm= 1:m
            for tabidx= 1:length(tracks_predict)
%                 if ~gate_flag || any(emm == tracks_predict{tabidx}.gate_meas) %if gating is off do all updates automatically, or if meas is validated do single target update
                if gate_flag && Zt(5,emm) > sc_threshold %if gating is off do all updates automatically, or if meas is validated do single target update
                    stoidx= length(tracks_predict)*emm + tabidx; %index of predicted track i updated with measurement j is (number_predicted+tracks*j + i)
                    [wtemp_update,mtemp_update,Ptemp_update] = kalman_update_sum(Zt(1:4,emm),1,model.C_posn,zeros(model.z_dim,1),model.R,tracks_predict{tabidx}.w,tracks_predict{tabidx}.m,tracks_predict{tabidx}.P);
                    tracks_temp{stoidx}= tracks_predict{tabidx};
                    tracks_temp{stoidx}.m = mtemp_update;
                    tracks_temp{stoidx}.P = Ptemp_update;
                    
                    iou_tr= iou(convert_x_to_bbox(Zt(:,emm)),convert_x_to_bbox(tracks_predict{tabidx}.m));
                    tracks_temp{stoidx}.iou = P_D*iou_tr;
                    tracks_temp{stoidx}.iou_tempd = P_D_tempd*iou_tr;

                    tracks_temp{stoidx}.qz = P_D*sum(wtemp_update); % bayes evidence for construction of cost matrix in data association step
                    tracks_temp{stoidx}.qz_tempd = P_D_tempd*sum(wtemp_update); % tempered bayes evidence for construction of cost matrix in data association step
                    
                    
                    tracks_temp{stoidx}.w = wtemp_update/sum(wtemp_update);
                    tracks_temp{stoidx}.ah= [tracks_predict{tabidx}.ah(:); emm];
                    tracks_temp{stoidx}.avpd_tempd= P_D_tempd;
                    tracks_temp{stoidx}.avqd_tempd= Q_D_tempd;
                    tracks_temp{stoidx}.used= 0;
                    tracks_temp{stoidx}.ref = tracks_predict{tabidx}.ref;
                    
                elseif Zt(5,emm) < sc_threshold
                    stoidx= length(tracks_predict)*emm + tabidx; %index of predicted track i updated with measurement j is (number_predicted+tracks*j + i)
                    
                    % Start KCF filtering
                    ref_model = tracks_predict{tabidx}.ref;
                    target_sz = ref_model.size;                    
                    pos =  tracks_predict{tabidx}.m(1:2)'; 

                    positions = kcf_tracker(img, ref_model, pos, target_sz, ...
                        padding, kernel, lambda, output_sigma_factor, interp_factor, ...
                        cell_size, features);        
                    Z_kcf = [positions, prod(target_sz), target_sz(2)/target_sz(1), Zt(5,emm)]';
                    
                    [wtemp_update,mtemp_update,Ptemp_update] = kalman_update_sum(Z_kcf(1:4),1,model.C_posn,zeros(model.z_dim,1),model.R,tracks_predict{tabidx}.w,tracks_predict{tabidx}.m,tracks_predict{tabidx}.P);
                    tracks_temp{stoidx}= tracks_predict{tabidx};
                    tracks_temp{stoidx}.m = mtemp_update;
                    tracks_temp{stoidx}.P = Ptemp_update;
                    
                    iou_tr= iou(convert_x_to_bbox(Z_kcf),convert_x_to_bbox(tracks_predict{tabidx}.m));
                    tracks_temp{stoidx}.iou = P_D*iou_tr;
                    tracks_temp{stoidx}.iou_tempd = P_D_tempd*iou_tr;

                    tracks_temp{stoidx}.qz = P_D*sum(wtemp_update); % bayes evidence for construction of cost matrix in data association step
                    tracks_temp{stoidx}.qz_tempd = P_D_tempd*sum(wtemp_update); % tempered bayes evidence for construction of cost matrix in data association step
                                      
                    tracks_temp{stoidx}.w = wtemp_update/sum(wtemp_update);
                    tracks_temp{stoidx}.ah= [tracks_predict{tabidx}.ah(:); emm];
                    tracks_temp{stoidx}.avpd_tempd= P_D_tempd;
                    tracks_temp{stoidx}.avqd_tempd= Q_D_tempd;
                    tracks_temp{stoidx}.used= 0;
                    tracks_temp{stoidx}.ref = tracks_predict{tabidx}.ref;
                    
                else
                    %meas is not validated, don't bother updating, set evidence to zero
                    stoidx= length(tracks_predict)*emm + tabidx; %index of predicted track i updated with measurement j is (number_predicted+tracks*j + i)
                    tracks_temp{stoidx}.qz= 0; %bayes evidence is identically zero
                    tracks_temp{stoidx}.iou= 0;
                    tracks_temp{stoidx}.used= 0;
                    tracks_temp{stoidx}.gate_meas= [];
                end
            end
        end

        %update hypotheses
        for n= 0:N_max %loop over updated cardinality
                nidx= n+1;
                numcmp= length(hyps_update{nidx});
                if n==0 
                   numcmp=1; %trick to force loop entry for 0 cardinality prediction and  update using same code 
                end
                    
                for cidx= 1:numcmp %loop over all components
                    hbes= hbescell{nidx}(cidx); %number of h-best to generate (use to allocate proportionally)
                    if hbes ~= 0
                        %cost matrix
                        nbirthtracks= length(model.bar_q);
                        nexisttracks= n;
                        ntotaltracks= nbirthtracks+nexisttracks; 
                        PSvec= zeros(ntotaltracks,1);
                        PDvec= zeros(ntotaltracks,1);
                        QSvec= zeros(ntotaltracks,1);
                        QDvec= zeros(ntotaltracks,1);
                        costm= zeros(ntotaltracks,m);
                        %tempered cost matrix
                        PSvec_tempd= zeros(ntotaltracks,1);
                        PDvec_tempd= zeros(ntotaltracks,1);
                        QSvec_tempd= zeros(ntotaltracks,1);
                        QDvec_tempd= zeros(ntotaltracks,1);
                        costm_tempd= zeros(ntotaltracks,m);
                        %calculate values for birth tracks
                        for bidx= 1:nbirthtracks
                            PSvec(bidx)= tracks_temp{bidx}.avps;
                            PDvec(bidx)= tracks_temp{bidx}.avpd;
                            QSvec(bidx)= tracks_temp{bidx}.avqs;
                            QDvec(bidx)= tracks_temp{bidx}.avqd;
                            PSvec_tempd(bidx)= tracks_temp{bidx}.avps_tempd;
                            PDvec_tempd(bidx)= tracks_temp{bidx}.avpd_tempd;
                            QSvec_tempd(bidx)= tracks_temp{bidx}.avqs_tempd;
                            QDvec_tempd(bidx)= tracks_temp{bidx}.avqd_tempd;
                            for emm= 1:m
                                linidx= length(tracks_predict)*emm+bidx;
                                if tracks_temp{linidx}.iou   %i.e. must be non-zero
                                    costm(bidx,emm)= PSvec(bidx)/QSvec(bidx)*tracks_temp{linidx}.iou/(model.lambda_c*model.clutterpdf*QDvec(bidx));
                                    costm_tempd(bidx,emm)= PSvec_tempd(bidx)/QSvec_tempd(bidx)*tracks_temp{linidx}.iou_tempd/(model.lambda_c*model.clutterpdf*QDvec_tempd(bidx));

%                                     costm(bidx,emm)= PSvec(bidx)/QSvec(bidx)*tracks_temp{linidx}.qz/(lambda_c*clutterpdf*QDvec(bidx));
%                                     costm_tempd(bidx,emm)= PSvec_tempd(bidx)/QSvec_tempd(bidx)*tracks_temp{linidx}.qz_tempd/(lambda_c*clutterpdf*QDvec_tempd(bidx));
                                end
                            end 
                        end
                        %calculate values for existing tracks
                        for tee= 1:n
                            PSvec(nbirthtracks+tee)= tracks_temp{nbirthtracks+hyps_update{nidx}{cidx}{tee}}.avps;
                            PDvec(nbirthtracks+tee)= tracks_temp{nbirthtracks+hyps_update{nidx}{cidx}{tee}}.avpd;
                            QSvec(nbirthtracks+tee)= tracks_temp{nbirthtracks+hyps_update{nidx}{cidx}{tee}}.avqs;
                            QDvec(nbirthtracks+tee)= tracks_temp{nbirthtracks+hyps_update{nidx}{cidx}{tee}}.avqd;
                            PSvec_tempd(nbirthtracks+tee)= tracks_temp{nbirthtracks+hyps_update{nidx}{cidx}{tee}}.avps_tempd;
                            PDvec_tempd(nbirthtracks+tee)= tracks_temp{nbirthtracks+hyps_update{nidx}{cidx}{tee}}.avpd_tempd;
                            QSvec_tempd(nbirthtracks+tee)= tracks_temp{nbirthtracks+hyps_update{nidx}{cidx}{tee}}.avqs_tempd;
                            QDvec_tempd(nbirthtracks+tee)= tracks_temp{nbirthtracks+hyps_update{nidx}{cidx}{tee}}.avqd_tempd;
                            for emm= 1:m
                                linidx= length(tracks_predict)*emm+nbirthtracks+hyps_update{nidx}{cidx}{tee};
                                if tracks_temp{linidx}.iou    %i.e. must be non-zero  
                                    costm(nbirthtracks+tee,emm)= PSvec(nbirthtracks+tee)/QSvec(nbirthtracks+tee)*tracks_temp{linidx}.iou/(model.lambda_c*model.clutterpdf*QDvec(nbirthtracks+tee));
                                    costm_tempd(nbirthtracks+tee,emm)= PSvec_tempd(nbirthtracks+tee)/QSvec_tempd(nbirthtracks+tee)*tracks_temp{linidx}.iou_tempd/(model.lambda_c*model.clutterpdf*QDvec_tempd(nbirthtracks+tee));                                    
%                                     costm(nbirthtracks+tee,emm)= PSvec(nbirthtracks+tee)/QSvec(nbirthtracks+tee)*tracks_temp{linidx}.qz/(lambda_c*clutterpdf*QDvec(nbirthtracks+tee));
%                                     costm_tempd(nbirthtracks+tee,emm)= PSvec_tempd(nbirthtracks+tee)/QSvec_tempd(nbirthtracks+tee)*tracks_temp{linidx}.qz_tempd/(lambda_c*clutterpdf*QDvec_tempd(nbirthtracks+tee));
                                end
                            end
                        end
                        costm= [diag(1./QDvec) diag(PSvec./QSvec) costm]; %[notsurvived survived_but_not_detected survived_detected_and_generated_measurement]
%                         costm_tempd= [diag(1./QDvec_tempd) diag(PSvec_tempd./QSvec_tempd) costm_tempd]; %[notsurvived survived_but_not_detected survived_detected_and_generated_measurement]
                        
                        neglogcostm= -log(costm); %DON'T transpose to leave tracks on rows instead of measurements (track to measurement assignment)
%                         neglogcostm_tempd= -log(costm_tempd); %DON'T transpose to leave tracks on rows instead of measurements (track to measurement assignment)
                         
                        %USE EITHER: gibbs sampling trick (one target per measurement, and one measurement per measurement)
%                         [assnmt,nlcost]= mbestwrap_updt_gibbsamp_tempered(neglogcostm,neglogcostm_tempd,hbes); rankscosts= exp(-nlcost);
                        [assnmt,nlcost]= mbestwrap_updt_gibbsamp(neglogcostm,hbes); rankscosts= exp(-nlcost);

                        
%                         %OR: optimal assignment trick (one target per measurement, and one measurement per measurement)
%                         [assnmt,nlcost]= mbestwrap_updt_joint(neglogcostm,hbes); rankscosts= exp(-nlcost);
                        
                        assnmt=assnmt-ntotaltracks; assnmt(assnmt<=0)= assnmt(assnmt<=0)-1; %set not born/not survived states to negative assignment
                        %meas update/clutter update for tracks
                        for hidx=1:min(hbes,length(rankscosts))    
                            nupdatetracks= sum(assnmt(hidx,:)>0);
                            if nupdatetracks <= N_max
                                nuidx= nupdatetracks+1;
                                comptpos= length(hyps_temp{nuidx})+1;
                                trackcounter=1;
                                for tidx=1:ntotaltracks
                                    asstmp= assnmt(hidx,tidx);
                                    if asstmp > 0
                                        %offset of current track is position in predicted track table
                                        if tidx > nbirthtracks
                                            newoffset= nbirthtracks+hyps_update{nidx}{cidx}{tidx-nbirthtracks};
                                        else
                                            newoffset= tidx;
                                        end
                                        
                                        %index of corresponding updated track
                                        if asstmp > ntotaltracks    %measurement assignment
                                            linidx= length(tracks_predict)*(asstmp-ntotaltracks)+newoffset;
                                        elseif asstmp > 0           %missed detection
                                            linidx= newoffset;
                                        end
                                        hyps_temp{nuidx}{comptpos}{trackcounter}= linidx; trackcounter=trackcounter+1;
                                        tracks_temp{linidx}.used = 1;
                                    end
                                end
                            end

                                
                            if nupdatetracks == 0
                                log_wtv_temp{nuidx}= logsumexp([log_wtv_temp{nuidx},log_cdn_update(nidx)+sum(log_wtv_update{nidx}(cidx))+(-model.lambda_c)+sum(log(QSvec))+sum(log(QDvec))+m*(model.log_lambda_c+model.log_clutterpdf)+(-nlcost(hidx))]);
                            elseif nupdatetracks <= N_max
                                log_wtv_temp{nuidx}= cat(1,log_wtv_temp{nuidx},log_cdn_update(nidx)+sum(log_wtv_update{nidx}(cidx))+(-model.lambda_c)+sum(log(QSvec))+sum(log(QDvec))+m*(model.log_lambda_c+model.log_clutterpdf)+(-nlcost(hidx)));
                            end
                                
                        end
                    end
                end               
        end


    
    %calc cardinality distribution and normalize component weights
    for n= 0:N_max
        nidx= n+1;
        log_wtv_temp{nidx}= log_wtv_temp{nidx}+eps(0);
        if isempty(log_wtv_temp{nidx}), log_cdn_temp(nidx)=-inf; else log_cdn_temp(nidx)= logsumexp(log_wtv_temp{nidx}); end
        log_wtv_temp{nidx}= log_wtv_temp{nidx}-logsumexp(log_wtv_temp{nidx});
    end
    log_wtv_temp{1}=0;
    wtv_temp= cell(N_max+1,1); for n=0:N_max, nidx= n+1; wtv_temp{nidx}= exp(log_wtv_temp{nidx}); end
    
    log_cdn_temp= log_cdn_temp - logsumexp(log_cdn_temp); cdn_temp= exp(log_cdn_temp); 
    cdn_temp= cdn_temp/sum(cdn_temp);
    
    cdn_update= cdn_temp; log_cdn_update= log_cdn_temp;
 
    tracks_update= tracks_temp; hyps_update= hyps_temp; log_wtv_update= log_wtv_temp; wtv_update= wtv_temp;
    
    
    %--- merge duplicate hypotheses from posterior
    %flatten posterior into component table
    hcell= cell(0,0);
    tcell= cell(0,0);
    nvect= zeros(0,0);
    wvect= zeros(0,0);
    hidx=1;
    for n=0:N_max
        nidx=n+1;
        if n==0
            hcell{hidx}= sprintf('%i*',[]);
            tcell{hidx}= {};
            nvect(hidx)= n;
            wvect(hidx)= log_wtv_update{nidx};
            hidx= hidx+1;
        else
            for cidx= 1:length(hyps_update{nidx});
                trackpointerstemp= zeros(n,1);
                for tidx=1:n
                    trackpointerstemp(tidx)= hyps_update{nidx}{cidx}{tidx};
                end
                hcell{hidx}= sprintf('%i*',sort(trackpointerstemp(:)'));
                tcell{hidx}= hyps_update{nidx}{cidx};
                nvect(hidx)= n;
                wvect(hidx)= log_wtv_update{nidx}(cidx);
                hidx= hidx+1;
            end
        end
    end
    
    %find unique components and preallocate memory
    [nc,ia,ic]= unique(hcell);
    hyps_temp= cell(N_max+1,1); log_wtv_temp= cell(N_max+1,1); wtv_temp= cell(N_max+1,1);
    for n=0:N_max
        nidx= n+1;
        ncomps= sum(nvect(ia)==n);
        if ncomps ~=0 
            log_wtv_temp{nidx}= NaN*zeros(ncomps,1); wtv_temp{nidx}= NaN*zeros(ncomps,1);
        end
    end
    %write new components
    for hidx= 1:length(nc)
        write_idx_n= length(tcell{ia(hidx)})+1;
        write_idx_c= find(isnan(log_wtv_temp{write_idx_n}),1);
        if ~isempty(tcell{ia(hidx)}), hyps_temp{write_idx_n}{write_idx_c}= tcell{ia(hidx)}; end
        log_wtv_temp{write_idx_n}(write_idx_c,1)= logsumexp(wvect(ic==hidx));
        wtv_temp{write_idx_n}(write_idx_c,1)= exp(log_wtv_temp{write_idx_n}(write_idx_c,1));
    end
    hyps_update= hyps_temp; log_wtv_update= log_wtv_temp; wtv_update= wtv_temp;
    
    %--- compact track table and reindex hypthesis cell
    tracks_temp= cell(0,1); hyps_temp= cell(N_max+1,1);
    trackcount= 0;
    for tabidx= 1:length(tracks_update)
        if tracks_update{tabidx}.used == 1
            trackcount= trackcount+1;
            tracks_update{tabidx}.newidx= trackcount;
            tracks_temp{trackcount,1}= tracks_update{tabidx};
        end
    end
    
    for n=0:N_max
        nidx= n+1;
        for cidx=1:length(hyps_update{nidx})
            for tidx= 1:n
                hyps_temp{nidx}{cidx}{tidx}= tracks_update{hyps_update{nidx}{cidx}{tidx}}.newidx;
            end
        end
    end
      
    for tabidx= 1:length(tracks_temp)
        tracks_temp{tabidx}= rmfield(tracks_temp{tabidx},{'used','newidx','qz'});
    end
    
    tracks_update= tracks_temp; hyps_update= hyps_temp;

    %--- chopping and truncating
 
    %rank hypothesis by weight and store cardinality and birth labels
    comparstack= [];
    weightstack= [];
    nlabelstack= [];
    clabelstack= [];
    for n= 1:N_max
        nidx= n+1;
        numcmp= length(hyps_update{nidx});
        comparstack= cat(1, comparstack, cdn_update(nidx)*wtv_update{nidx});
        weightstack= cat(1, weightstack, wtv_update{nidx});
        nlabelstack= cat(1, nlabelstack, nidx*ones(numcmp,1));
        clabelstack= cat(1, clabelstack, [1:numcmp]');
    end
    totcompraw= length(clabelstack);
    
    [idxkeep]= find(comparstack > chop_threshold);
    comparstack= comparstack(idxkeep);
    weightstack= weightstack(idxkeep); 
    nlabelstack= nlabelstack(idxkeep); 
    clabelstack= clabelstack(idxkeep);
    
    [~,idxsort]= sort(-comparstack); totcomp= length(comparstack); idxsort= idxsort(1:min(totcomp,Hbes));
    comparstack= comparstack(idxsort); weightstack= weightstack(idxsort); nlabelstack= nlabelstack(idxsort); clabelstack= clabelstack(idxsort);
    
    wtv_temp= cell(N_max+1,1); hyps_temp= cell(N_max+1,1);
    
    for u=1:min(totcomp,Hbes) %copy best hypotheses into temp variables
        nidx= nlabelstack(u);
        cidx= clabelstack(u);
        
        wtv_temp{nidx}= cat(1,wtv_temp{nidx},wtv_update{nidx}(cidx)); 
        hyps_temp{nidx}{end+1}= hyps_update{nidx}{cidx};
    end 
    wtv_temp{1}= wtv_update{1}; hyps_temp{1}= hyps_update{1}; %copy zero cardinality
    
    for n= 1:N_max %enforce minimum number per cardinality
        nidx= n+1;
        if length(wtv_temp{nidx}) < Cmin && length(wtv_update{nidx}) ~= 0
            [~,idxcomp]= sort(-wtv_update{nidx});
            idxcomp= idxcomp(1:min(Cmin,length(wtv_update{nidx})));
            for cidx=1:length(idxcomp)
                wtv_temp{nidx}(cidx)= wtv_update{nidx}(idxcomp(cidx));
                hyps_temp{nidx}{cidx}= hyps_update{nidx}{idxcomp(cidx)};
            end
        end
        wtv_temp{nidx}= wtv_temp{nidx}/sum(wtv_temp{nidx});
    end
    wtv_update= wtv_temp; hyps_update= hyps_temp;
    for n=0:N_max 
        nidx=n+1; 
        log_wtv_update{nidx}= log(wtv_update{nidx}); 
        cdn_update(nidx)= cdn_update(nidx)*sum(wtv_update{nidx}); %resets the cdn to zero if all components in this cardinality were truncated 
        log_cdn_update(nidx)= log(cdn_update(nidx));
    end
    %--- state extraction

    %use MAP estimate for cardinality dn and pick highest weights
    [~,mode] = max(cdn_update);
    hat_N_MAP(k) = mode-1;
    hat_N(k) = hat_N_MAP(k);
    [~,idx1]= max(wtv_update{mode});
    hat_X{k}= []; 
    hat_T{k}= [];
    for n=1:hat_N(k)
        [~,idx2]= max(tracks_update{hyps_update{mode}{idx1}{n}}.w);
        
        est_m= tracks_update{hyps_update{mode}{idx1}{n}}.m(:,idx2);
        est_l= tracks_update{hyps_update{mode}{idx1}{n}}.l(3);
        
        hat_X{k} = [hat_X{k}, est_m];
        hat_T{k} = [hat_T{k}, est_l];
    end
    

    %compute diagnoistics for cardinality distribution at update
    [~,mode] = max(cdn_update);
    cdn_update_mean(k) = sum(nvector .* cdn_update);
    cdn_update_mode(k) = mode-1;
    cdn_update_var(k)  = sum(nvector.^2 .* cdn_update) - cdn_update_mean(k)^2;
%          ' N(k)=',num2str(N_true(k)), ...
    if ~strcmp(model.run_flag,'silence')
        disp(['Updt:',...
         ' time= ',num2str(k),...
         ' EN(k)=',num2str(cdn_update_mean(k),4),...
         ' MODN(k)=',num2str(cdn_update_mode(k),4),...
         ' VARN(k)=',num2str(cdn_update_var(k),4),...
         ' #trak updt=' num2str(length(tracks_update)),...
         ' #comp updt=',num2str(totcompraw),...
         ' #comp merg=',num2str(min(totcomp,Hbes))   ]);
    end

    %--- stats
    cpu_time_lg(k)= cputime-time_start;

    %--- store
    cdn_update_stack(:,k) = cdn_update; 
%     if strcmp(run_flag,'disp')
%        imshow(img);
%        hold on;
%        draw_box_tr(k,img,hat_X{k},hat_T{k});        
%     end
    
    if ~isempty(hat_X{k})
        W = sqrt(hat_X{k}(3,:).*hat_X{k}(4,:));
        H = hat_X{k}(3,:)./W;
        X = hat_X{k}(1,:) - W/2;
        Y = hat_X{k}(2,:) - H/2;
        id = hat_T{k};
    else
        W=[]; H=[]; X=[]; Y=[]; id=[];
    end
    
    track.w = [W'; track.w];
    track.h = [H'; track.h];
    track.xhat = [X'; track.x];
    track.yhat = [Y'; track.y];
    track.fr = [k'; track.fr];
    track.id = [id'; track.id];
    
    
end

if strcmp(model.disp_flag,'vis')
    for k=1:K
        try
            if strcmp(other_param.seq,'PETS2009') || strcmp(other_param.seq,'KITTI_train') || strcmp(other_param.seq,'KITTI_test')
                img = imread(sprintf(input_frames, k-1+offset)); %% read an image
            else
                img = imread(sprintf(input_frames, k+offset)); %% read an image
            end
        catch
            break
        end
        imshow(img);
        hold on;
        draw_box_tr(k,img,hat_X{k},hat_T{k}); 
        pause(.01);
        drawnow
        hold off;
     
    end
  
end

%put the measurements back once we're done
Z_gate = Z; Z = Z_orig;
lambda_c= origclutrate; log_lambda_c= log_origclutrate; clutterpdf= origclutpdf; log_clutterpdf= log_origclutpdf;

end