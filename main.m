
% Demo code for VGLMB

addpath(genpath(pwd));
addpath(genpath('/Users/duyongkim/Documents/MATLAB/contracking-v1.0/'));
addpath(genpath('/Users/duyongkim/Documents/MATLAB/contracking-v1.0/utils'));
addpath(genpath('/Users/duyongkim/Documents/MATLAB/devkit'));

setOtherParameters
setPathVariables;

% for i = 1:length(det_input_path) 
for i = 11
    % set a null hypothesis likelihood
    adjustOtherParameters;
    setFilterParameters;
        
    % load detections
    det = loadDet(det_input_path{i}, other_param);
    
    % run MHT
%     track = MHT(det, kalman_param, other_param);
    track = VGLMB(det,img_input_path{i},model,other_param);
    % save tracking output into images
    visTracks(track, other_param, img_input_path{i}, img_output_path{i}, max(det.fr));                     
end


% Seq_Name = '/data/Crowd_PETS09/S2/L1/Time_12-34/View_001';
% load 'PETS09_S2L1_View1_det';
% Seq_Name = '/data/MOT16/test/MOT16-07/img1';

% Seq_Name = '/data/MOT17Det/train/MOT17-11/img1';
% Seq_Name = '/data/MOT17Det/test/MOT17-12/img1';
% K = length(detection);
% load '/Users/duyongkim/Documents/MATLAB/devkit/data/MOT16/test/MOT16-07/det/det.txt';
% load '/Users/duyongkim/Documents/MATLAB/devkit/data/MOT17Det/train/MOT17-11/det/det.txt';
% load '/Users/duyongkim/Documents/MATLAB/devkit/data/MOT17Labels/test/MOT17-12-DPM/det/det.txt';
% K = max(det(:,1));
% 
% detection = cell(K,1);
% det_threshold = 0;
% for t=1:max(det(:,1))
%     detection{t}= det(det(:,1)==t & det(:,7) > det_threshold,[3:7])';
% end

% 
% filedirec = ['/Users/duyongkim/Documents/MATLAB/devkit',Seq_Name];
% %length_of_frames = length(dir(fullfile(filedirec, '*png')));
% length_of_frames = length(dir(fullfile(filedirec, '*jpg')));
% 
% %Filtering stage
% declare_problem;
% K = 200;
% bf_lin_logs10;
% % extract_labels1;
%     % clearvars -except K title filedirec Detections hat_X hat_T hat_N idx Z c X tmpl_stack;
% % clearvars -except K filedirec detection hat_X hat_T hat_N idx Z model cpu_time_lg cdn_update_mean cdn_update_mode cdn_update_var;
% clearvars -except K track_count filedirec detection hat_X hat_T hat_N Z Z_birt model cpu_time_lg cdn_update_mean cdn_update_mode cdn_update_var;
% % save_result;
% 
% %Display
% figure(1);
% for f=1:K
% %     FileDirec = [filedirec, '/frame_', repmat('0',1,4-length(num2str(f-1))) num2str(f-1) '.jpg'];
%     FileDirec = [filedirec, '/', repmat('0',1,6-length(num2str(f))) num2str(f) '.jpg'];
%     
%     Img = imread(FileDirec);
% 
%     subplot(211);
%     subaxis(2,1,1,'Spacing',0.03,'Padding',0,'Margin',0);
%     imshow(Img);
%     Det{f} = Z{f}';
%     state= Z{f};    
%     draw_box(f,Img,state);
%     drawnow;
%     hold off;
%     
%     subplot(212);
%     subaxis(2,1,2,'Spacing',0.03,'Padding',0,'Margin',0);
%     imshow(Img);
% %     hold on;
%     draw_box_tr(f,Img,hat_X{f},hat_T{f})
%     pause(.01);
%     drawnow
%     hold off;
% 
% end



        
