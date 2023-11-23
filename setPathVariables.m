Sequences = {'PETS2009','MOT_Challenge_train','MOT_Challenge_test'};

det_input_dir = {'input/PETS2009/','input/MOT_Challenge/train/','input/MOT_Challenge/test/'};
                    
% set your output directory here
img_output_dir = {'/output/temp_img/','/output/temp_img/',...
                        '/output/temp_img/'};

% set your dataset directory here                       
img_input_dir = {'//rmit.internal/USRHome/el0/e75710/Configuration/Desktop/vglmb_code/external/devkit/data/Crowd_PETS09/', '//rmit.internal/USRHome/el0/e75710/MATLAB/devkit/data/2DMOT2015/train/',...
                        '//rmit.internal/USRHome/el0/e75710/MATLAB/devkit/data/2DMOT2015/test/'};

% set your image path here. The entire image path should be [img_input_dir img_input_subdir].                   
img_input_subdir = {'S2/L1/Time_12-34/View_001/frame_%04d.jpg', 'S2/L2/Time_14-55/View_001/frame_%04d.jpg', 'S2/L3/Time_14-41/View_001/frame_%04d.jpg',...
                        'S1/L1/Time_13-59/View_001/frame_%04d.jpg', 'S1/L2/Time_14-06/View_001/frame_%04d.jpg',...
                        'ADL-Rundle-6/img1/%06d.jpg', 'ADL-Rundle-8/img1/%06d.jpg', 'ETH-Bahnhof/img1/%06d.jpg', 'ETH-Pedcross2/img1/%06d.jpg',...
                        'ETH-Sunnyday/img1/%06d.jpg', 'KITTI-13/img1/%06d.jpg', 'KITTI-17/img1/%06d.jpg', 'PETS09-S2L1/img1/%06d.jpg', 'TUD-Campus/img1/%06d.jpg', 'TUD-Stadtmitte/img1/%06d.jpg', 'Venice-2/img1/%06d.jpg',...
                        'ADL-Rundle-1/img1/%06d.jpg','ADL-Rundle-3/img1/%06d.jpg','AVG-TownCentre/img1/%06d.jpg','ETH-Crossing/img1/%06d.jpg','ETH-Jelmoli/img1/%06d.jpg','ETH-Linthescher/img1/%06d.jpg',...
                        'KITTI-16/img1/%06d.jpg','KITTI-19/img1/%06d.jpg','PETS09-S2L2/img1/%06d.jpg','TUD-Crossing/img1/%06d.jpg','Venice-1/img1/%06d.jpg'};
                    
det_input_name = {'PETS2009-S2L1-c1-app_pca','PETS2009-S2L2-c1-app_pca','PETS2009-S2L3-c1-app_pca','PETS2009-S1L1-2-c1-app_pca','PETS2009-S1L2-1-c1-app_pca',...
                    'ADL-Rundle-6', 'ADL-Rundle-8', 'ETH-Bahnhof', 'ETH-Pedcross2', 'ETH-Sunnyday', 'KITTI-13', 'KITTI-17', 'PETS09-S2L1', 'TUD-Campus', 'TUD-Stadtmitte', 'Venice-2',...
                    'ADL-Rundle-1','ADL-Rundle-3','AVG-TownCentre','ETH-Crossing','ETH-Jelmoli','ETH-Linthescher','KITTI-16','KITTI-19','PETS09-S2L2','TUD-Crossing','Venice-1'};         


% find an index of the query sequence
seq_idx = 1;
while 1 
    if strcmp(other_param.seq,Sequences{seq_idx})        
        break;
    end
    seq_idx = seq_idx + 1;
end

% select input image indices based on the query sequence
input_idx = 0;
switch seq_idx
    case 1       
        input_idx = 1:5;
    case 2
        input_idx = 6:16;
    case 3
        input_idx = 17:27;
%     case 4
%         input_idx = 38:48;
%     case 5
%         input_idx = 49:77;
    otherwise
        error('unexpected sequence index');    
end

% set the input detection path 
det_input_path = cell(1,length(input_idx));
for i = 1:length(input_idx)
    det_input_path{i} = [det_input_dir{seq_idx} det_input_name{input_idx(i)} '.mat'];
end

% set the input/output imgage path 
img_output_path = cell(1,length(input_idx));
img_input_path = cell(1,length(input_idx));
for i = 1:length(input_idx)
    if isempty(other_param.appSel)        
        img_output_path{i} = [img_output_dir{seq_idx} Sequences{seq_idx} '/' det_input_name{input_idx(i)} '/mot/'];
    else
        img_output_path{i} = [img_output_dir{seq_idx} Sequences{seq_idx} '/' det_input_name{input_idx(i)} '/app/'];
    end        
    img_input_path{i} = [img_input_dir{seq_idx} img_input_subdir{input_idx(i)}];
end

% load camera parameters when PETS is selected
if strcmp(other_param.seq,'PETS2009')
    % get the camera parameters
    load input/PETS2009/PETS2009_S2L1_camera_parameters.mat;    
    other_param.camParam = camParam;
end