
% Demo code for VGLMB

addpath(genpath(pwd));
addpath(genpath('\\rmit.internal\USRHome\el0\e75710\MATLAB\contracking-v1.0'));
addpath(genpath('\\rmit.internal\USRHome\el0\e75710\MATLAB\contracking-v1.0\utils'));

setOtherParameters
setPathVariables;
% set a null hypothesis likelihood
adjustOtherParameters;
setFilterParameters;

mot15_seqs_train = ["ADL-Rundle-6", "ADL-Rundle-8", "ETH-Bahnhof", "ETH-Pedcross2", "ETH-Sunnyday", "KITTI-13", "KITTI-17", "PETS09-S2L1", "TUD-Campus", "TUD-Stadtmitte", "Venice-2"];
mot15_seqs_test = ["ADL-Rundle-1", "ADL-Rundle-3", "AVG-TownCentre", "ETH-Crossing", "ETH-Jelmoli", "ETH-Linthescher", "KITTI-16", "KITTI-19", "PETS09-S2L2", "TUD-Crossing", "Venice-1"];
mot16_seqs_train = ["MOT16-02", "MOT16-04", "MOT16-05", "MOT16-09", "MOT16-10", "MOT16-11", "MOT16-13"];
mot16_seqs_test = ["MOT16-01", "MOT16-03", "MOT16-06", "MOT16-07", "MOT16-08", "MOT16-12", "MOT16-14"];

root_mot16_train_imgs = "D:/dataset/tracking/mot/MOT16/train/";
root_mot16_test_imgs = "D:/dataset/tracking/mot/MOT16/test/";
mot16_det_dir = "mot16_fairmot128_conf0.5_tlbr/";

root_mot15_train_imgs = "D:/dataset/tracking/mot/2DMOT2015/train/";
root_mot15_test_imgs = "D:/dataset/tracking/mot/2DMOT2015/test/";
mot15_det_dir = "mot15_fairmot128_conf0.3_tlbr/"; 

for i = 1:size(mot15_seqs_train, 2) % mot16_seqs_train
    % run VGLMB
    seq = mot15_seqs_train(i); % mot16_seqs_train
    root_img = strcat(root_mot15_train_imgs, seq, "\img1\"); % root_mot15_train_imgs
    
    img_input_path = dir(root_img);
    img_input_path = string({img_input_path.name});
    img_input_path = sort(img_input_path(1, 3:end));
    img_input_path = root_img + img_input_path;
    
    track = VGLMB(seq,mot15_det_dir,img_input_path,model,other_param); % mot16_det_dir
    % save tracking output into images
    % visTracks(track, other_param, img_input_path{i}, img_output_path{i}, max(det.fr));                     
end



        
