
% Demo code for VGLMB

addpath(genpath(pwd));
addpath(genpath('\\rmit.internal\USRHome\el0\e75710\MATLAB\contracking-v1.0'));
addpath(genpath('\\rmit.internal\USRHome\el0\e75710\MATLAB\contracking-v1.0\utils'));

setOtherParameters
setPathVariables;

for i = 1                 
    % set a null hypothesis likelihood
    adjustOtherParameters;
    setFilterParameters;
        
    % load detections
    det = loadDet(det_input_path{i}, other_param);
    
    % run VGLMB
    track = VGLMB(det,img_input_path{i},model,other_param);
    % save tracking output into images
    visTracks(track, other_param, img_input_path{i}, img_output_path{i}, max(det.fr));                     
end



        
