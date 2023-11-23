

base_path = '\\rmit.internal\USRHome\el0\e75710\MATLAB\devkit\data\Crowd_PETS09\S2\L1\Time_12-34\View_001\';

video_path = base_path;
img_files = dir([video_path '*.jpg']);

% [img_files, pos, target_sz, ground_truth, video_path] = load_video_info(base_path, video);


		
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
		
%call tracker function with all the relevant parameters
[positions, time] = tracker(video_path, img_files, pos, target_sz, ...
    padding, kernel, lambda, output_sigma_factor, interp_factor, ...
    cell_size, features, show_visualization);