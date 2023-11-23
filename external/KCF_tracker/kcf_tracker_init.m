% function [positions, time] = kcf_tracker(video_path, img_files, pos, target_sz, ...
% 	padding, kernel, lambda, output_sigma_factor, interp_factor, cell_size, ...
% 	features, show_visualization)

function ref = kcf_tracker_init(img, pos, target_sz, ...
	padding, kernel, lambda, output_sigma_factor, interp_factor, cell_size, ...
	features)


%TRACKER Kernelized/Dual Correlation Filter (KCF/DCF) tracking.
%   This function implements the pipeline for tracking with the KCF (by
%   choosing a non-linear kernel) and DCF (by choosing a linear kernel).
%
%   It is meant to be called by the interface function RUN_TRACKER, which
%   sets up the parameters and loads the video information.
%
%   Parameters:
%     VIDEO_PATH is the location of the image files (must end with a slash
%      '/' or '\').
%     IMG_FILES is a cell array of image file names.
%     POS and TARGET_SZ are the initial position and size of the target
%      (both in format [rows, columns]).
%     PADDING is the additional tracked region, for context, relative to 
%      the target size.
%     KERNEL is a struct describing the kernel. The field TYPE must be one
%      of 'gaussian', 'polynomial' or 'linear'. The optional fields SIGMA,
%      POLY_A and POLY_B are the parameters for the Gaussian and Polynomial
%      kernels.
%     OUTPUT_SIGMA_FACTOR is the spatial bandwidth of the regression
%      target, relative to the target size.
%     INTERP_FACTOR is the adaptation rate of the tracker.
%     CELL_SIZE is the number of pixels per cell (must be 1 if using raw
%      pixels).
%     FEATURES is a struct describing the used features (see GET_FEATURES).
%     SHOW_VISUALIZATION will show an interactive video if set to true.
%
%   Outputs:
%    POSITIONS is an Nx2 matrix of target positions over time (in the
%     format [rows, columns]).
%    TIME is the tracker execution time, without video loading/rendering.
%
%   Joao F. Henriques, 2014


	%if the target is large, lower the resolution, we don't need that much
	%detail
	resize_image = (sqrt(prod(target_sz)) >= 100);  %diagonal size >= threshold
	if resize_image,
		pos = floor(pos / 2);
		target_sz = floor(target_sz / 2);
	end


	%window size, taking padding into account
	window_sz = floor(target_sz * (1 + padding));
	
% 	%we could choose a size that is a power of two, for better FFT
% 	%performance. in practice it is slower, due to the larger window size.
% 	window_sz = 2 .^ nextpow2(window_sz);

	
	%create regression labels, gaussian shaped, with a bandwidth
	%proportional to target size
	output_sigma = sqrt(prod(target_sz)) * output_sigma_factor / cell_size;
	yf = fft2(gaussian_shaped_labels(output_sigma, floor(window_sz / cell_size)));

	%store pre-computed cosine window
	cos_window = hann(size(yf,1)) * hann(size(yf,2))';	

	
	%note: variables ending with 'f' are in the Fourier domain.
    im = img;
    if size(im,3)>1
        im = rgb2gray(im);
    end
    if resize_image
        im = imresize(im, 0.5);
    end
    
    %obtain a subwindow for training at newly estimated target position
    patch = get_subwindow(im, pos, window_sz);
	xf = fft2(get_features(patch, features, cell_size, cos_window));

	%Kernel Ridge Regression, calculate alphas (in Fourier domain)
	switch kernel.type
        case 'gaussian',
			kf = gaussian_correlation(xf, xf, kernel.sigma);
		case 'polynomial',
			kf = polynomial_correlation(xf, xf, kernel.poly_a, kernel.poly_b);
		case 'linear',
			kf = linear_correlation(xf, xf);
		end
	alphaf = yf ./ (kf + lambda);   %equation for fast training
    
    ref.alphaf = alphaf;
    ref.xf = xf;
    ref.size = target_sz;
%  
% 
% 	if resize_image,
% 		positions = positions * 2;
% 	end
end

