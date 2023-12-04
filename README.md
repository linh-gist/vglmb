# vglmb
Source code for the paper "**A labeled random finite set online multi-object tracker for video data**".

_Kim, D. Y., Vo, B. N., Vo, B. T., & Jeon, M. (2019). A labeled random finite set online multi-object tracker for video data. Pattern Recognition, 90, 377-389._

# Detection
FairMOT detector with format **_{frame}, {top},{left},{bottom},{right},{confidence score}_**, frame index starting from zero

**MOT15**: _mot15_fairmot128_conf0.3_tlbr_, loaded "fairmot_dla34.pth" pre-trained weight, obtain detection with a confidence score higher than 0.3

**MOT16**: _mot16_fairmot128_conf0.5_tlbr_, loaded "fairmot_dla34.pth" pre-trained weight, obtain detection with a confidence score higher than 0.5


_Zhang, Y., Wang, C., Wang, X., Zeng, W., & Liu, W. (2021). FairMOT: On the fairness of detection and re-identification in multiple object tracking. International Journal of Computer Vision, 129, 3069-3087._


# Result
Saved with MOTChallenge format **_{frame}, {id}, {bb_left}, {bb_top}, {bb_width}, {bb_height}, 1,-1,-1,-1_**
