
%compute IUO between two bboxes in the form [x1,y1,x2,y2]

function o = iou(bb_test,bb_gt)

xx1 = max(bb_test(1), bb_gt(1));
yy1 = max(bb_test(2), bb_gt(2));
xx2 = min(bb_test(3), bb_gt(3));
yy2 = min(bb_test(4), bb_gt(4));
w = max(0, xx2 - xx1);
h = max(0, yy2 - yy1);
wh = w * h;
% o = wh / ((bb_test(3)-bb_test(1))*(bb_test(4)-bb_test(2)) * (bb_gt(3)-bb_gt(1)) * (bb_gt(4)-bb_gt(2)) - wh);

o = wh / ((bb_test(3)-bb_test(1))*(bb_test(4)-bb_test(2))+ (bb_gt(3)-bb_gt(1))*(bb_gt(4)-bb_gt(2)) - wh);