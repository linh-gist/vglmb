function det = NMS(detection)

ov_threshold = 0.4;
% ov_threshold = 0.3;
observation.x = detection(1,:)';
observation.y = detection(2,:)';
observation.w = detection(3,:)';
observation.h = detection(4,:)';
observation.sc = detection(5,:)';

% calculate overlap
obsNo = length(observation.x);
obsInd = 1:obsNo;
overlap_mat = zeros(obsNo,obsNo);
for i = 1:obsNo
    obsIndSel = setdiff(obsInd,i);
    overlap = calc_overlap(observation,i,observation,obsIndSel);
    overlap_mat(i,obsIndSel) = overlap;
end


overlap_mat(overlap_mat < ov_threshold) = 0;
[Matching cost] = Hungarian(-overlap_mat);
[row,col] = find((Matching==1) & (overlap_mat >= ov_threshold));

% when overlapped, pick more confident detection
obsIndDel = [];
for i = 1:length(row)
   if observation.sc(row(i)) > observation.sc(col(i))
       obsIndDel = [obsIndDel; col(i)];
   else
       obsIndDel = [obsIndDel; row(i)];
   end
end

observation.x(obsIndDel) = [];
observation.y(obsIndDel) = [];
% observation.bx(obsIndDel) = [];
% observation.by(obsIndDel) = [];
observation.w(obsIndDel) = [];
observation.h(obsIndDel) = [];
observation.sc(obsIndDel) = [];

obs = [observation.x'; observation.y'; observation.w'; observation.h'; observation.sc'];

det = [obs(1:2,:); obs(1:2,:)+obs(3:4,:);obs(5,:)];
% observation.fr(obsIndDel) = [];
% 
% if other_param.isAppModel
%     observation.app(obsIndDel,:) = [];
% end

