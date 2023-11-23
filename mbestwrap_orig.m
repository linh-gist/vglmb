function [assignments,costs]= mbestwrap(P0,m)
% MBEST wrapper for Murty's algorithm
% used for multi-target tracking with clutter and missed detections

n1= size(P0,1);
n2= size(P0,2);

blk1= -log(ones(n1,n1));
blk2= -log(ones(n2,n2));
blkr= -log(ones(n2,n1));

% blk1= -log(eye(n1));
% blk2= -log(eye(n2));
% blkr= -log(zeros(n2,n1));

P0= [P0 blk1; blk2 blkr];

[assignments, costs]= murty(P0,m); %or use murty_extra if alternative implementation is desired

assignments= assignments(:,1:n1);

assignments(assignments>n2)= 0;

[b,i,j]= unique(assignments,'rows');

assignments= assignments(i,:);
costs= costs(i);

