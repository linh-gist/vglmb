function [assignments,costs]= mbestwrap(P0,m)
% MBEST wrapper for Murty's algorithm
% used for multi-target tracking with clutter and missed detections

n1= size(P0,1);
n2= size(P0,2);

%padding blocks for dummy variables
blk1= -log(ones(n1,n1));
blk2= -log(ones(n2,n2));
blkr= -log(ones(n2,n1));

% blk1= -log(eye(n1));
% blk2= -log(eye(n2));
% blkr= -log(zeros(n2,n1));

P0= [P0 blk1; blk2 blkr];

%murty
[assignments, costs]= murty(P0,m); %or use murty_extra if alternative implementation is desired

%strip dummy variables
assignments= assignments(:,1:n1);

%dummy assignments are births
assignments(assignments>n2)= 0;

%keep only unique solutions for original variables
[b,i,j]= unique(assignments,'rows');
assignments= assignments(i,:);
costs= costs(i);

%keep only unique combinations (prediction only!)
assignments= sort(assignments,2,'ascend');
[b,i,j]= unique(assignments,'rows');
assignments= assignments(i,:);
costs= costs(i);
