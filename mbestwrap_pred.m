function [assignments,costs]= murty_wrapper_predict(P0,m)

  n1 = size(P0,1);
  n2 = size(P0,2);
      
  % Padding blocks for dummy variables
  blk1 = -log(ones(n1,n1));
  blk2 = -log(ones(n2,n2));
  blkr = -log(ones(n2,n1));

  % blk1= -log(eye(n1));
  % blk2= -log(eye(n2));
  % blkr= -log(zeros(n2,n1));

  P0 = [P0 blk1; blk2 blkr];

  % Make costs non-negative (required by 'assignmentoptimal')
  x = min(min(P0));
  P0 = P0 - x;
      
  % Murty
  [assignments, costs] = murty(P0,m);
      
  % Restore the correct costs to assignments
  costs = costs + (x.*sum(assignments>0,2))';
      
  % Strip dummy variables
  assignments = assignments(:,1:n1);

  % Dummy assignments are births
  assignments(assignments>n2)= 0;

  % Keep only unique solutions for original variables
  [b,i,j] = unique(assignments,'rows');
  assignments = assignments(i,:);
  costs = costs(i);

  % Keep only unique combinations (prediction only!)
  assignments = sort(assignments,2,'ascend');
  [b,i,j] = unique(assignments,'rows');
  assignments = assignments(i,:);
  costs = costs(i);
      
end
