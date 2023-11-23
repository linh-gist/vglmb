function [assignments,costs]= mbestwrap_updt_gibbsamp(P0,m)
    
n1 = size(P0,1);
n2 = size(P0,2);

assignments= zeros(m,n1);
costs= zeros(m,1)';

%[tempinit,~]= mofo(P0); currsoln= tempinit(:)';
%[tempinit,~]= munkres(P0-min(min(P0))); currsoln= tempinit(:)';
[tempinit,~]= assignmentoptimal(P0-min(min(P0))); currsoln= tempinit(:)';

assignments(1,:)= currsoln;
costs(1)=sum(P0(sub2ind(size(P0),1:n1,currsoln)));

% %do in log domain
% for sol= 2:m
%     for var= 1:n1
%         tempsamp= P0(var,:); %grab row of costs for current association variable
%         tempsamp(currsoln([1:var-1,var+1:end]))= -log(0); %lock out current and previous iteration step assignments except for the one in question
%         tempsamp= -tempsamp; %restore positive weights
%         tempsamp(isfinite(tempsamp))= tempsamp(isfinite(tempsamp))-logsumexp(tempsamp(isfinite(tempsamp))); %normalize weights
%         currsoln(var)= resample(exp(tempsamp),1); %generate next sample from conditional
%     end
%     assignments(sol,:)= currsoln;
%     costs(sol)= sum(P0(sub2ind(size(P0),1:n1,currsoln)));
% end

%do in real domain
for sol= 2:m
    for var= 1:n1
        tempsamp= exp(-P0(var,:)); %grab row of costs for current association variable
        tempsamp(currsoln([1:var-1,var+1:end]))= 0; %lock out current and previous iteration step assignments except for the one in question
        [~,currsoln(var)]= histc(rand(1,1),[0;cumsum(tempsamp(:))/sum(tempsamp)]); %currsoln(var)= resample(tempsamp/sum(tempsamp),1);
    end
    assignments(sol,:)= currsoln;
    costs(sol)= sum(P0(sub2ind(size(P0),1:n1,currsoln)));
end

[C,I,~]= unique(assignments,'rows','stable');
assignments= C;
costs= costs(I);
