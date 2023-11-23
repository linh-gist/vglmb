function X= gen_gms(w,m,P,num_par)
% for Gaussian mixture intensity

x_dim=size(m,1);
X=zeros(x_dim,num_par);

w= w/sum(w);
w= sort(w,'descend');

nc= length(w);

comps= randsample(1:nc,num_par,true,w);

ns= zeros(nc,1);
for c=1:nc
    ns(c)= nnz(comps==c);
end

startpt= 1;
for i=1:nc
    endpt= startpt+ns(i)-1;
    X(:,startpt:endpt)= mvnrnd(m(:,i)',P(:,:,i),ns(i))';
    startpt= endpt+1;
end


