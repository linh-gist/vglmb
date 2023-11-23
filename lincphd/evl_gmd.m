function wts= evl_gmd(X,w,m,P)
% for Gaussian mixture intensity

x_dim=size(m,1);

nc= length(w);
ns= size(X,2);

wts= zeros(ns,nc);

for i=1:nc
    wts(:,i)= w(i)*mvnpdf(X',m(:,i)',P(:,:,i))';
end

wts= sum(wts,2);
