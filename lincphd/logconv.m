function logw = logconv(logu,logv)

m = length(logu);
n = length(logv);
logw = zeros(m+n-1,1);

for k=1:m+n-1
    idxs = max(1,k+1-n):min(k,m);
    logw(k) = logsumexp(logu(idxs)+logv(k+1-idxs));
end