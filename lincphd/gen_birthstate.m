function X= gen_birthstate(model)
% this initial target state generation function uses a Gaussian mixture dn
%  sum lambda_b(j)*N(x;bar_x(j), bar_B(j)*bar_B(j)')
% j=1:J
% and respective discrete cardinality distribution model.cdn_b

L_b= length(model.lambda_b); nx= length(model.bar_x(:,1));
N_max = length(model.cdn_b)-1;

nbirths = randsample(0:N_max,1,true,model.cdn_b);
if nbirths == 0
    X = [];
else
    X = [];
    for targetnum=1:nbirths
        mixtureid = randsample(1:L_b,1,true,model.lambda_b);
        X = [X model.bar_x(:,mixtureid)+model.bar_B(:,:,mixtureid)*randn(nx,1)];
    end
end
