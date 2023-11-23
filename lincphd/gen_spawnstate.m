function X= gen_spawnstate(model,x_old)
% this initial target state generation function uses a Gaussian mixture distribution
% N( chk_x+ chk_A*x_old, chk_Q) with average number of targets lambda_s for each
% component

X= [];
L_s= length(model.lambda_s); 
if L_s== 0, return; end;
nx= length(model.chk_x(:,1)); nv= size(model.chk_B(:,:,1),2);
for i=1:L_s
    num= poissrnd(model.lambda_s(i));
    X= [ X repmat(model.chk_x(:,i)+model.chk_A(:,:,i)*x_old,[ 1 num ])+ ...
            model.sigma_s(i)*model.chk_B(:,:,i)*randn(nv,num) ];
end;
