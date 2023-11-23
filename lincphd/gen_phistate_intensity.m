function X= gen_birthstate_intensity(model,old_X,P_S)
% for Gaussian mixture intensity

L_s= length(model.lambda_s); 
[nx,num_par]= size(old_X);
sum_term= sum(model.lambda_s)+ P_S;

temp= rand(num_par,1);
threshold= 0;
X= zeros(nx,num_par);
%spawn
for i=1:L_s
    threshold= threshold+ model.lambda_s(i)/sum_term;
    idx= find(temp <= threshold);
    X(:,idx)= repmat(model.chk_x(:,i),[1 length(idx)])+ ...
        model.chk_A(:,:,i)*old_X(:,idx)+ ...
        model.sigma_s(i)*model.chk_B(:,:,i)*randn(size(model.chk_B(:,:,i),2),length(idx));
    temp(idx)= 1.1;
end;
%survived
idx= find(temp <= 1);
X(:,idx)= model.A*old_X(:,idx)+ model.sigma_v*model.B*randn(size(model.B,2),length(idx));

