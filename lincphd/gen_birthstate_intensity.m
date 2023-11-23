function X= gen_birthstate_intensity(model,num_par)
% for Gaussian mixture intensity

L_b= length(model.lambda_b); nx= length(model.bar_x(:,1));
sum_lambda_b= sum(model.lambda_b);
temp= rand(num_par,1);
threshold= 0;
X=[];
for i=1:L_b
    threshold= threshold+ model.lambda_b(i)/sum_lambda_b;
    idx= find(temp > threshold);
    X= [ X repmat(model.bar_x(:,i),[1 length(temp)- length(idx)])+ ...
            model.bar_B(:,:,i)*randn(nx,length(temp)- length(idx)) ];
    temp= temp(idx);
end;


%%% modification note: seems to be called by smc_phd_codes routes
%%% may need modification to account for cardinality dn parameters
%%% missing third input argument? phdfilter code calls with 3rd arg p_death