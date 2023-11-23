%=== signal 're-generation' using the same X



Z= cell(K,1); %state and observation sets

for k=1:K
    if N_true(k)> 0,
        tmp_P_D_vec= P_D*ones(N_true(k),1); if any(tbirth==k), tmp_P_D_vec(find(track_list{k}==(find(tbirth==k))))=1; end
        idx= find( rand(N_true(k),1) <= tmp_P_D_vec ); if length(idx) ~= N_true(k), disp([num2str(k) ' MD']); end
        Z{k}= gen_observation(model,X{k}(:,idx));
    end;
    N_c= poissrnd(lambda_c);    %no. of clutter points
    C= repmat(range_c(:,1),[1 N_c])+ diag(range_c*[ -1; 1 ])*rand(z_dim,N_c);  %clutter generation
    Z{k}= [ Z{k} C ];
end;
    
