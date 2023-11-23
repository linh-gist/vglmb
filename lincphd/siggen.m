%=== signal generation

%--- this one is important; it defines the problem and model
declare_problem
%---

X= cell(K,1); Z= cell(K,1); %state and observation sets
N_true= zeros(K,1);         %true no. of targets
track_list= cell(K,1);         %record the track number list at time k
total_tracks= 0;
%--- time step 1
X_birth= gen_birthstate(model);
N_b= size(X_birth,2);
if N_b > 0,
    X{1}= X_birth;
else    %enforce 1 birth
    X{1}= model.bar_x(:,1)+ model.bar_B(:,:,1)*randn(x_dim,1);
    N_b= 1;
end;
track_list{1}= 1:N_b; total_tracks= N_b; N_true(1)= N_b;
idx= find( rand(N_b,1) <= P_D );
Z{1}= gen_observation(model,X{1}(:,idx));
N_c= poissrnd(lambda_c);    %no. of clutter points
C= repmat(range_c(:,1),[1 N_c])+ diag(range_c*[ -1; 1 ])*rand(z_dim,N_c);  %clutter generation
Z{1}= [ Z{1} C ];

%--- time step 2 to K
for k=2:K
    %target death/survival
    for i=1:N_true(k-1)
        if rand < P_S,
            x_temp= gen_newstate(model,X{k-1}(:,i));
            z_temp= model.C_posn*x_temp;
            if (sum(z_temp >= range_c(:,1))+ sum(z_temp<= range_c(:,2))) == 2*z_dim, %in range?
                X{k}= [ X{k} x_temp ];
                N_true(k)= N_true(k)+ 1;
                track_list{k}= [ track_list{k} track_list{k-1}(i) ];
            end;
        end;
    end;
    %target birth
    if N_true(k) < M_max,
        X_birth= gen_birthstate(model);
        if size(X_birth,2) > M_max- N_true(k),
            X_birth= X_birth(:,M_max- N_true(k));
        end;
        X{k}= [ X{k} X_birth ];
        N_b= size(X_birth,2);
        N_true(k)= N_true(k)+ N_b;
        track_list{k}= [ track_list{k} total_tracks+(1:N_b) ];
        total_tracks= total_tracks+N_b;
    end;
    
    %cphd filter does not cater for spawning

    %observation
    if N_true(k)> 0,
        idx= find( rand(N_true(k),1) <= P_D );
        Z{k}= gen_observation(model,X{k}(:,idx));
    end;
    N_c= randsample(0:N_max,1,true,cdn_clutter);    %no. of clutter points
    C= repmat(range_c(:,1),[1 N_c])+ diag(range_c*[ -1; 1 ])*rand(z_dim,N_c);  %clutter generation
    Z{k}= [ Z{k} C ];
end;
    
