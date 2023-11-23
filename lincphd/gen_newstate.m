function X= gen_newstate(model,X_old)
%--- this new state generation function follows the linear state spc eqn.
% x= Ax_old + Bv

if isempty(X_old),
    X= [];
else
    X= model.A*X_old+ model.sigma_v*model.B*randn(size(model.B,2),size(X_old,2));
end;