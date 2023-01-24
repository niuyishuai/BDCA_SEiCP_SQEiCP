function [fval,dfval] = fobj_eval_QP(x,A,opt)
if opt==1 % evaluate function only
    fval = -x'*A*x;
else % evaluate function and gradient
    Ax = A*x;
    fval = -x'*Ax;
    dfval = -2*Ax;
end
end