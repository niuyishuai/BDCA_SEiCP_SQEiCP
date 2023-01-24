function [fval,dfval] = fobj_eval_LnP(x,A,B,opt)
if opt==1 % evaluate function only
    fval = log((x'*B*x)/(x'*A*x));
else % evaluate function and gradient
    Bx=B*x;
    Ax=A*x;
    xBx=x'*Bx;
    xAx=x'*Ax;
    fval = log(xBx/xAx);
    dfval = 2*((Bx/xBx) - (Ax)/(xAx));
end
end