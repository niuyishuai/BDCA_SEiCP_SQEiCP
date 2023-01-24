function [ui1,iter] = FISTA(u0,t0,tol,eta,A,B,Lg)
% FISTA for solving (LnPk)
ti=t0;
yi=u0;
ui=u0;
Bu0=B*u0;
gradh = eta*u0 - (2/(u0'*Bu0))*(Bu0);

Du = inf;
iter=0;
if Lg <= 1e+12 % Lg is not too big
    while Du>tol
        Li=Lg;
        Ayi = A*yi;
        gradphi = eta*yi - (2/(yi'*Ayi))*(Ayi) - gradh;
        ui1 = simplexproj(yi-gradphi/Li);
        ti1 = (1+sqrt(1+4*ti^2))/2;
        du = ui1-ui;
        yi = ui1 + ((ti-1)/ti1)*du; % compute yi+1
        Du = norm(du)/(1+norm(ui1));
        ui=ui1;
        ti=ti1;
        iter=iter+1;
    end
else % Lg is too big
    t=10;
    funcphi=@(x)eta*x'*x/2 - log(x'*A*x) - x'*gradh;
    while Du>tol
        s=1;
        Ayi=A*yi;
        gradphi = eta*yi - 2*Ayi/(yi'*Ayi) - gradh;
        while true % pick a suitable L
            Li=s;
            ui1 = simplexproj(yi-gradphi/Li);
            du = ui1-ui;
            if 0 <= (Li/2)*norm(du)^2 + gradphi'*du + funcphi(yi) - funcphi(ui1)
                break;
            end
            s=s*t;
        end
        ti1 = (1+sqrt(1+4*ti^2))/2;
        yi = ui1 - ((ti-1)/ti1)*du; % compute yi+1
        Du = norm(du)/(1+norm(ui1));
        ui=ui1;
        ti=ti1;
        iter=iter+1;
    end
end
end
