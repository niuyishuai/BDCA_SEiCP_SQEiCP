function xopt = simplexproj(y)
% check feasibility
if (abs(sum(y)-1)<1e-8 && all(y>=0))
    xopt=y;
    return;
end
% n=1
n=length(y);
if (n==1)
    xopt=1;
    return;
end
% find projection
v = sort(y,'descend');
s=v(1);
bestmu=v(1)-1;
for i=2:n
    s=s+v(i);
    mu=(s-1)/i;
    if mu>=v(i)
        xopt = max(y-bestmu,0);
        return;
    end
    bestmu=mu;
end
xopt = max(y-bestmu,0);
end