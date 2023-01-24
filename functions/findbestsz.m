function bestz = findbestsz(coef,ubz)
% coef = [a1,b1,c1,a2,b2,c2]; ubz>0.
a1=coef(1);
b1=coef(2);
c1=coef(3);
a2=coef(4);
b2=coef(5);
c2=coef(6);

rr=roots([a1*b2-a2*b1,2*(a1*c2-a2*c1),b1*c2-b2*c1]);
r1=rr(rr>=0);
r2=r1(r1<=ubz);
L=[0,ubz,r2];
nom=polyval([a1,b1,c1],L);
denom=polyval([a2,b2,c2],L);
[~,idx]=min(nom./denom);
bestz=min(L(idx));
%bestz=max(L(idx)); % for largest stepsize
end