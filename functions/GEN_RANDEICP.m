%% generate RANDEICP

%randbnd = [-1,1];
randbnd = [-10,10];

for n=[50,100,200,400,600,800]
A=rand(n)*(randbnd(2)-randbnd(1))+randbnd(1);
A=(A'+A)/2;
B=eye(n);
varmu=sdpvar(1);
optimize([A+varmu*B>=0],varmu,sdpsettings('solver','mosek','verbose',0));
mu=value(varmu)+1;
A = A + mu*B;
x0 = rand(n,1);
%x0 = x0/sqrt(x0'*B*x0);
filename = sprintf('RANDEICP_%d_%d_%d',randbnd(1),randbnd(2),n);
save(['RANDEICP/',filename],'n','A','B','mu','x0');
end