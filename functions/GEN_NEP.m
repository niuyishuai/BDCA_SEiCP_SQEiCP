NEPlst=dir('NEP\\*.mtx');
for i=1:length(NEPlst)
    fname=NEPlst(i).name;
    fprintf('Loading model from %s\n',fname);
    [A,n]=mmread(sprintf('NEP\\%s',fname));
    A=(A+A')/2;
    B=eye(n);
    % compute best mu by solving an SDP
    varmu=sdpvar(1);
    optimize(A+varmu*B>=0,varmu,sdpsettings('solver','mosek','verbose',0));
    mu=value(varmu)+1;
    A = A + mu*B;
    x0 = rand(n,1);
    save(sprintf('NEP\\NEP_%s.mat',fname),'n','A','B','mu','x0');
end