%% well-conditioned sparse RANDQEICP

for n=[50,100,200,400,600]
    for d=[0.05,0.1,0.5,0.7,0.9]
        A=eye(n);
        B=sprandsym(n,0.01);
        C=-matgen_welldefspd(n,0.01);
        D = [A,zeros(n);zeros(n),-C];
        G = [-B,-C;-C,zeros(n)];
        varmu=sdpvar(1);
        optimize([G+varmu*D>=0],varmu,sdpsettings('solver','mosek','verbose',0));
        mu = value(varmu)+1;
        G = G + mu*D;
        x0 = rand(n,1);
        lambda0 = (-x0'*B*x0 + sqrt((x0'*B*x0)^2 - 4*(x0'*A*x0)*(x0'*C*x0)))/(2*x0'*A*x0);
        z0 = [lambda0*x0;x0]/(1+lambda0);
        filename = sprintf('RANDQEICP_%.2f_%d',d,n);
        fprintf(' * the problem %s is generated\n',filename);
        save(sprintf('RANDQEICP\\%s.mat',filename),'n','A','B','C','D','G','mu','z0');
    end
end

%% ill-conditioned dense
randbnd = [-1,1];
%randbnd = [-10,10];
%randbnd = [0,100];

for n=[50,100,200,400,600]
    A=eye(n);
    B=rand(n)*(randbnd(2)-randbnd(1))+randbnd(1);
    B=(B'+B)/2;
    F=rand(n);
    C=-eye(n)-F'*F;
    D = [A,zeros(n);zeros(n),-C];
    G = [-B,-C;-C,zeros(n)];
    varmu=sdpvar(1);
    optimize([G+varmu*D>=0],varmu,sdpsettings('solver','mosek','verbose',0));
    %mineigsG=min(eig(G));
    %mineigsD=min(eig(D));
    %mu = abs(mineigsG)/mineigsD + 1;
    mu = value(varmu)+1;
    G = G + mu*D;
    x0 = rand(n,1);
    lambda0 = (-x0'*B*x0 + sqrt((x0'*B*x0)^2 - 4*(x0'*A*x0)*(x0'*C*x0)))/(2*x0'*A*x0);
    z0 = [lambda0*x0;x0]/(1+lambda0);
    filename = sprintf('RANDQEICP_%d_%d_%d',randbnd(1),randbnd(2),n);
    fprintf(' * the problem %s is generated\n',filename);
    save(['RANDQEICP/',filename],'n','A','B','C','D','G','mu','z0');
end