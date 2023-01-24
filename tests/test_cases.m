%% test random data sample
n=200;
randbnd=[-10,10];
A=rand(n)*(randbnd(2)-randbnd(1))+randbnd(1);
A=(A'+A)/2;
B=eye(n);
varmu=sdpvar(1);
optimize([A+varmu*B>=0],varmu,sdpsettings('solver','mosek','verbose',0));
mu=value(varmu)+1;
A = A + mu*B;
x0 = rand(n,1);
%x0 = x0/sum(x0);

%% test NEP data sample
[A,n]=mmread('NEP\\dwa512.mtx');
A=(A+A')/2;
B=eye(n);
varmu=sdpvar(1);
optimize([A+varmu*B>=0],varmu,sdpsettings('solver','mosek','verbose',0));
mu=value(varmu)+1;
A = A + mu*B;
x0 = rand(n,1);

%% Test BDCA for (QP)
% create a dc function object
dcf=dcfunc;
dcf.f=@(x,opt)fobj_eval_QP(x,A,opt);

% create a dc problem object
mydcp=dcp(dcf,[]);
fprintf('End of initialization.\n');

% Test BDCA solver
fprintf('Solving EiCP using BDCA.\n');

% create a dca object
mydca = dca(mydcp,x0);
mydca.A=A;
mydca.B=B;
mydca.plot=0;
mydca.tolf=0;
mydca.tolx=1e-6;
mydca.maxiter=1e+4;
mydca.verbose = 1;
mydca.linesearch=1;
mydca.linesearch_type='exact';
mydca.model='QP';

% solve model using BDCA
status=mydca.optimize();

% get results
bdcatime=status.time;
bdca_times = bdcatime;
bdca_objs = mydca.fopt;
bdca_iters = status.iter;
xopt=mydca.xopt;

fprintf('Solution for BDCA: time %.3f sec, obj %.6f iters %d\n',bdca_times,bdca_objs,status.iter);

%% Test BDCA for (LnP)
% create a dc function object
dcf=dcfunc;
dcf.f=@(x,opt)fobj_eval_LnP(x,A,B,opt);

% create a dc problem object
mydcp=dcp(dcf,[]);
fprintf('End of initialization.\n');

% Test BDCA solver
fprintf('Solving EiCP using BDCA.\n');
% create a dca object
%x0 = rand(n,1);
%x0 = x0/sum(x0);
%x0 = x0/sqrt(x0'*BB*x0);

mydca = dca(mydcp,x0);
mydca.A=A;
mydca.B=B;
condA=cond(A);
condB=cond(B);
mydca.etabar=n;
mydca.Lg=mydca.etabar; 
mydca.plot=1;
mydca.tolf=0;
mydca.tolx=1e-8;
mydca.fistatol=1e-6;
mydca.maxiter=1e+4;
mydca.verbose = 1;
mydca.linesearch=1;
mydca.linesearch_type='exact';
mydca.model='LnP';

% solve model using BDCA
status=mydca.optimize();

% get results
bdcatime=status.time;
bdca_times = bdcatime;
bdca_objs = mydca.fopt;
bdca_iters = status.iter;
xopt=mydca.xopt;

fprintf('Solution for BDCA: time %.3f sec, obj %.4e iters %d\n',bdca_times,bdca_objs,status.iter);


%% other solvers for QP
if 0
    % with yalmip
    x=sdpvar(n,1);
    assign(x,x0);
    optimize([x'*B*x<=1, x>=0],-x'*A*x,sdpsettings('solver','fmincon','usex0',1,'verbose',0))
    xopt=value(x);
    value(-x'*A*x)
else
    % without yalmip
    fobj=@(x)-x'*A*x;
    gradfobj=@(x)-2*A*x;
    nlcon=@(x)x'*B*x-1;
    nljac=@(x)2*(B*x)';
    lb=zeros(n,1);
    ub=[];
    fmincon_nlcon=@(x)deal(nlcon(x),0);    
        
    tic
    %[xopt,fval,exitflag,info]=fmincon(fobj,x0,[],[],[],[],lb,ub,fmincon_nlcon,optimset('MaxFunEvals',1e+5))
    [xopt,fval,exitflag,info]=knitromatlab(fobj,x0,[],[],[],[],lb,ub,fmincon_nlcon)
    %[xopt,fval,exitflag,info] = filtersd(fobj,gradfobj,x0,lb,ub,nlcon,nljac,-inf,0)
    toc
end

%% other solvers for LnP
if 0
    % with yalmip
    x=sdpvar(n,1);
    assign(x,x0);
    optimize([sum(x)==1, x>=0],log((x'*B*x)/(x'*A*x)),sdpsettings('solver','fmincon','usex0',1,'verbose',0));
    xopt=value(x);
    value(log((x'*B*x)/(x'*A*x)))
else
    % without yalmip
    fobj=@(x)log((x'*B*x)/(x'*A*x));
    gradfobj=@(x)(2*B*x)/(x'*A*x) - (2*A*x*(x'*B*x))/(x'*A*x)^2;
    lcon=@(x)ones(1,n)*x;
    ljac=@(x)ones(1,n);
    Aeq=ones(1,n);
    beq=1;
    lb=zeros(n,1);
    ub=[];
    tic
    %[xopt,fval,exitflag,output]=fmincon(fobj,x0,[],[],Aeq,beq,lb,ub,[],optimset('MaxFunEvals',1e+5))
    [xopt,fval,exitflag,info]=knitromatlab(fobj,x0,[],[],Aeq,beq,lb,ub)
    %[xopt,fval,exitflag,info] = filtersd(fobj,gradfobj,x0,lb,ub,lcon,ljac,1,1)
    toc
end

%% verification for EiCP

lamb = (xopt'*A*xopt)/(xopt'*B*xopt);
w=lamb*B*xopt - A*xopt;
fprintf(' * w positivity error: %e\n',norm(min(w,0)));
fprintf(' * x positivity error: %e\n',norm(min(xopt,0)));
fprintf(' * <x,w> : %e\n',xopt'*w);
fprintf(' * lambda : %f\n',lamb);
