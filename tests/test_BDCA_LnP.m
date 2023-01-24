% require n,A,B,mu,x0
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
%x0 = x0/sqrt(x0'*B*x0);

mydca = dca(mydcp,x0);
mydca.A=A;
mydca.B=B;
mydca.etabar=n;
%condA=cond(A);
%condB=cond(B);
%mydca.etabar=9*max([condA,condB]);
mydca.Lg=mydca.etabar;
mydca.plot=0;
mydca.tolf=0;
mydca.tolx=1e-6;
mydca.fistatol=1e-6;
mydca.maxiter=1e+4;
mydca.verbose = 0;
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

% verification for EiCP

lamb = (xopt'*A*xopt)/(xopt'*B*xopt);
w=lamb*B*xopt - A*xopt;
fprintf(' * w positivity error: %e\n',norm(min(w,0)));
fprintf(' * x positivity error: %e\n',norm(min(xopt,0)));
fprintf(' * <x,w> : %e\n',xopt'*w);
fprintf(' * lambda : %f\n',lamb-mu);