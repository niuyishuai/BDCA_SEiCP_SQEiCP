% require A,B,C,D,G,mu,n,z0
%% Test BDCA for SEiCP(G,D)
% create a dc function object
dcf=dcfunc;
dcf.f=@(x,opt)fobj_eval_LnP(x,G,D,opt);

% create a dc problem object
mydcp=dcp(dcf,[]);

% create a dca object
mydca = dca(mydcp,z0);
mydca.A=G;
mydca.B=D;
condA=condest(G);
condB=condest(D);
mydca.etabar=9*max([condA,condB]);
mydca.Lg=mydca.etabar;
mydca.plot=0;
mydca.tolf=0;
mydca.tolx=1e-8;
mydca.fistatol=1e-6;
mydca.maxiter=1e+4;
mydca.verbose = 1;
mydca.linesearch=1;
mydca.model='LnP';

% solve model using BDCA
status=mydca.optimize();

% get results
cputime=status.time;
iters = status.iter;
% results for EiCP
xopt = mydca.xopt;
fopt = mydca.fopt;
lambda = (xopt'*mydca.A*xopt)/(xopt'*mydca.B*xopt);
w=lambda*mydca.B*xopt - mydca.A*xopt;
err = [norm(min(w,0)),norm(min(xopt,0)),abs(xopt'*w)];
fprintf('Solution for EiCP formulation by BDCA (time %.3f sec, obj %.5f iters %d, lambda %.5f, err %.3e)\n',cputime,fopt,iters,lambda,sum(err));
% from EICP to QEICP 
%Gorg = G - mu*D;
lambda_qeicp=lambda-mu;
xopt_qeicp=(1+lambda_qeicp)*xopt(n+1:2*n); 
w_qeicp = lambda_qeicp^2*A*xopt_qeicp + lambda_qeicp*B*xopt_qeicp + C*xopt_qeicp;
err = [norm(min(w_qeicp,0)),norm(min(xopt_qeicp,0)),xopt_qeicp'*w_qeicp];
fprintf('Solution for QEiCP by DCA (time %.3f sec, obj %.5f iters %d, lambda %.5f, err %.3e)\n',cputime,fopt,iters,lambda_qeicp,sum(err));            

