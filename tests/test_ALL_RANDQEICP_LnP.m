%clear;
range_n=[50,100,200,400,600]; 
range_d=[0.05,0.1,0.5,0.7,0.9];

%% initialize some parameters for tests
% flag setting
testdca=1;
testbdca=1;
testknitro=1;
testfmincon=1;
testfiltersd=1;
mkdir('RESULT');

%% test DCA, BDCA, fmincon, knitro and filtersd for (LnP) on RANDEICP dataset
for d=range_d
    for n=range_n
        % read data
        fname=sprintf('RANDQEICP_%.2f_%d',d,n);
        fprintf('Loading model from %s.mat\n',fname);
        load(['RANDQEICP//',fname,'.mat']);
        
        %% Test DCA for (LnP)
        if testdca
            fprintf('Solving model %s using DCA\n',fname);

            % create a dc function object
            dcf=dcfunc;
            dcf.f=@(x,opt)fobj_eval_LnP(x,G,D,opt);
            
            % create a dc problem object
            mydcp=dcp(dcf,D);
            
            % create a dca object            
            mydca = dca(mydcp,z0);
            mydca.A=G;
            mydca.B=D;
            condA=condest(G);
            condB=condest(D);
            mydca.etabar=10*max([condA,condB]);
            mydca.Lg=mydca.etabar;
            mydca.plot=0;
            mydca.tolf=0;
            mydca.tolx=1e-8;
            mydca.fistatol=1e-6;
            mydca.maxiter=1e+4;
            mydca.verbose = 0;
            mydca.linesearch=0;
            mydca.model='LnP';
            
            % solve model using dca
            status=mydca.optimize();
            
            % get results
            cputime=status.time;
            iters = status.iter;
            % solution of EICP
            xopt = mydca.xopt;
            fopt = mydca.fopt;
            lambda = (xopt'*mydca.A*xopt)/(xopt'*mydca.B*xopt);
            % from a solution of EICP to QEICP
            lambda=lambda-mu;
            xopt=(1+lambda)*xopt(n+1:2*n);
            w = lambda^2*A*xopt + lambda*B*xopt + C*xopt;
            err = [norm(min(w,0)),norm(min(xopt,0)),abs(xopt'*w)];
            fprintf('Solution for DCA (time %.3f sec, obj %.5f iters %d, lambda %.5f, err %.3e)\n',cputime,fopt,iters,lambda,sum(err));
            save(sprintf('RESULT//DCA_LnP_%s.mat',fname),'xopt','fopt','cputime','iters','lambda','w','err','mu');
            
        end
        
        %% Test BDCA for (LnP)
        if testbdca
            fprintf('Solving model %s using BDCA\n',fname);
            
            % create a dc function object
            dcf=dcfunc;
            dcf.f=@(x,opt)fobj_eval_LnP(x,G,D,opt);
            
            % create a dc problem object
            mydcp=dcp(dcf,D);
            
            % create a dca object            
            mydca = dca(mydcp,z0);
            mydca.A=G;
            mydca.B=D;
            condA=condest(G);
            condB=condest(D);
            mydca.etabar=10*max([condA,condB]);
            mydca.Lg=mydca.etabar;
            mydca.plot=0;
            mydca.tolf=0;
            mydca.tolx=1e-8;
            mydca.fistatol=1e-6;
            mydca.maxiter=1e+4;
            mydca.verbose = 0;
            mydca.linesearch=1;
            mydca.linesearch_type='exact';
            mydca.model='LnP';
            
            % solve model using dca
            status=mydca.optimize();
            
            % get results
            cputime=status.time;
            iters = status.iter;
            % solution of EICP
            xopt = mydca.xopt;
            fopt = mydca.fopt;
            lambda = (xopt'*mydca.A*xopt)/(xopt'*mydca.B*xopt);
            % from a solution of EICP to QEICP
            lambda=lambda-mu;
            xopt=(1+lambda)*xopt(n+1:2*n);
            w = lambda^2*A*xopt + lambda*B*xopt + C*xopt;
            err = [norm(min(w,0)),norm(min(xopt,0)),abs(xopt'*w)];
            fprintf('Solution for BDCA (time %.3f sec, obj %.5f iters %d, lambda %.5f, err %.3e)\n',cputime,fopt,iters,lambda,sum(err));
            save(sprintf('RESULT//BDCA_LnP_%s.mat',fname),'xopt','fopt','cputime','iters','lambda','w','err','mu');
            
        end
        
        %% Test knitro for (LnP)
        if testknitro
            fprintf('Solving model %s using KNITRO\n',fname);
            fobj=@(x)log((x'*D*x)/(x'*G*x));
            Aeq=ones(1,2*n);
            beq=1;
            lb=zeros(2*n,1);
            ub=[];
            tic
            [xopt,fopt,exitflag,info]=knitromatlab(fobj,z0,[],[],Aeq,beq,lb,ub,[],[],optimset('Display','off'));
            
            % get results
            cputime=toc;
            iters = info.iterations;
            lambda = (xopt'*G*xopt)/(xopt'*D*xopt);
            lambda=lambda-mu;
            xopt=(1+lambda)*xopt(n+1:2*n);
            w = lambda^2*A*xopt + lambda*B*xopt + C*xopt;
            err = [norm(min(w,0)),norm(min(xopt,0)),abs(xopt'*w)];
            fprintf('Solution for KNITRO (time %.3f sec, obj %.5f iters %d, lambda %.5f, err %.3e)\n',cputime,fopt,iters,lambda,sum(err));
            save(sprintf('RESULT//KNITRO_LnP_%s.mat',fname),'xopt','fopt','cputime','iters','lambda','w','err','mu');
        end
        
        %% Test FILTERSD for (LnP)
        if testfiltersd
            fprintf('Solving model %s using FLISTERSD\n',fname);
            fobj=@(x)log((x'*D*x)/(x'*G*x));
            gradfobj=@(x)(2*D*x)/(x'*G*x) - (2*G*x*(x'*D*x))/(x'*G*x)^2;
            lcon=@(x)ones(1,2*n)*x;
            ljac=@(x)ones(1,2*n);
            Aeq=ones(1,2*n);
            beq=1;
            lb=zeros(2*n,1);
            ub=[];
            tic
            
            [xopt,fopt,exitflag,info] = filtersd(fobj,gradfobj,z0,lb,ub,lcon,ljac,1,1);
            
            % get results
            cputime=toc;
            iters = info.niter;
            lambda = (xopt'*G*xopt)/(xopt'*D*xopt);
            lambda=lambda-mu;
            xopt=(1+lambda)*xopt(n+1:2*n);
            w = lambda^2*A*xopt + lambda*B*xopt + C*xopt;
            err = [norm(min(w,0)),norm(min(xopt,0)),abs(xopt'*w)];
            fprintf('Solution for FILTERSD (time %.3f sec, obj %.5f iters %d, lambda %.5f, err %.3e)\n',cputime,fopt,iters,lambda,sum(err));
            save(sprintf('RESULT//FILTERSD_LnP_%s.mat',fname),'xopt','fopt','cputime','iters','lambda','w','err','mu');
        end
        
        %% Test FMINCON for (LnP)
        if testfmincon
            fprintf('Solving model %s using FMINCON\n',fname);
            fobj=@(x)log((x'*D*x)/(x'*G*x));
            Aeq=ones(1,2*n);
            beq=1;
            lb=zeros(2*n,1);
            ub=[];
            tic
            [xopt,fopt,exitflag,info]=fmincon(fobj,z0,[],[],Aeq,beq,lb,ub,[],optimset('Display','off','MaxFunEvals',1e+6));
            
            % get results
            cputime=toc;
            iters = info.iterations;
            lambda = (xopt'*G*xopt)/(xopt'*D*xopt);
            lambda=lambda-mu;
            xopt=(1+lambda)*xopt(n+1:2*n);
            w = lambda^2*A*xopt + lambda*B*xopt + C*xopt;
            err = [norm(min(w,0)),norm(min(xopt,0)),abs(xopt'*w)];
            fprintf('Solution for FMINCON (time %.3f sec, obj %.5f iters %d, lambda %.5f, err %.3e)\n',cputime,fopt,iters,lambda,sum(err));
            save(sprintf('RESULT//FMINCON_LnP_%s.mat',fname),'xopt','fopt','cputime','iters','lambda','w','err','mu');
        end
    end
end

fprintf('All test finished!\n');