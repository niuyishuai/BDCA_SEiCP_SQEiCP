%% initialize some parameters for tests
% flag setting
testdca=1;
testbdca=1;
testknitro=1;
testfmincon=1;
testfiltersd=1;
mkdir('RESULT');

%% test DCA, BDCA, fmincon, knitro and filtersd for (LnP) on NEP dataset
NEPlst=dir('NEP\\*.mat');
for i=1:length(NEPlst)
    % read data
    fname=NEPlst(i).name;
    fprintf('Loading model from %s\n',fname);
    load(sprintf('NEP\\%s',fname));
            
    %% Test DCA for (LnP)
    if testdca
        fprintf('Solving model %s using DCA\n',fname);
        
        % create a dc function object
        dcf=dcfunc;
        dcf.f=@(x,opt)fobj_eval_LnP(x,A,B,opt);
        
        % create a dc problem object
        mydcp=dcp(dcf,B);
        
        % create a dca object
        %x0 = rand(n,1);
        %x0 = x0/sqrt(x0'*BB*x0);
        
        mydca = dca(mydcp,x0);
        mydca.A=A;
        mydca.B=B;
        condA=cond(A);
        condB=cond(B);
        mydca.etabar=n;
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
        xopt = mydca.xopt;
        fopt = mydca.fopt;
        lambda = (xopt'*A*xopt)/(xopt'*B*xopt);
        w=lambda*B*xopt - A*xopt;
        err = [norm(min(w,0)),norm(min(xopt,0)),xopt'*w];
        fprintf('Solution for DCA (time %.3f sec, obj %.5f iters %d, lambda %.5f, err %.3e)\n',cputime,fopt,iters,lambda,sum(err));
        save(sprintf('RESULT//DCA_LnP_%s',fname),'xopt','fopt','cputime','iters','lambda','w','err','mu');
        
    end
    
    %% Test BDCA for (LnP)
    if testbdca
        fprintf('Solving model %s using BDCA\n',fname);
        
        % create a dc function object
        dcf=dcfunc;
        dcf.f=@(x,opt)fobj_eval_LnP(x,A,B,opt);
        
        % create a dc problem object
        mydcp=dcp(dcf,B);
        
        % create a dca object
        %x0 = rand(n,1);
        %x0 = x0/sqrt(x0'*BB*x0);
        
        mydca = dca(mydcp,x0);
        mydca.A=A;
        mydca.B=B;
        condA=cond(A);
        condB=cond(B);
        mydca.etabar=n;
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
        xopt = mydca.xopt;
        fopt = mydca.fopt;
        lambda = (xopt'*A*xopt)/(xopt'*B*xopt);
        w=lambda*B*xopt - A*xopt;
        err = [norm(min(w,0)),norm(min(xopt,0)),xopt'*w];
        fprintf('Solution for BDCA (time %.3f sec, obj %.5f iters %d, lambda %.5f, err %.3e)\n',cputime,fopt,iters,lambda,sum(err));
        save(sprintf('RESULT//BDCA_LnP_%s',fname),'xopt','fopt','cputime','iters','lambda','w','err','mu');
        
    end
    
    %% Test knitro for (LnP)
    if testknitro
        fprintf('Solving model %s using KNITRO\n',fname);
        fobj=@(x)log((x'*B*x)/(x'*A*x));
        Aeq=ones(1,n);
        beq=1;
        lb=zeros(n,1);
        ub=[];
        tic
        [xopt,fopt,exitflag,info]=knitromatlab(fobj,x0,[],[],Aeq,beq,lb,ub,[],[],optimset('Display','off'));
        
        % get results
        cputime=toc;
        iters = info.iterations;
        lambda = (xopt'*A*xopt)/(xopt'*B*xopt);
        w=lambda*B*xopt - A*xopt;
        err = [norm(min(w,0)),norm(min(xopt,0)),xopt'*w];
        fprintf('Solution for KNITRO (time %.3f sec, obj %.5f iters %d, lambda %.5f, err %.3e)\n',cputime,fopt,iters,lambda,sum(err));
        save(sprintf('RESULT//KNITRO_LnP_%s',fname),'xopt','fopt','cputime','iters','lambda','w','err','mu');
    end
    
    %% Test FILTERSD for (LnP)
    if testfiltersd
        fprintf('Solving model %s using FLISTERSD\n',fname);
        fobj=@(x)log((x'*B*x)/(x'*A*x));
        gradfobj=@(x)(2*B*x)/(x'*A*x) - (2*A*x*(x'*B*x))/(x'*A*x)^2;
        lcon=@(x)ones(1,n)*x;
        ljac=@(x)ones(1,n);
        lb=zeros(n,1);
        ub=[];
        tic
        
        [xopt,fopt,exitflag,info] = filtersd(fobj,gradfobj,x0,lb,ub,lcon,ljac,1,1);
        
        % get results
        cputime=toc;
        iters = info.niter;
        lambda = (xopt'*A*xopt)/(xopt'*B*xopt);
        w=lambda*B*xopt - A*xopt;
        err = [norm(min(w,0)),norm(min(xopt,0)),xopt'*w];
        fprintf('Solution for FILTERSD (time %.3f sec, obj %.5f iters %d, lambda %.5f, err %.3e)\n',cputime,fopt,iters,lambda,sum(err));
        save(sprintf('RESULT//FILTERSD_LnP_%s',fname),'xopt','fopt','cputime','iters','lambda','w','err','mu');
    end
    
    %% Test FMINCON for (LnP)
    if testfmincon
        fprintf('Solving model %s using FMINCON\n',fname);
        fobj=@(x)log((x'*B*x)/(x'*A*x));
        Aeq=ones(1,n);
        beq=1;
        lb=zeros(n,1);
        ub=[];
        tic
        [xopt,fopt,exitflag,info]=fmincon(fobj,x0,[],[],Aeq,beq,lb,ub,[],optimset('Display','off','MaxFunEvals',1e+6));
        
        % get results
        cputime=toc;
        iters = info.iterations;
        lambda = (xopt'*A*xopt)/(xopt'*B*xopt);
        w=lambda*B*xopt - A*xopt;
        err = [norm(min(w,0)),norm(min(xopt,0)),xopt'*w];
        fprintf('Solution for FMINCON (time %.3f sec, obj %.5f iters %d, lambda %.5f, err %.3e)\n',cputime,fopt,iters,lambda,sum(err));
        save(sprintf('RESULT//FMINCON_LnP_%s',fname),'xopt','fopt','cputime','iters','lambda','w','err','mu');
    end
end

fprintf('All test finished!\n');