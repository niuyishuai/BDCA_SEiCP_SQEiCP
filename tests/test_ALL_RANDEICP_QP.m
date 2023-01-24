%clear;
range_n=[50,100,200,400,600,800];
bndlst={[-1,1],[-10,10]};

%% initialize some parameters for tests
% flag setting
testdca=1;
testbdca=1;
testknitro=1;
testfmincon=1;
testfiltersd=1;
mkdir('RESULT');

%% test DCA, BDCA, fmincon, knitro and filtersd for (QP) on RANDEICP dataset
for bnd=bndlst
    for n=range_n
        % read data
        fname=sprintf('RANDEICP_%d_%d_%d',bnd{1}(1),bnd{1}(2),n);
        fprintf('Loading model from %s.mat\n',fname);
        load(['RANDEICP//',fname,'.mat']);
        
        %% Test DCA for (QP)
        if testdca
            fprintf('Solving model %s using DCA\n',fname);
            
            % create a dc function object
            dcf=dcfunc;
            dcf.f=@(x,opt)fobj_eval_QP(x,A,opt);
            
            % create a dc problem object
            mydcp=dcp(dcf,B);
            
            % create a dca object
            %x0 = rand(n,1);
            %x0 = x0/sqrt(x0'*BB*x0);
            
            mydca = dca(mydcp,x0);
            mydca.A=A;
            mydca.B=B;
            mydca.plot=0;
            mydca.tolf=0;
            mydca.tolx=1e-5;
            mydca.verbose = 0;
            mydca.linesearch=0;
            mydca.model='QP';
            
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
            save(sprintf('RESULT//DCA_QP_%s',fname),'xopt','fopt','cputime','iters','lambda','w','err','mu');
            
        end
        
        %% Test BDCA for (QP)
        if testbdca
            fprintf('Solving model %s using BDCA\n',fname);
            
            % create a dc function object
            dcf=dcfunc;
            dcf.f=@(x,opt)fobj_eval_QP(x,A,opt);
            
            % create a dc problem object
            mydcp=dcp(dcf,B);
            
            % create a dca object
            %x0 = rand(n,1);
            %x0 = x0/sqrt(x0'*BB*x0);
            
            mydca = dca(mydcp,x0);
            mydca.A=A;
            mydca.B=B;
            mydca.plot=0;
            mydca.tolf=0;
            mydca.tolx=1e-5;
            mydca.verbose = 0;
            mydca.linesearch=1;
            mydca.linesearch_type='exact';
            mydca.model='QP';
            
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
            save(sprintf('RESULT//BDCA_QP_%s',fname),'xopt','fopt','cputime','iters','lambda','w','err','mu');
            
        end
        
        %% Test knitro for (QP)
        if testknitro
            fprintf('Solving model %s using KNITRO\n',fname);
            fobj=@(x)-x'*A*x;
            nlcon=@(x)deal(x'*B*x-1,0);    
            lb=zeros(n,1);
            ub=[];
            tic
            [xopt,fopt,exitflag,info]=knitromatlab(fobj,x0,[],[],[],[],lb,ub,nlcon,[],optimset('Display','off'));
            
            % get results
            cputime=toc;
            iters = info.iterations;
            lambda = (xopt'*A*xopt)/(xopt'*B*xopt);
            w=lambda*B*xopt - A*xopt;
            err = [norm(min(w,0)),norm(min(xopt,0)),xopt'*w];
            fprintf('Solution for KNITRO (time %.3f sec, obj %.5f iters %d, lambda %.5f, err %.3e)\n',cputime,fopt,iters,lambda,sum(err));
            save(sprintf('RESULT//KNITRO_QP_%s',fname),'xopt','fopt','cputime','iters','lambda','w','err','mu');
        end
        
        %% Test FILTERSD for (QP)
        if testfiltersd
            fprintf('Solving model %s using FLISTERSD\n',fname);
            fobj=@(x)-x'*A*x;
            gradfobj=@(x)-2*A*x;
            nlcon=@(x)x'*B*x-1;
            nljac=@(x)2*(B*x)';
            lb=zeros(n,1);
            ub=[];
            tic
            [xopt,fopt,exitflag,info]=filtersd(fobj,gradfobj,x0,lb,ub,nlcon,nljac,-inf,0);
            
            % get results
            cputime=toc;
            iters = info.niter;
            lambda = (xopt'*A*xopt)/(xopt'*B*xopt);
            w=lambda*B*xopt - A*xopt;
            err = [norm(min(w,0)),norm(min(xopt,0)),xopt'*w];
            fprintf('Solution for FILTERSD (time %.3f sec, obj %.5f iters %d, lambda %.5f, err %.3e)\n',cputime,fopt,iters,lambda,sum(err));
            save(sprintf('RESULT//FILTERSD_QP_%s',fname),'xopt','fopt','cputime','iters','lambda','w','err','mu');
        end
        
        %% Test FMINCON for (QP)
        if testfmincon
            fprintf('Solving model %s using FMINCON\n',fname);
            fobj=@(x)-x'*A*x;
            nlcon=@(x)deal(x'*B*x-1,0);    
            lb=zeros(n,1);
            ub=[];
            tic
            [xopt,fopt,exitflag,info]=fmincon(fobj,x0,[],[],[],[],lb,ub,nlcon,optimset('Display','off','MaxFunEvals',1e+6));
            
            % get results
            cputime=toc;
            iters = info.iterations;
            lambda = (xopt'*A*xopt)/(xopt'*B*xopt);
            w=lambda*B*xopt - A*xopt;
            err = [norm(min(w,0)),norm(min(xopt,0)),xopt'*w];
            fprintf('Solution for FMINCON (time %.3f sec, obj %.5f iters %d, lambda %.5f, err %.3e)\n',cputime,fopt,iters,lambda,sum(err));
            save(sprintf('RESULT//FMINCON_QP_%s',fname),'xopt','fopt','cputime','iters','lambda','w','err','mu');
        end
    end
end

fprintf('All test finished!\n');