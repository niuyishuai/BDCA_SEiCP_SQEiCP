classdef dca < matlab.mixin.Copyable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % general dc programming solver (dca class)
    %
    % Author: yi-shuai niu
    % 2019-4: Initial coding
    % 2023-3: Modified for solving SEiCP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        dcp % dc program
        x0 % initial point
    end
    properties(GetAccess = public,SetAccess = private)  % read only
        fopt = inf % objective value
        xopt = [] % optimal solution
        iter = 0 % iterations of dca
    end
    properties
        tolf = 0 % tolerence for objective function
        tolx = 1e-6 % tolerence for iterative point
        maxiter = 1e+4 % max iterations for dca
        plot = 0  %1: draw iterations of dca, 0: otherwise
        verbose = 1  %1: display iterations, 0: otherwise
        linesearch = true % use line search for acceleration
        linesearch_type = 'exact' %'armijo|exact'
        plotlinetype='b-s'; % if plot = 1, this option set the plot line type
        plotinnewfig=true; % if true, we flot in a new fig, otherwise, we plot in fig 1
        % add some personal attributes here
        A;
        B;
        model; % QP|LnP
        % parameters for (LnP) model only
        etabar;
        Lg;
        fistatol=1e-6;
    end
    
    methods
        function obj = dca(dcp,x0)
            % dca constructor
            % obj = dca(dcp,x0)
            % where dcp is a dc program object
            % x0 is an initial point. It will be a random point if x0 is not given.
            if nargin==1 % if x0 is not given
                obj.dcp = dcp;
                x0 = rand(size(dcp.X));
                obj.x0 = x0(:);
            elseif nargin==2 % if x0 is given
                obj.dcp = dcp;
                obj.x0 = x0(:);
            else
                error('wrong input arguments.');
            end
        end
        function xopt = get.xopt(obj)
            % get optimal solution
            xopt = obj.xopt;
        end
        function fopt = get.fopt(obj)
            % get objective value
            fopt = obj.fopt;
        end
        function set.tolf(obj,val)
            % set tolerence of objective function
            obj.tolf = val;
        end
        function set.tolx(obj,val)
            % set tolerence of iterative point
            obj.tolx = val;
        end
        function set.maxiter(obj,val)
            % set max iterations of dca
            obj.maxiter = val;
        end
        function set.plot(obj,val)
            % set plotting option of dca
            obj.plot = val;
        end
        function set.verbose(obj,val)
            % set verbose option of dca
            obj.verbose = val;
        end
        function set.linesearch(obj,yn)
            % set line search
            obj.linesearch = yn;
        end
        % dca algorithm for solving dcp with starting point x0
        function status = optimize(obj)
            % dca optimizer
            % status.flag : 0 dca converges with tolx or tolf
            %               1 maxiter exceed
            %               2 problem infeasible or unbounded
            % status.info : solution informations.
            % status.iter : number of iterations of dca.
            % status.time : cpu time for dca (sec.)
            % status.avgt : average time for each iteration (sec.)
            
            % initialization
            obj.iter = 0;
            xk = obj.x0;
            n = numel(xk);
            vecz = zeros(n,1);
            veczt = vecz';
            fk=obj.dcp.F.f(xk,1); % 0 for computing both f(xk) and df(xk), 1 for computing f(xk) only
            mosekcmd='minimize info echo(0) statuskeys(1) symbcon';
            % plotting if active
            if (obj.plot==1)
                if (obj.plotinnewfig)
                    figure
                else
                    figure(1);
                end
            end
            % display if active
            if (obj.verbose == 1)
                fprintf('------------------------------------------------------------\n');
                if obj.linesearch
                    fprintf('BDCA for solving SEiCP (version 1.0) \n');
                else
                    fprintf('DCA for solving SEiCP (version 1.0) \n');
                end
                if obj.linesearch
                    switch obj.linesearch_type
                        case 'exact'
                            fprintf('* activate exact linesearch acceleration\n');
                    end
                end
                fprintf('------------------------------------------------------------\n');
                fprintf('Iterations | Objective values |   Delta x   |   Delta f \n');
            end
            cputime=tic;
            % BDCA/DCA for solving QP or LnP models
            switch obj.model
                case 'QP'
                    while obj.iter < obj.maxiter
                        obj.iter = obj.iter+1;
                        yk = -obj.A*xk; % compte dh(xk), a factor 2 is ignored
                        % initialize mosek problem
                        subprob.c=yk;
                        subprob.a=veczt;
                        subprob.blc=-inf;
                        subprob.buc=1;
                        subprob.blx=vecz;
                        subprob.bux=[];
                        subprob.sol.itr.xx=xk;
                        [subprob.qcsubi,subprob.qcsubj,subprob.qcval] = find(tril(obj.B));
                        subprob.qcval = 2*subprob.qcval; %MOSEK includes 1/2x'Qx
                        subprob.qcsubk=ones(size(subprob.qcval));
                        
                        % solve convex subproblem via mosek
                        [~,res] = mosekopt(mosekcmd,subprob);
                        xk1 = res.sol.itr.xx;
                        %xk1 = mosekqcqp([],yk,[],[],[],obj.B,vecz,-inf,1,vecz,[],xk);
                        
                        [fk1,dfk1] = obj.dcp.F.f(xk1,0);
                        % accelerate with line search
                        if obj.linesearch == true
                            activesetxk = checkactiveset(xk);
                            activesetxk1 = checkactiveset(xk1);
                            dk = xk1 - xk;
                            
                            if (min(activesetxk - activesetxk1)>=0 && dk'*dfk1 < 0 )
                                switch obj.linesearch_type
                                    case 'exact'
                                        Ik = dk<0;
                                        if sum(Ik)==0 % Ik is empty
                                            alphak = (-xk1'*obj.B*dk + sqrt((xk1'*obj.B*dk)^2 - dk'*obj.B*dk*(xk1'*obj.B*xk1-1)))/(dk'*obj.B*dk);
                                        else
                                            alphak = min([min(-xk1(Ik)./dk(Ik)),(-xk1'*obj.B*dk + sqrt((xk1'*obj.B*dk)^2 - dk'*obj.B*dk*(xk1'*obj.B*xk1-1)))/(dk'*obj.B*dk)]);
                                        end
                                        if alphak > 1e-8 % if alphak is not too small, update xk1
                                            xacc = xk1 + alphak*dk;
                                            facc = obj.dcp.F.f(xacc,1);
                                            if facc < fk1
                                                if obj.verbose == 1
                                                    fprintf('accelerated: obj-reduced %17.3e  x-moved %17.3e \n',facc-fk1,norm(xacc-xk1));
                                                end
                                                xk1 = xacc;
                                                fk1 = facc;
                                            end
                                        end
                                end
                                
                            end
                        end
                        
                        % compute errors
                        normx = norm(xk1-xk);
                        normf = abs(fk1-fk);
                        % display iterations of dca if active
                        if (obj.verbose == 1)
                            fprintf('%5d %19.5e %17.3e %13.3e\n',obj.iter,fk1,normx,normf);
                        end
                        % plotting if active
                        if (obj.plot==1)
                            myplotf(fk,fk1,obj.iter,obj.plotlinetype);
                        end
                        % check stopping
                        if (normx < obj.tolx*(1+norm(xk1)) || normf < obj.tolf*(1+abs(fk1)))
                            if (obj.verbose == 1)
                                fprintf('------------------------------------------------------------\n');
                            end
                            obj.fopt = fk1;
                            obj.xopt = xk1;
                            status = setstatus(toc(cputime),obj.iter,0,'Successfully solved.');
                            return;
                        end
                        
                        xk = xk1;
                        fk = fk1;
                    end
                case 'LnP'
                    while obj.iter < obj.maxiter
                        obj.iter = obj.iter+1;
                        
                        xk1 = FISTA(xk,1,obj.fistatol,obj.etabar,obj.A,obj.B,obj.Lg);
                        
                        % accelerate with line search
                        if obj.linesearch == true
                            [fk1,dfk1] = obj.dcp.F.f(xk1,0);
                            activesetxk = checkactiveset(xk);
                            activesetxk1 = checkactiveset(xk1);
                            dk = xk1 - xk;
                            
                            if (min(activesetxk - activesetxk1)>=0 && dk'*dfk1 < 0)
                                switch obj.linesearch_type
                                    case 'exact'
                                        Ik = dk<0;
                                        Bdk=obj.B*dk;
                                        Adk=obj.A*dk;
                                        ploycoefs = zeros(1,6);
                                        ploycoefs(1)=dk'*Bdk; %a1
                                        ploycoefs(2)=2*xk1'*Bdk; %b1
                                        ploycoefs(3)=xk1'*obj.B*xk1; %c1
                                        ploycoefs(4)=dk'*Adk; %a2
                                        ploycoefs(5)=2*xk1'*Adk; %b2
                                        ploycoefs(6)=xk1'*obj.A*xk1; %c2
                                        if sum(Ik)==0 % Ik is empty
                                            alphak = findbestsz(ploycoefs,inf);
                                        else
                                            alphak = findbestsz(ploycoefs,min(-xk1(Ik)./dk(Ik)));
                                        end
                                        if alphak > 1e-8 % if alphak is not too small, update xk1
                                            xacc = xk1 + alphak*dk;
                                            facc = obj.dcp.F.f(xacc,1);
                                            if facc < fk1
                                                if obj.verbose == 1
                                                    fprintf('accelerated: reduced %17.3e  moved %17.3e \n',facc-fk1,norm(xacc-xk1));
                                                end
                                                xk1 = xacc;
                                                fk1 = facc;
                                            end
                                        end
                                end
                            end
                        else
                            fk1 = obj.dcp.F.f(xk1,1);
                        end
                        
                        % compute errors
                        normx = norm(xk1-xk);
                        normf = abs(fk1-fk);
                        % display iterations of dca if active
                        if (obj.verbose == 1)
                            fprintf('%5d %19.5e %17.3e %13.3e\n',obj.iter,fk1,normx,normf);
                        end
                        % plotting if actived
                        if (obj.plot==1)
                            myplotf(fk,fk1,obj.iter,obj.plotlinetype);
                        end
                        % check stopping
                        if (normx < obj.tolx*(1+norm(xk1)) || normf < obj.tolf*(1+abs(fk1)))
                            if (obj.verbose == 1)
                                fprintf('------------------------------------------------------------\n');
                            end
                            obj.fopt = fk1;
                            obj.xopt = xk1;
                            status = setstatus(toc(cputime),obj.iter,0,'Successfully solved.');
                            return;
                        end
                        
                        xk = xk1;
                        fk = fk1;
                    end
            end
            % maxiter exceed
            if (obj.verbose == 1)
                fprintf('------------------------------------------------------------\n');
            end
            obj.fopt = fk1;
            obj.xopt = xk1;
            status = setstatus(toc(cputime),obj.iter,1,'Max interation exceed.');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function myplotf(fk,fk1,iter,plotline)
hold on;
if iter == 1
    title('DCA iterations');
    xlabel('Iterations');
    ylabel('Objectives');
end
if iter>1
    plot([iter-1,iter], [fk,fk1],plotline);
    drawnow
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting dca solution status
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = setstatus(timer,iter,flag,info)
status.time = timer;
status.iter = iter;
status.avgt = timer/iter;
status.flag = flag;
status.info = info;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checking active set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function activeset=checkactiveset(x)
activeset = abs(x)<=1e-7;
end
