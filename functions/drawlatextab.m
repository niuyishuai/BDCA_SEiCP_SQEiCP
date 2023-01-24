%% format 1 for EiCP to QP & LnP models (RAND and NEP)
% Prob | DCA (lambda,CPU,IT,c) | BDCA (lambda,CPU,IT,c) | FMINCON
% (lambda,CPU,c) | KNITRO (lambda,CPU,c) | FILTERSD (lambda,CPU,c)
% average | CPU IT c | CPU IT c | CPU c | CPU c| CPU c
algolst={'DCA','BDCA','FMINCON','KNITRO','FILTERSD'};
model='LnP'; % QP|LnP

% RANDEICP dataset
range_n=[50,100,200,400,600,800];
bndlst={[-1,1],[-10,10]};
nbalgos=length(algolst);
avgdata=zeros(nbalgos,3);
nbdata = 0;
for bnd=bndlst
    for n=range_n
        fprintf('RANDEICP($%d,%d,%d$)',bnd{1}(1),bnd{1}(2),n);
        for algo=algolst
            load(sprintf('RESULT\\%s_%s_RANDEICP_%d_%d_%d.mat',algo{1},model,bnd{1}(1),bnd{1}(2),n));
            c = -fix(log10(sum(err)));
            switch algo{1}
                case 'DCA'
                    fprintf(' & $%.4f$ & $%.3f$ & $%d$ & $%d$',lambda-mu,cputime,iters,c);
                    avgdata(1,:) = avgdata(1,:) + double([cputime,iters,c]);
                case 'BDCA'
                    fprintf(' & $%.4f$ & $%.3f$ & $%d$ & $%d$',lambda-mu,cputime,iters,c);
                    avgdata(2,:) = avgdata(2,:) + double([cputime,iters,c]);
                case 'FMINCON'
                    fprintf(' & $%.4f$ & $%.3f$ & $%d$',lambda-mu,cputime,c);
                    avgdata(3,:) = avgdata(3,:) + double([cputime,c,0]);
                case 'KNITRO'
                    fprintf(' & $%.4f$ & $%.3f$ & $%d$',lambda-mu,cputime,c);
                    avgdata(4,:) = avgdata(4,:) + double([cputime,c,0]);
                case 'FILTERSD'
                    fprintf(' & $%.4f$ & $%.3f$ & $%d$',lambda-mu,cputime,c);
                    avgdata(5,:) = avgdata(5,:) + double([cputime,c,0]);
            end
        end
        nbdata=nbdata+1;
        fprintf('\\\\\n');
    end
    fprintf('\\midrule\n');
end
% NEP dataset
NEPlst=dir('NEP\\*.mtx');
for i=1:length(NEPlst)
        [~,fname]=fileparts(NEPlst(i).name);
        fprintf('NEP-%s',fname);
        for algo=algolst
            load(sprintf('RESULT\\%s_%s_NEP_%s.mtx.mat',algo{1},model,fname));
            c = -fix(log10(sum(err)));
            switch algo{1}
                case 'DCA'
                    fprintf(' & $%.4f$ & $%.3f$ & $%d$ & $%d$',lambda-mu,cputime,iters,c); 
                    avgdata(1,:) = avgdata(1,:) + double([cputime,iters,c]);
                case 'BDCA'
                    fprintf(' & $%.4f$ & $%.3f$ & $%d$ & $%d$',lambda-mu,cputime,iters,c);
                    avgdata(2,:) = avgdata(2,:) + double([cputime,iters,c]);
                case 'FMINCON'
                    fprintf(' & $%.4f$ & $%.3f$ & $%d$',lambda-mu,cputime,c);
                    avgdata(3,:) = avgdata(3,:) + double([cputime,c,0]);
                case 'KNITRO'
                    fprintf(' & $%.4f$ & $%.3f$ & $%d$',lambda-mu,cputime,c);
                    avgdata(4,:) = avgdata(4,:) + double([cputime,c,0]);
                case 'FILTERSD'
                    fprintf(' & $%.4f$ & $%.3f$ & $%d$',lambda-mu,cputime,c);
                    avgdata(5,:) = avgdata(5,:) + double([cputime,c,0]);
            end
        end
        nbdata=nbdata+1;
        fprintf('\\\\\n');
end
fprintf('\\midrule\n');
avgdata=avgdata/nbdata;
fprintf('avg & & $%.3f$ & $%d$ & $%d$ & & $%.3f$ & $%d$ & $%d$ & & $%.3f$ & $%d$ & & $%.3f$ & $%d$ & & $%.3f$ & $%d$\\\\\n',...
avgdata(1,1),int32(avgdata(1,2)),int32(avgdata(1,3)),...
avgdata(2,1),int32(avgdata(2,2)),int32(avgdata(2,3)),...
avgdata(3,1),int32(avgdata(3,2)),...
avgdata(4,1),int32(avgdata(4,2)),...
avgdata(5,1),int32(avgdata(5,2)));
fprintf('\\bottomrule\n');

%% format 2 for QEiCP to QP & LnP models (RAND)
% Prob | DCA (lambda,CPU,IT,c) | BDCA (lambda,CPU,IT,c) | FMINCON
% (lambda,CPU,c) | KNITRO (lambda,CPU,c) | FILTERSD (lambda,CPU,c)
% average | CPU IT c | CPU IT c | CPU c | CPU c| CPU c
algolst={'DCA','BDCA','FMINCON','KNITRO','FILTERSD'};
model='QP'; % QP|LnP

% RANDQEICP dataset
range_n=[50,100,200,400,600]; 
range_d=[0.05,0.1,0.5,0.7,0.9];
nbalgos=length(algolst);
avgdata=zeros(nbalgos,3);
nbdata = 0;
for d=range_d
    for n=range_n
        fprintf('RANDQEICP($%d\\%%,%d$)',100*d,n);
        for algo=algolst
            load(sprintf('RESULT\\%s_%s_RANDQEICP_%.2f_%d.mat',algo{1},model,d,n));
            c = -fix(log10(sum(err)));
            switch algo{1}
                case 'DCA'
                    fprintf(' & $%.4f$ & $%.3f$ & $%d$ & $%d$',lambda-mu,cputime,iters,c);
                    avgdata(1,:) = avgdata(1,:) + double([cputime,iters,c]);
                case 'BDCA'
                    fprintf(' & $%.4f$ & $%.3f$ & $%d$ & $%d$',lambda-mu,cputime,iters,c);
                    avgdata(2,:) = avgdata(2,:) + double([cputime,iters,c]);
                case 'FMINCON'
                    fprintf(' & $%.4f$ & $%.3f$ & $%d$',lambda-mu,cputime,c);
                    avgdata(3,:) = avgdata(3,:) + double([cputime,c,0]);
                case 'KNITRO'
                    fprintf(' & $%.4f$ & $%.3f$ & $%d$',lambda-mu,cputime,c);
                    avgdata(4,:) = avgdata(4,:) + double([cputime,c,0]);
                case 'FILTERSD'
                    fprintf(' & $%.4f$ & $%.3f$ & $%d$',lambda-mu,cputime,c);
                    avgdata(5,:) = avgdata(5,:) + double([cputime,c,0]);
            end
        end
        nbdata=nbdata+1;
        fprintf('\\\\\n');
    end
    fprintf('\\midrule\n');
end
avgdata=avgdata/nbdata;
fprintf('avg & & $%.3f$ & $%d$ & $%d$ & & $%.3f$ & $%d$ & $%d$ & & $%.3f$ & $%d$ & & $%.3f$ & $%d$ & & $%.3f$ & $%d$\\\\\n',...
avgdata(1,1),int32(avgdata(1,2)),int32(avgdata(1,3)),...
avgdata(2,1),int32(avgdata(2,2)),int32(avgdata(2,3)),...
avgdata(3,1),int32(avgdata(3,2)),...
avgdata(4,1),int32(avgdata(4,2)),...
avgdata(5,1),int32(avgdata(5,2)));
fprintf('\\bottomrule\n');


