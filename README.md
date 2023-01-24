# BDCA_SEiCP_SQEiCP
An accelerated DC Algorithm (BDCA) for solving Symmetric Eigenvalue Complementarity Problem (SEiCP) and Symmetric Quadratic Eigenvalue Complementarity Problem (SQEiCP).

This project is supported by the National Natural Science Foundation of China (Grant No: 11601327).

## Instruction

Run 'install' to install the toolbox on MATLAB.

This toolbox is developed based on the DCAM toolbox, see [here](https://github.com/niuyishuai/DCAM).

Two SEiCP models (the Logarithmic model - LnP and the Quadratic model - QP) are considered. BDCA and DCA are implemented for solving these two models. See the article [here](https://arxiv.org/abs/2301.09098) for more details.


## Dependencies

This toolbox depends on MOSEK for solving the convex subproblems of the QP model. The LnP model does not depend on any external solver.

The compared optimization solvers are KNITRO, FILTERSD, and MATLAB FMINCON. Please make sure that you have installed the corresponding solvers for comparison.

## Citation

```
@Misc{DCAM,
	title = {An Accelerated DC Programming Approach with Exact Line Search for The Symmetric Eigenvalue Complementarity Problem.},
	author = {Yi-Shuai Niu},	
	year = {2023},
	url = {https://arxiv.org/abs/2301.09098}
}
```

## Samples

* A simple example for generating a SEiCP is given by:
``` Matlab
n=200;
randbnd=[-1,1];
A=rand(n)*(randbnd(2)-randbnd(1))+randbnd(1);
A=(A'+A)/2;
B=eye(n);
varmu=sdpvar(1);
optimize([A+varmu*B>=0],varmu,sdpsettings('solver','mosek','verbose',0));
mu=value(varmu)+1;
A = A + mu*B;
x0 = rand(n,1);
```

Here, we use Yalmip + Mosek to find a suitable parameter $\mu$ for converting the symmetric matrix $A$ into a positive definite matrix by solving the SDP: 
$$\text{min}_{\mu} \lbrace \mu : A + \mu B \succeq 0\rbrace$$

* A simple example for solving LnP model using BDCA
```
% Create a dc function object
dcf=dcfunc;
dcf.f=@(x,opt)fobj_eval_LnP(x,A,B,opt);

% Create a dc problem object
mydcp=dcp(dcf,[]);

% Test BDCA solver
mydca = dca(mydcp,x0);
mydca.A=A;
mydca.B=B;
condA=cond(A);
condB=cond(B);
mydca.etabar=n;
mydca.Lg=mydca.etabar; 
mydca.tolf=0;
mydca.tolx=1e-8;
mydca.fistatol=1e-6;
mydca.maxiter=1e+4;
mydca.verbose = 1;
mydca.linesearch=1;
mydca.linesearch_type='exact';
mydca.model='LnP';

% Solve model using BDCA
status=mydca.optimize();

% Access solution
mydca.xopt
mydca.fopt
status.time
status.iter
status.avgt
```

`mydca.xopt`: optimal solution

`mydca.fopt`: optimal value

`status.time`: computing time

`status.iter`: number of iterations

`status.avgt`: average CPU time per iteration

See more test samples in the folder 'tests'.

## Available Dataset
Two datasets of EiCP: RANDEICP, NEP. 

One dataset of QEiCP: RANDQEICP.

See `GEN_NEP.m`, `GEN_RANDEICP.m`, `GEN_RANDQEICP.m` to generate more datasets.

## License

Released under MIT license

## Contact

niuyishuai82@hotmail.com
