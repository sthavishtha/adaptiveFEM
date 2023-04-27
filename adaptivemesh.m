% finite element code with P1-C0 linear basis functions for adaptive mesh
% based on Long Chen's IFEM MATLAB package
% (c) Sthavishtha

close all; 
clear variables;

%% Parameters 
maxIt = 4; 
N = zeros(maxIt,1);
h = zeros(maxIt,1);
errL2 = zeros(maxIt,1); 
errH1 = zeros(maxIt,1); 
erruIuh = zeros(maxIt,1);

%% Generate an initial mesh 
[node,elem] = squaremesh([0 1 0 1], 0.25);
bdFlag = setboundary(node,elem,'Dirichlet','x==0','Neumann','(x == 1) | (y == 0) | (y == 1)');

mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
option.L0 = 1;
option.maxIt = 200;
option.maxN = 1e4;
option.theta = 0.4;
option.plotflag = 1;
option.rateshift = 20;
option.elemType = 'P1';
option.refType = 'red'; %red %bisect
option.estType = 'residual'; %residual %recovery
option.dispspace = 1;
option.markType = 'L2';

%% Get the data of the pde
pde = powerlaw;

%% Finite Element Method        
err = afemPoisson(mesh,pde,option);

%% Plot convergence rates and display error table
showrate2(err.N,err.H1,20,'k-*','||Du-Du_h||',err.N,err.eta,20,'-+','eta');

function pde = powerlaw

% change the value of r
r = 1.6;

pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'g_N',@g_N,'Du',@Du);

    % load data (right hand side function)
    function rhs =  f(p)
    x = p(:,1); y = p(:,2);
    rhs =  r*(1. - r)*x.^(r-2.);
    end
    % exact solution
    function u =  exactu(p)
    x = p(:,1); y = p(:,2);
    u =  x.^r;
    end
    % Dirichlet boundary condition
    function u =  g_D(p)
    u =  exactu(p);
    end
    % Neumann boundary condition
    function f =  g_N(p) 
    f = zeros(size(p,1),1);
    x = p(:,1); y = p(:,2);
    downbd = (abs(y - 0)<eps);
    upbd = (abs(y - 1)<eps);
    rightbd = (abs(x - 1)<eps);
    f(rightbd) = r*ones(size(p(rightbd),1),1);
    f(downbd) = zeros(size(p(downbd),1),1);
    f(upbd) = zeros(size(p(upbd),1),1);
    end    
    % Derivative of the exact solution
    function uprime =  Du(p)
    x = p(:,1); y = p(:,2);
    uprime(:,1) = r*x.^(r-1.);
    uprime(:,2) = 0.;
    end
end