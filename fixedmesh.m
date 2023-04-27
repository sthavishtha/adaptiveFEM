% finite element code with P1-C0 linear basis functions for uniform mesh
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
errL2_ratio = zeros(maxIt,1);
errH1_semi_ratio = zeros(maxIt,1);

%% Generate an initial mesh 
[node,elem] = squaremesh([0 1 0 1], 0.25);
bdFlag = setboundary(node,elem,'Dirichlet','x==0','Neumann','(x == 1) | (y == 0) | (y == 1)');

%% Get the data of the pde
pde = powerlaw;

%% Finite Element Method        
for k = 1:maxIt
    % refine mesh
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    % solve the equation
    [soln,eqn,info] = Poisson(node,elem,bdFlag,pde);
    uh = soln.u;
    % record and plot
    N(k) = size(elem,1);
    h(k) = 1./(sqrt(size(node,1))-1);
    if N(k) < 2e3 % show mesh and solution for small size
        figure(1);  showresult(node,elem,uh);    
    end
    % compute error
    uI = pde.exactu(node); % nodal interpolation
    errL2(k) = getL2error(node,elem,pde.exactu,uh);
    errH1(k) = getH1error(node,elem,pde.Du,soln.Du);
    erruIuh(k) = sqrt((uh-uI)'*eqn.Lap*(uh-uI));

    % extra computations
    L2_exact = getL2error(node,elem,pde.exactu,zeros(size(uh)), 1);
    H1_exact = getH1error(node,elem,pde.Du,zeros(size(soln.Du)), 1);
    errL2_ratio(k) = errL2(k)/L2_exact;
    H1_semi = sqrt(errH1(k)*errH1(k) - errL2(k)*errL2(k));
    H1_semi_exact = sqrt(H1_exact*H1_exact - L2_exact*L2_exact);
    errH1_semi_ratio(k) = H1_semi/H1_semi_exact;
end

%% Plot convergence rates and display error table
figure(2);
showrateh3(h,errH1,1,'-*','||Du-Du_h||',...
           h,errL2,1,'k-+','||u-u_h||', ...
           h,erruIuh,1,'m-+','||DuI-Du_h||');

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','||Du-Du_h||','||DuI-Du_h||','||u-u_h||/||u||','||Du-Du_h||/||Du||'};
disptable(colname,N,[],h,'%0.3e',errL2,'%0.5e',errH1,'%0.5e',erruIuh,'%0.5e',errL2_ratio,'%0.5e',errH1_semi_ratio,'%0.5e');

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