


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need a model written as:
%
% snext=G(s,x,e)
% E F(s,x,e,snext,xnext)=0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('D:\Users\Utilisateur\Documents\MATLAB\COMPECON\CEtools');

%% Parameters
none=[];
model=rbc_matlab('model',none,none,none,none,none,none,none);

%[beta, sigma, eta, chi, delta, alpha, rho, zbar, nsurk]
% model.params;

%model.params(1) = 0.9;

beta=model.params(1);
sigma=model.params(2);
eta=model.params(3);
chi=model.params(4);
delta=model.params(5);
alpha=model.params(6);
rho=model.params(7);
zbar=model.params(8);

params = model.params;

% phi= 1;
% chi = 4;
% delta = 0.025;     
% alpha = 0.33;     
% rho_z = 0.95;
% sigma= 1;
% eta= 1;
% zbar= 1;


nbar=0.33;

%% Defined the grid

    nz=40;
    nk=40;


    %For capital stock k
    kbar= nbar.*((1/beta-1+delta)/(alpha))^(1/(alpha-1));
    kpas=0.1;


    %k=[kbar-(floor(nk/2):-1:1).*kpas, kbar, kbar+(1:floor(nk/2)).*kpas];
    %nk=size(k,2);

    
    
    k = linspace(kbar*0.95,kbar*1.05,nk);
    kmin=min(k);
    kmax=max(k);
    %For aggregate productivity
    
    zpas=0.00025;

    z = linspace(0.95, 1.05, nz);
    zmin=min(z);
    zmax=max(z);
    
    
    %Compute the grid
    grid=gridmake(z',k');
    
    n_s = size(grid,1);
    
%% Iterate
%Definied interpolators
order=[nz nk];

gridmin=[zmin kmin ];
gridmax=[zmax kmax ];
cdef=fundefn('lin',order,gridmin,gridmax);

%Deterministic case
e=zeros(n_s,2);


% Convergence criterion
tol=1e-10;
maxiteration=1000;
options=optimset('MaxIter',1E5,'MaxFunEvals',1E5,'Display','iter');

%Initialisation
x_ss = model.x_ss;
s_ss = model.s_ss;
X_s = model.X{2};
iinit = x_ss(1) + (grid(:,1)-s_ss(1)) * X_s(1,1) + (grid(:,2)-s_ss(2)) * X_s(1,2);
ninit = x_ss(2) + (grid(:,1)-s_ss(1)) * X_s(2,1) + (grid(:,2)-s_ss(2)) * X_s(2,2);

xinit=[iinit, ninit];


x=xinit;


iteration=1;
converge=0;

tic;

integration='hermite';
switch integration;
    case 'naive';
    N_shocks = 10;
    Sigma = [[0.0025^2]];
    epsilons = normrnd(0,Sigma, N_shocks,1);
    weights = ones(N_shocks,1)/N_shocks;

    case 'hermite';
    N_shocks = 5;
    [epsi,w]=hernodes(N_shocks);
    Sigma = [[0.0025]];
    epsilons = 0+Sigma*sqrt(2)*epsi;
    weights=w;
end;


disp('Starting iterations');
while converge==0 && iteration < maxiteration
    
    [coeff,B]=funfitxy(cdef, grid, x);
    
    fobj = @(xt) step_residuals(grid, xt, epsilons, weights, model.params, @F_rbc, @G_rbc, model, coeff, cdef);
    
    %x_up=fsolve(fobj, x, options);
    x_up = newton_solver(fobj, x);
    
    err=sum(sum(abs(x-x_up)))
    
    if (err < tol);
        converge=1;
    end;
    
    x=x_up;
    iteration = iteration+1;
    %display([iteration err]);
end;

toc;
