


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

beta=model.params(1);
sigma=model.params(2);
eta=model.params(3);
chi=model.params(4);
delta=model.params(5);
alpha=model.params(6);
rho=model.params(7);
zbar=model.params(8);

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
    %For capital stock k
    kbar= nbar.*((1/beta-1+delta)/(alpha))^(1/(alpha-1));
    kpas=0.01;
    nk=5;


    k=[kbar-(floor(nk/2):-1:1).*kpas, kbar, kbar+(1:floor(nk/2)).*kpas];
    nk=size(k,2);
    kmin=min(k);
    kmax=max(k);
    %For aggregate productivity
    
    zpas=0.025;
    nz=5;

    z=[zbar-(floor(nz/2):-1:1).*zpas, zbar, zbar+(1:floor(nz/2)).*zpas];
    nz=size(z,2);
    zmin=min(z);
    zmax=max(z);
    
    %Compute the grid
    grid=gridmake(z',k');
    
    
%% Iterate
%Definied interpolators
order=[nz nk];
gridmin=[zmin kmin];
gridmax=[zmax kmax];
cdef=fundefn('lin',order,gridmin,gridmax);

%Deterministic case
e=zeros(25,2);


% Convergence criterion
tol=1e-3;
maxiteration=100;
options=optimset('MaxIter',1E5,'MaxFunEvals',1E5);

%Initialisation
iinit=delta*kbar;
ninit=0.33;
xinit=[iinit*ones(25,1) ninit*ones(25,1)];

x=xinit;


iteration=1;
converge=0;

while converge==0,
    
    
    
    
    snext=G_rbc(grid,x,e,model);
    
    [coeff,B]=funfitxy(cdef,grid,x);
    xnext=funeval(coeff, cdef, snext);

    
    
    %x_up=fsolve(@(xt) F_rbc(grid,xt,e,G_rbc(grid,xt,e,model),xnext,model),x,options);
    x_up=fsolve(@(xt) F_rbc(grid,xt,e,snext,xnext,model),x,options);
    
    err=sum(sum(abs(x-x_up)))
    
   
    
    if (err < tol && iteration < maxiteration);
        converge=1;
    end;
    
    x=x_up
    iteration = iteration+1;
    display([iteration err]);
end;

