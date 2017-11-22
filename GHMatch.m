% Author:       Soledad Villar
% Filename:     GHMatch.m
% Last edited:  May 22nd 2017
% Description:  GHMatch algorithm from [1]
% Parameters:       
%               D1: pairwise distance matrix               
%               D2: pairsise distance matrix
%               T: number of itereations 
%               sigma: Optimization parameter for augmented Lagrangiang, see [1]
%               mu: Optimization parameter for augmented Lagrangiang, see [1]
% Outputs:
%               map: correspondence between D1 and D2
%               Z: intermediate correspondences of the algorithm
%               feas: feasibilty values
%               obj: optimization objective values
% 
% References: 
% 
% [1] Villar, Bandeira, Blumberg, Ward. A polynomial-time relaxation of the 
%     Gromov-Hausdorff distance (https://arxiv.org/pdf/1610.05214.pdf)
% -------------------------------------------------------------------------
function [map, Z, feas, obj]=GHMatch(D1, D2,T, sigma, mu)
tic
n=size(D1,2);

%construct the objective function
gamma=zeros(n^2,n^2);
for i=1:n
    for j=1:n
        for k=1:n
            for L=1:n
                gamma((i-1)*n+j,(k-1)*n+L) = abs(D1(i,k)-D2(j,L));
            end
        end
    end
end
C=gamma;

%construct the constraints
[A,b]=constraints(n);
m=size(A,1); %number of constraints

%initial guess for Y
Y=ones(n^2,1)/2;

%initial multiplier
lambda=ones(m,1);
N=n^2;

%store the iterations
Z=zeros(N,T);
Z(:,1)=Y;

%feasibility
feas=zeros(T,1);
feas(1)= norm(A*Y-b);
obj=zeros(T,1);
obj(1,1)=trace(C*Y*Y');

%construct the optimization problem
options = optimoptions('fmincon','Algorithm','trust-region-reflective','GradObj','On', 'Hessian','user-supplied','Display','Off');
options.FunctionTolerance=1e-4;
problem.options = options;
problem.lb=zeros(N,1);
problem.ub=ones(N,1);
problem.solver = 'fmincon';

for i=2:T
    %update problem's parameters
    problem.x0 = Y;
    problem.objective = @(Y)lagrangian(C, A,b, Y, lambda, sigma);
    Y = fmincon(problem);
    
    %update feasibility
    feas(i,1)=norm(A*Y-b);
    
    %update multipliers
    lambda=lambda-sigma*(A*Y - b);
    sigma=mu*sigma;
    Z(:,i)=Y;
    obj(i,1)=trace(C*Y*Y');
end

%find the corresponding map
map=zeros(n,1);
for i=1:n
    [a,map(i,1)]=max(Z(1+(i-1)*n:i*n,T));
end

%number of unmatched nodes
if length(unique(map))<n
    error=[length(unique(map)),n]
end

%plot(feas(4:T))
time=toc
end


function [A,b] = constraints(N)
A=zeros(2*N,N^2);
b=ones(2*N,1);
for i=1:N
    A(i,(i-1)*N+1:i*N)=ones(1,N);
    for j=1:N
        A(i+N, (j-1)*N+i)=1;
    end
end
end


function [f, g, H] = lagrangian(C, A,b, Y, lambda, sigma)
f=trace(C*(Y*Y'))-lambda'*(A*Y-b)+ 0.5*sigma*norm(A*Y-b,2)^2;

if nargout > 1 % gradient required
    grad=2*C*Y-A'*lambda + sigma*A'*(A*Y-b);
    [n,k]=size(Y);
    [m,N]=size(A);
    g=grad;
    if nargout>2
        H=2*C + sigma*(A'*A);
    end
end
end