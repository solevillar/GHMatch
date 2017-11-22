% Author:       Soledad Villar
% Filename:     main.m
% Last edited:  May 22nd 2017
% Description:  Renders examples of algorithm GHMAtch from [1] on 3d  
%               surfaces from [2], available in [3].
% Parameters:       
%               n=25: number of points to sample              
%               T=15: Number of GHMatch iterations
%               sigma=4: Optimization parameter, see [1]
%               mu=8: Optimization parameter, see [1]
% Outputs:
%               renders a plot
% 
% References:
% 
% [1] Villar, Bandeira, Blumberg, Ward. A polynomial-time relaxation of the 
%     Gromov-Hausdorff distance (https://arxiv.org/pdf/1610.05214.pdf)
% [2] Alexander Bronstein, Michael Bronstein, and Ron Kimmel. Numerical 
%     Geometry of Non-Rigid Shapes. Springer 2008.
% [3] http://tosca.cs.technion.ac.il/book/resources_data.html
% -------------------------------------------------------------------------

%parameters
n=25; %Number of points to match
T=15; %Number of GHMatch iterations
sigma=4; %Optimization parameter, see [1]
mu=8;    %Optimization parameter, see [1]

rng(1) %seed (for randomness)

%EXAMPLE, see other options at ./ronrigid3d
s1=load('./nonrigid3d/david0.mat');
s2=load('./nonrigid3d/david7.mat');
s1=s1.surface;
s2=s2.surface;

%simplify the mesh for faster rendering
surf1=reducepatch(s1.TRIV,[s1.X,s1.Y, s1.Z], 0.2);
surf2=reducepatch(s2.TRIV,[s2.X,s2.Y, s2.Z], 0.2);


%plot surfaces
figure;
trisurf(surf1.faces, surf1.vertices(:,1),surf1.vertices(:,2), surf1.vertices(:,3), 'FaceColor', [.27 .27 .27]); 
hold on
trisurf(surf2.faces, surf2.vertices(:,1),surf2.vertices(:,2)-100, surf2.vertices(:,3), 'FaceColor', [.27 .27 .27]); 
axis equal 
axis off

%triangulation to adjacency matrix
G1=tri2graph(s1);
G2=tri2graph(s2);


%here we sample from the mesh
%we choose this sampling strategy to get uniform points from the surface
%other sampling strategies may be better
v1=zeros(n,1);
aux1=conncomp(graph(G1));
for i=1:n
    while v1(i)==0
        x=randsample(floor(size(s1.X,1)/n),1);
        if and(aux1(x)==1, sum(aux1==x)==0)
            v1(i)=x+ (i-1)*floor(size(s1.X,1)/n);
        end
    end
end

v2=v1; %same indices for simplicity


scatter3(s1.X(v1),s1.Y(v1),s1.Z(v1), 'fill', 'black')
scatter3(s2.X(v2),s2.Y(v2)-100,s2.Z(v2), 'fill', 'black');


D1=zeros(n,n);
D2=zeros(n,n);
for i=1:n
    for j=i+1:n
        D1(i,j)=graphshortestpath(G1,v1(i), v1(j), 'directed', false);
        D2(i,j)=graphshortestpath(G2,v2(i), v2(j), 'directed',false);
    end
end
D1=D1+D1';
D2=D2+D2';

%call GHMatch
[map12, Z12, feas12, obj]=GHMatch(D1, D2,T, sigma, mu);

%plot the correspondence
CM=spring(n);
[a,b]=sort(s1.Z(v1));
for i=1:n
    plot3( [s1.X(v1(b(i))); s2.X(v2(map12(b(i))))] , [s1.Y(v1(b(i))); s2.Y(v2(map12(b(i))))-100], [s1.Z(v1(b(i))); s2.Z(v2(map12(b(i))))], 'color', CM(i,:), 'LineWidth',4 );
end

