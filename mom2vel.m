%% mom2vel
%  Calculate the velocity on cell surfaces from the 
%  momentum at the cell centre.

function [u,v]=mom2vel(rhou,rhov,rho)

% Initialisation
[Ima,Jma]=size(rhou);
u=zeros(Ima,Jma);
v=zeros(Ima,Jma);

% Calculation
u(1:Ima-1,:) = 0.5/rho * (rhou(1:Ima-1,:) + rhou(2:Ima,:));
v(:,1:Jma-1) = 0.5/rho * (rhov(:,1:Jma-1) + rhov(:,2:Jma));

