%% Caculate the convective fluxes
%  fluxConX, fluxConY over the cell surfaces as a function of
%  a conserved scalar field Phi from the velocities u and v 
%  for a grid-spacing of deltaX

function [fluxConX,fluxConY]=calcFluxConCDS(Phi,u,v,deltaX)

% Initialisation
[Ima,Jma]=size(Phi);
fluxConX=zeros(Ima,Jma);
fluxConY=zeros(Ima,Jma);

% Calculation
fluxConX(1:Ima-1,:) = deltaX^2/2 * u(1:Ima-1,:)...
                    .* (Phi(1:Ima-1,:)+Phi(2:Ima,:));
fluxConY(:,1:Jma-1) = deltaX^2/2 * v(:,1:Jma-1)...
                    .* (Phi(:,1:Jma-1)+Phi(:,2:Jma));
                