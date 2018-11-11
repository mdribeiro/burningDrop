%% Caculate the convective fluxes
%  fluxConX, fluxConY over the cell surfaces as a function of
%  a conserved scalar field Phi from the velocities u and v 
%  for a grid-spacing of deltaX

function [fluxConX,fluxConY]=calcFluxConUSCDS(Phi,u,v,deltaX,dt)

% Initialisation
[Ima,Jma]=size(Phi);
fluxConX=zeros(Ima,Jma);
fluxConY=zeros(Ima,Jma);

% Calculate weigths
ww(1:Ima-1,:)=0.5*(u(1:Ima-1,:)*dt/deltaX + 1);
we(2:Ima,:)  =1.0-ww(1:Ima-1,:);
ws(:,1:Jma-1)=0.5*(v(:,1:Jma-1)*dt/deltaX + 1);
wn(:,2:Jma)  =1.0-ws(:,1:Jma-1);

% Calculation
fluxConX(1:Ima-1,:) = deltaX^2 * u(1:Ima-1,:)...
                    .* ( ww(1:Ima-1,:).*Phi(1:Ima-1,:)...
                       + we(2:Ima,:)  .*Phi(2:Ima,:) );
fluxConY(:,1:Jma-1) = deltaX^2 * v(:,1:Jma-1)...
                    .* ( ws(:,1:Jma-1).*Phi(:,1:Jma-1)...
                       + wn(:,2:Jma)  .*Phi(:,2:Jma) );
                