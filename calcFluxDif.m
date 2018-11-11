%% Caculate the Diffusive fluxes
%  fluxDifX, fluxDifY over the cell surfaces as a function of
%  a conserved scalar field Phi, the diffusivity D and a 
%  grid-spacing of deltaX

function [fluxDifX,fluxDifY]=calcFluxDif(Phi,deltaX,D)

% Initialisation
[Ima,Jma]=size(Phi);
fluxDifX=zeros(Ima,Jma);
fluxDifY=zeros(Ima,Jma);

% Calculation
fluxDifX(1:Ima-1,:) = D*deltaX*(Phi(2:Ima,:)-Phi(1:Ima-1,:));
fluxDifY(:,1:Jma-1) = D*deltaX*(Phi(:,2:Jma)-Phi(:,1:Jma-1));
                