%% calcPresVel
%  calculates the pressure induced change in velocity
%  DeltaU and DeltaV at the cell surfaces after a 
%  timestep of length DeltaT

function [DeltaU, DeltaV]=calcPresVel(P,DeltaT,Rho,DeltaX)

% Initialisation
[Ima,Jma]=size(P);
DeltaU=zeros(Ima,Jma);
DeltaV=zeros(Ima,Jma);

% Calculation
DeltaU(2:Ima-1,:)=DeltaT/(2*Rho^2*DeltaX)...
                 *(P(1:Ima-2,:)-P(3:Ima,:));
DeltaV(:,2:Ima-1)=DeltaT/(2*Rho^2 *DeltaX)...
                 *(P(:,1:Ima-2)-P(:,3:Ima));