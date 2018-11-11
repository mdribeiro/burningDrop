%% calcPresMom
%  calculates the pressure induced chang in momentnum
%  DeltaRhoU and DeltaRhoV at the cell center after a 
%  timestep of length DeltaT

function [DeltaRhoU, DeltaRhoV]=calcPresMom(P,DeltaT,Rho,DeltaX)

% Initialisation
[Ima,Jma]=size(P);
DeltaRhoU=zeros(Ima,Jma);
DeltaRhoV=zeros(Ima,Jma);

% Calculation
DeltaRhoU(2:Ima-1,:)=DeltaT/(2*Rho*DeltaX)...
                    *(P(1:Ima-2,:)-P(3:Ima,:));
DeltaRhoV(:,2:Ima-1)=DeltaT/(2*Rho*DeltaX)...
                    *(P(:,1:Ima-2)-P(:,3:Ima));