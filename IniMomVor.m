%% IniMomVor
%  Initialises the momentum Field RhoU, RhoV of a clipped potential Vortex

function [RhoU,RhoV]=IniMomVor(Rho,Imap2,Jmap2);

% Strength factor for potential vortex
VFac=   100;
% Subtractor for vortex velocity to get velocity to zero at the boundaries
VSub=   200.0/min(Imap2,Jmap2);

% Parameters
Imidm=Imap2/2; Jmidm=Jmap2/2;

% Generate geometry arrays. Generate X, Y such that R is not zero at any
% point to avoid divistion by zero later on!
[Y,X]=meshgrid(1-Jmidm-0.5:Jmap2-Jmidm-0.5,1-Imidm-0.5:Imap2-Imidm-0.5); 
R=sqrt(X.^2+Y.^2);

% Calculate absolute velocity for a potential vortex with a solid body
% rotation inside.
AbsV=max(0., min(VFac./(R)-VSub, 2*R));

% Calculate velocity components
RhoU=Rho*AbsV.*(-Y./R);
RhoV=Rho*AbsV.*( X./R);
