%% IniMomDiv
%  Initialises the momentum Field RhoU, RhoV to create a divergence.

function [RhoU,RhoV]=IniMomDiv(Rho,Imap2,Jmap2);

% Momentum in axial direction
RhoU(1:Imap2,1:Jmap2) = 0.0;
RhoU(ceil(Imap2/3):floor(2*Imap2/3),1:Jmap2) = Rho*1.0;

% Momentum in lateral direction
RhoV(1:Imap2,1:Jmap2) =0.0;

