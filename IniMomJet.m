%% IniMomJet
%  Initialises the momentum Field RhoU, RhoV of a jet

function [RhoU,RhoV]=IniMomJet(Rho,Imap2,Jmap2,JJet)

% Parameters
Imidm=Imap2/2; Jmidm=Jmap2/2; Jmidp=Jmidm+1;

Mag = 1.0;

% Momentum in axial direction
RhoU(1:Imap2,1:Jmap2) = 0;
RhoU(1:Imap2,Jmidm-floor(JJet/2)-1:Jmidp+ceil(JJet/2)+1) = Rho*Mag/2;
RhoU(1:Imap2,Jmidm-floor(JJet/2):Jmidp+ceil(JJet/2)) = Rho*Mag;
% RhoU(1:Imap2,34:47) = Rho*0;

RhoU(1:Imap2,38:43) = Rho*0;
RhoV(1:Imap2,38:43) = Rho*0;

RhoU(50:Imap2,1:Jmap2) = Rho*Mag/2;

% Momentum in lateral direction. 
RhoV(1:Imap2,1:Jmap2) =0.0;

