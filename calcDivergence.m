%% calcDivergence
%  calculates the divergence of the velocity field, div = du/dx + dv/dy

function [div]=calcDivergence(U,V,Dx,Imap2,Jmap2)

% Initialisation
div(1:Imap2,1:Jmap2) = 0.0;

% Calculate divergence inside domain
div(2:Imap2-1,2:Jmap2-1) = 1/1 ...
                         * ( U(2:Imap2-1,2:Jmap2-1)-U(1:Imap2-2,2:Jmap2-1)...
                           + V(2:Imap2-1,2:Jmap2-1)-V(2:Imap2-1,1:Jmap2-2) );