%% PoissonSolver
%  Solve (approximately) the Poisson equation for the pressure correction.
%  Find a field PHI such that div(grad(PHI)) = div
%  This is a primitive, iterative solver - better solvers are available!

function [PHI]=PoissonSolver(div,Imap2,Jmap2,Dx)

% Number of iterations. More is more accurate, less is faster...
nIt = 20;

% Initialisation
PHI(1:Imap2,1:Jmap2)=0;
Ul(1:Imap2,1:Jmap2)=0.0; % Local velocities
Vl(1:Imap2,1:Jmap2)=0.0; % Local velocities

% Loop over iteration steps. 
for n=1:nIt
    % Calculate velocities induced by Phi
    Ul(1:Imap2-1,1:Jmap2)=(PHI(2:Imap2,1:Jmap2)-PHI(1:Imap2-1,1:Jmap2));
    Vl(1:Imap2,1:Jmap2-1)=(PHI(1:Imap2,2:Jmap2)-PHI(1:Imap2,1:Jmap2-1));
    % Calculate local divergence from Ul, Vl
    [divl]=calcDivergence(Ul,Vl,Dx,Imap2,Jmap2);
    
%    disp('Debug-Out');
%    imagesc(divl); pause(0.1); imagesc(div); pause(0.1);
    
    % Update PHI accordingly
    PHI = PHI + 0.2499*(div + divl);
    PHI(1,:)=PHI(2,:); PHI(Imap2,:)=PHI(Imap2-1,:);
    PHI(:,1)=PHI(:,2); PHI(:,Jmap2)=PHI(:,Jmap2-1);
end