%% mainPassive
%  is a simple solver for the transport of a passive scalar (enthalpy,
%  mass-fraction, momentum component in third direction, ...) in 2D
%  (C) Andreas Kempf, 21.5.2013, code for "Numerics and Fluid Flow Simulation"
%  course at the University of Duisburg-Essen

%% Preparation
clc; clear; hold off;

%% Set parameters
CFL =     0.5    ;% [  --  ] CFL number to calculate timestep from
Rho =     1.2    ;% [kg/m^3] Density of the fluid (constant throughout flow-field)
Mu  =     1.0E-5 ;% [ Pa.s ] Dynamic viscosity
Dx  =     1.0E-3 ;% [  m   ] Size of computational cells in x and y direction
Nt  = 20000      ;% [  --  ] Number of time-steps
Ima =   200      ;% [  --  ] Number of real cells (not ghost-cells) in x-direction
Jma =   200      ;% [  --  ] Number of real cells (not ghost-cells) in y-direction
Case=     1      ;% [      ] 1: Jet, 2: Clipped potential vortex
JJet=    16      ;% [  --  ] Width of jet in cells (for case 1)

%% Initialise useful helper parameters
Imap2 = Ima+2; Jmap2 = Jma+2;       % Ima (I maximum), Imap2 (Ima plus 2)
Imidm = Imap2/2; Imidp = Imap2/2+1; % Jmidm (J middle minus)
Jmidm = Jmap2/2; Jmidp = Jmap2/2+1; % Jmidm (J middle minus)
Ifim=1; Ifi=2; Ifip=3; Ilam=Imap2-2; Ila=Imap2-1; Ilap=Imap2; % Ifim (I first minus)
Jfim=1; Jfi=2; Jfip=3; Jlam=Jmap2-2; Jla=Jmap2-1; Jlap=Jmap2; % Jlap (J last plus)

%% Initialise momentum field

% Set up momentum field for cases
if (Case == 1); [RhoU,RhoV]=IniMomJet(Rho,Imap2,Jmap2,JJet); 
elseif (Case == 2); [RhoU,RhoV]=IniMomVor(Rho,Imap2,Jmap2);
end
% Initialise velocity field
[U,V]=mom2vel(RhoU,RhoV,Rho);

%% Initialise scalar field
RhoPhi(1:Imap2,1:Jmap2)=0.0;
for j=1:Jmap2;
    RhoPhi(:,j)=Rho*0.5*( 1.0+sin((1:Imap2)/(2.0*pi)) ); 
end

%% Time-integration
disp('Start time-stepping')
for n=1:Nt
    
    % Calculate time-step width
    det = CFL * Dx / max(max(max(U)),max(max(V)));

    % Calculate fluxes
    [FluxConX, FluxConY]=calcFluxConUSCDS(RhoPhi,U,V,Dx,det);
    [FluxDifX, FluxDifY]=calcFluxDif(RhoPhi,Dx,Mu/Rho);

    % Apply fluxes
    RhoPhi(Ifi:Ila,Jfi:Jla) = RhoPhi(Ifi:Ila,Jfi:Jla) ...
                             + det/(Dx^3)*( FluxConX(Ifim:Ilam,Jfi:Jla)-FluxConX(Ifi:Ila,Jfi:Jla)...
                                          + FluxConY(Ifi:Ila,Jfim:Jlam)-FluxConY(Ifi:Ila,Jfi:Jla)) ...
                             - det/(Dx^3)*( FluxDifX(Ifim:Ilam,Jfi:Jla)-FluxDifX(Ifi:Ila,Jfi:Jla)...
                                          + FluxDifY(Ifi:Ila,Jfim:Jlam)-FluxDifY(Ifi:Ila,Jfi:Jla));
                                      
    % Boundary treatment
    if (Case == 1)
      RhoPhi(Ilap,:) = RhoPhi(Ila,:);
      RhoPhi(:,Jfim) = RhoPhi(:,Jfi);
      RhoPhi(:,Jlap) = RhoPhi(:,Jla);
    end
    
    % Output during calculation, only every N timesteps to save time
    if (mod(n,100) == 0); pause(1.0);
        imagesc(RhoPhi/Rho);
    end
end
disp('Time stepping finished');


%% Post-processing
%  The imagesc routine sets the origin of the coordinate system at the
%  top-left corner, x pointing downwards, y to the right.
imagesc(U); pause();
imagesc(V); pause();
imagesc(RhoPhi/Rho); pause();