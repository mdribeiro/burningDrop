%% mainPassive
%  is a simple solver for the burning of a single fuel drop in 2D
%  (C) Mateus Ribeiro, 2016, based on (C) Andreas Kempf, 21.5.2013, code for "Numerics and Fluid Flow Simulation"
%  course at the University of Duisburg-Essen

%% Preparation
clc; clear; %hold off;

get(0,'Factory');
set(0,'defaultfigurecolor',[1 1 1]);
FigHandle = figure;
set(FigHandle,'Position',[20,100,1550,600]);

%% Set parameters
CFL =     0.2    ;% [  --  ] CFL number to calculate timestep from
Rho =     1.2    ;% [kg/m^3] Density of the fluid (constant throughout flow-field)
Mu  =     4*1.0E-5 ;% [ Pa.s ] Dynamic viscosity  
%Sc  =     0.7   ;% [  --  ] Schmidt-number; ratio of viscosity / diffusivity
Sc  =     0.2    ;% [  --  ] Reduced for visualisation. Correct value: 0.7
Dx  =     1.0E-3 ;% [  m   ] Size of computational cells in x and y direction
Nt  = 80000      ;% [  --  ] Number of time-steps
Ima =   250      ;% [  --  ] Number of real cells (not ghost-cells) in x-direction  Org = 128
Jma =    80      ;% [  --  ] Number of real cells (not ghost-cells) in y-direction  Org = 80
Case=     1      ;% [  --  ] 1: Jet, 2: Clipped potential vortex, 3: Divergence field
JJet=    Jma-2      ;% [  --  ] Width of jet in cells (for case 1)

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
elseif (Case == 3); [RhoU,RhoV]=IniMomDiv(Rho,Imap2,Jmap2);
end
% Initialise velocity field
[U,V]=mom2vel(RhoU,RhoV,Rho);

%% Droplet initial properties and fuel propeties
r0 = 10; %[m]
d0 = 2*r0*0.001;
V0 = (4/3)*pi*(r0*0.001)^3; %[m3]
rho0 = 690; %[kg/m3]
rhog = 3.5; %[kg/m3]
m0 = rho0*V0; %[kg]
Re = Rho*d0*max(max(U))/Mu;
B = 9.70;       % Spalding Transfer Number for Iso-octane
%B = 4;
cp = 2.145;     % [kJ/(kg.K)] Sensible heat of Iso-octane

%% Initialise scalar field
RhoPhi(1:Imap2,1:Jmap2)=700*Rho;
RhoO2(1:Imap2,1:Jmap2)=0.232*Rho;
RhoN2(1:Imap2,1:Jmap2)=0.768*Rho;
RhoFuel(1:Imap2,1:Jmap2)=0*Rho;

% In order to conform the the geometry (droplet part)
x0 = 31;
y0 = Jma/2+1;
for i=21:41
    for j=(Jma/2+1-d0*1000/2):(Jma/2+1+d0*1000/2)
        rij = floor(  sqrt(  (i-x0)^2 + (j-y0)^2  )  );
        if (rij <= (r0-1))
            RhoU(i,j) = 0;
            RhoV(i,j) = 0;
            RhoFuel(i,j) = 0;
            RhoO2(i,j) = 0;
            RhoN2(i,j) = 0;
        end

    end
end 
% Flow channel part
RhoU(1:x0,38:43) = Rho*0;
RhoV(1:x0,38:43) = Rho*0;
RhoFuel(1:x0,38:43) = Rho*0;
RhoO2(1:x0,38:43) = Rho*0;
RhoN2(1:x0,38:43) = Rho*0;
RhoPhi(1:x0,38:43) = Rho*400;
%imagesc(RhoU'/Rho);

% Waves
%for j=1:Jmap2; RhoPhi(:,j)=Rho*0.5*( 1.0+sin((1:Imap2)/(2.0*pi)) ); end
% For jet:
% RhoPhi(1:Imap2,Jmidm-floor(JJet/2)-1:Jmidp+ceil(JJet/2)+1) = Rho*0.5;
% RhoPhi(1:Imap2,Jmidm-floor(JJet/2):Jmidp+ceil(JJet/2)) = Rho*1.0;

%% Time-integration
disp('Start time-stepping')
cont1 = 0;
T = 0;
for n=1:Nt
    %% I. Initial work: calculate time-step width
    det = CFL * Dx / max(max(max(U)),max(max(V)));
    
    T = T + det;
    
    cont1 = cont1 + 1;

    %% 1. Calculate and apply fluxes
    % Scalar RhoPhi Phi = Temperature
    %[FluxConX, FluxConY]=calcFluxConCDS(RhoPhi,U,V,Dx);      % Calculate convective flux
    [FluxConX, FluxConY]=calcFluxConUSCDS(RhoPhi,U,V,Dx,det); % Calculate convective flux
    [FluxDifX, FluxDifY]=calcFluxDif(RhoPhi,Dx,Mu/(Rho*Sc));  % Calculate diffusive flux
    
      RhoPhi(Ifi:Ila,Jfi:Jla) = RhoPhi(Ifi:Ila,Jfi:Jla) ...     % Apply flux
                                                 - det/(Dx^3)*( FluxDifX(Ifim:Ilam,Jfi:Jla)-FluxDifX(Ifi:Ila,Jfi:Jla)...
                                     + FluxDifY(Ifi:Ila,Jfim:Jlam)-FluxDifY(Ifi:Ila,Jfi:Jla));
    
    %% Scalar RhoO2
    %[FluxConX, FluxConY]=calcFluxConCDS(RhoPhi,U,V,Dx);      % Calculate convective flux
    [FluxConX, FluxConY]=calcFluxConUSCDS(RhoO2,U,V,Dx,det); % Calculate convective flux
    [FluxDifX, FluxDifY]=calcFluxDif(RhoO2,Dx,Mu/(Rho*Sc));  % Calculate diffusive flux
                                     

           RhoO2(Ifi:Ila,Jfi:Jla) = RhoO2(Ifi:Ila,Jfi:Jla) ...     % Apply flux
                 - det/(Dx^3)*( FluxDifX(Ifim:Ilam,Jfi:Jla)-FluxDifX(Ifi:Ila,Jfi:Jla)...
     + FluxDifY(Ifi:Ila,Jfim:Jlam)-FluxDifY(Ifi:Ila,Jfi:Jla));
    %% Scalar RhoN2
    %[FluxConX, FluxConY]=calcFluxConCDS(RhoPhi,U,V,Dx);      % Calculate convective flux
    [FluxConX, FluxConY]=calcFluxConUSCDS(RhoN2,U,V,Dx,det); % Calculate convective flux
    [FluxDifX, FluxDifY]=calcFluxDif(RhoN2,Dx,Mu/(Rho*Sc));  % Calculate diffusive flux
                                                                         
               RhoN2(Ifi:Ila,Jfi:Jla) = RhoN2(Ifi:Ila,Jfi:Jla) ...     % Apply flux
                 - det/(Dx^3)*( FluxDifX(Ifim:Ilam,Jfi:Jla)-FluxDifX(Ifi:Ila,Jfi:Jla)...
     + FluxDifY(Ifi:Ila,Jfim:Jlam)-FluxDifY(Ifi:Ila,Jfi:Jla));
    %% Scalar RhoFuel
    %[FluxConX, FluxConY]=calcFluxConCDS(RhoPhi,U,V,Dx);      % Calculate convective flux
    [FluxConX, FluxConY]=calcFluxConUSCDS(RhoFuel,U,V,Dx,det); % Calculate convective flux
    [FluxDifX, FluxDifY]=calcFluxDif(RhoFuel,Dx,Mu/(Rho*Sc));  % Calculate diffusive flux
                                                                         
                   RhoFuel(Ifi:Ila,Jfi:Jla) = RhoFuel(Ifi:Ila,Jfi:Jla) ...     % Apply flux
                 - det/(Dx^3)*( FluxDifX(Ifim:Ilam,Jfi:Jla)-FluxDifX(Ifi:Ila,Jfi:Jla)...
     + FluxDifY(Ifi:Ila,Jfim:Jlam)-FluxDifY(Ifi:Ila,Jfi:Jla));
                                     
    RhoFuel(1:x0,38:43) = Rho*0;
    RhoO2(1:x0,38:43) = Rho*0;
    RhoN2(1:x0,38:43) = Rho*0;
    RhoPhi(1:x0,38:43) = Rho*400;
    %% Correcting convective errors at the surface of the droplet                                         
    for i=Ifim:Ilap 
        for j=Jfim:Jlap 
            if (RhoPhi(i,j) < Rho*400) 
                RhoPhi(i,j) = 700*Rho; 
            end
            if (RhoFuel(i,j) < 0) 
                RhoFuel(i,j) = 0*Rho; 
            end
            if (RhoO2(i,j) < 0) 
                RhoO2(i,j) = 0*Rho; 
            end
            if (RhoN2(i,j) < 0) 
                RhoO2(i,j) = 0*Rho; 
            end
            if (RhoU(i,j) < 0) 
                RhoU(i,j) = 0*Rho; 
            end
        end
    end
        
    %% Evaporation Model by Spalding and Godsave
    
    [m0,r0,RhoO2,RhoN2,RhoFuel]=evapModel(Jma,m0,d0,r0,rho0,RhoO2,RhoN2,RhoFuel,det,B,Rho,Dx,Re,Sc,Mu);
    
    %% Scalar RhoFuel Re-evaluation
    %[FluxConX, FluxConY]=calcFluxConCDS(RhoPhi,U,V,Dx);      % Calculate convective flux
    [FluxConX, FluxConY]=calcFluxConUSCDS(RhoFuel,U,V,Dx,det); % Calculate convective flux
    [FluxDifX, FluxDifY]=calcFluxDif(RhoFuel,Dx,Mu/(Rho*Sc));  % Calculate diffusive flux
                                                                         
                   RhoFuel(Ifi:Ila,Jfi:Jla) = RhoFuel(Ifi:Ila,Jfi:Jla) ...     % Apply flux
                 - det/(Dx^3)*( FluxDifX(Ifim:Ilam,Jfi:Jla)-FluxDifX(Ifi:Ila,Jfi:Jla)...
     + FluxDifY(Ifi:Ila,Jfim:Jlam)-FluxDifY(Ifi:Ila,Jfi:Jla));  
    
    RhoFuel(1:x0,38:43) = Rho*0;
         
    %% Drop Combustion Model by Spalding 1953
    
    [RhoPhi,RhoFuel,RhoN2]=combModel(Ifi,Ila,Jfi,Jla,RhoPhi,RhoFuel,RhoN2,RhoO2,Jma,d0,x0,y0,B,cp,Rho,r0);
                                      
    %% Momentum RhoU in x-direction - without pressure
    %[FluxConX, FluxConY]=calcFluxConCDS(RhoU,U,V,Dx);        % Calculate convective flux
    [FluxConX, FluxConY]=calcFluxConUSCDS(RhoU,U,V,Dx,det);   % Calculate convective flux
    [FluxDifX, FluxDifY]=calcFluxDif(RhoU,Dx,Mu/Rho);         % Calculate diffusive flux
    RhoU(Ifi:Ila,Jfi:Jla) = RhoU(Ifi:Ila,Jfi:Jla) ...         % Apply flux
                          + det/(Dx^3)*( FluxConX(Ifim:Ilam,Jfi:Jla)-FluxConX(Ifi:Ila,Jfi:Jla)...
                                       + FluxConY(Ifi:Ila,Jfim:Jlam)-FluxConY(Ifi:Ila,Jfi:Jla)) ...
                          - det/(Dx^3)*( FluxDifX(Ifim:Ilam,Jfi:Jla)-FluxDifX(Ifi:Ila,Jfi:Jla)...
                                       + FluxDifY(Ifi:Ila,Jfim:Jlam)-FluxDifY(Ifi:Ila,Jfi:Jla));
                                      
    %% Momentum RhoV in y-direction - without pressure
    %[FluxConX, FluxConY]=calcFluxConCDS(RhoV,U,V,Dx);        % Calculate convective flux
    [FluxConX, FluxConY]=calcFluxConUSCDS(RhoV,U,V,Dx,det);   % Calculate convective flux
    [FluxDifX, FluxDifY]=calcFluxDif(RhoV,Dx,Mu/Rho);         % Calcualte diffusive flux
    RhoV(Ifi:Ila,Jfi:Jla) = RhoV(Ifi:Ila,Jfi:Jla) ...         % Apply flux
                          + det/(Dx^3)*( FluxConX(Ifim:Ilam,Jfi:Jla)-FluxConX(Ifi:Ila,Jfi:Jla)...
                                       + FluxConY(Ifi:Ila,Jfim:Jlam)-FluxConY(Ifi:Ila,Jfi:Jla)) ...
                          - det/(Dx^3)*( FluxDifX(Ifim:Ilam,Jfi:Jla)-FluxDifX(Ifi:Ila,Jfi:Jla)...
                                       + FluxDifY(Ifi:Ila,Jfim:Jlam)-FluxDifY(Ifi:Ila,Jfi:Jla));
                                   
    %% Re-evaluation of RhoU and RhoV                               
    for i=21:41
        for j=(Jma/2+1-d0*1000/2):(Jma/2+1+d0*1000/2)
            rij = floor(  sqrt(  (i-x0)^2 + (j-y0)^2  )  );
            if (rij <= (r0-1))
                RhoU(i,j) = 0;
                RhoV(i,j) = 0;
            end
        end
    end     
    RhoU(1:x0,38:43) = Rho*0;
    RhoV(1:x0,38:43) = Rho*0;
                                      
    %% 2. Boundary treatment
    if (Case == 1) % Free jet
      RhoPhi(Ilap,:) = RhoPhi(Ila,:);
      RhoPhi(:,Jfim) = RhoPhi(:,Jfi);
      RhoPhi(:,Jlap) = RhoPhi(:,Jla);
      RhoU(Ilap,:) = max(0,RhoU(Ila,:)); % Inhibit inflow
      RhoU(:,Jfim) = RhoU(:,Jfi);
      RhoU(:,Jlap) = RhoU(:,Jla);
      RhoV(Ilap,:) = RhoV(Ila,:);
      RhoV(:,Jfim) = min(0,RhoV(:,Jfi)); % Inhibit inflow
      RhoV(:,Jlap) = max(0,RhoV(:,Jla)); % Inhibit inflow
      %turb = 1.5;   % org = 0.5
      turb = 0.5;
      RhoV(Ifim,:) = turb*(rand(1,Jmap2)-0.5);     % Perturbations at inlet to trigger pseudo turbulence
      
      % More pertubations at inlet
%       Mag = 1.0;
%       freq = 0.03;
%       Mag = Mag*(sin(freq*(T/det)*pi)+1)/2;
%       RhoU(1:Imap2,1:Jmap2) = 0.0;
%       RhoU(1:Imap2,Jmidm-floor(JJet/2)-1:Jmidp+ceil(JJet/2)+1) = Rho*Mag/2;
%       RhoU(1:Imap2,Jmidm-floor(JJet/2):Jmidp+ceil(JJet/2)) = Rho*Mag;
      
    elseif (Case == 2) % Potential vortex
      RhoU(Ifim,:) = 0.0; RhoU(Ilap,:) = 0.0; RhoU(:,Jfim) = 0.0; RhoU(:,Jlap) = 0.0;
      RhoV(Ifim,:) = 0.0; RhoV(Ilap,:) = 0.0; RhoV(:,Jfim) = 0.0; RhoV(:,Jlap) = 0.0;
    end
    
    %% 3. Update velocity from new momentum fields
    [U,V]=mom2vel(RhoU,RhoV,Rho);
    
    %% 4. Calculate effect of pressure such, that div(U,V) = 0.0
    al=1.0;
    % Calculate divergence
    [div]=calcDivergence(U,V,Dx,Imap2,Jmap2);
    % Calculate the potential field PHI
    [PHI]=PoissonSolver(div,Imap2,Jmap2,Dx);
    % Adjust velocities to consider pressure forces
    U(1:Imap2-1,1:Jmap2)=U(1:Imap2-1,1:Jmap2) + al*(PHI(2:Imap2,1:Jmap2)-PHI(1:Imap2-1,1:Jmap2));
    V(1:Imap2,1:Jmap2-1)=V(1:Imap2,1:Jmap2-1) + al*(PHI(1:Imap2,2:Jmap2)-PHI(1:Imap2,1:Jmap2-1));
    RhoU(2:Imap2-1,2:Jmap2-1)=RhoU(2:Imap2-1,2:Jmap2-1) + al*Rho/2*(PHI(3:Imap2,2:Jmap2-1)-PHI(1:Imap2-2,2:Jmap2-1));
    RhoV(2:Imap2-1,2:Jmap2-1)=RhoV(2:Imap2-1,2:Jmap2-1) + al*Rho/2*(PHI(2:Imap2-1,3:Jmap2)-PHI(2:Imap2-1,1:Jmap2-2));
    
    %RhoUmag(Ifim:Imap2,Jfim:Jmap2) = sqrt( RhoU(Ifim:Imap2,Jfim:Jmap2).^2 + RhoV(Ifim:Imap2,Jfim:Jmap2).^2  );

    %% P. Postprocessing: output during calculation, only every 100 timesteps to save time
    if (mod(n,100) == 0); pause(1.0); 
        %contourf(RhoPhi/Rho); %, set(gca,'Ydir','reverse');
        imagesc(RhoPhi'/Rho); 
        %imagesc(RhoU'/Rho);
        %imagesc(RhoO2'/Rho); 
        %imagesc(RhoFuel'/Rho);
        colormap(hot);
        c = colorbar;
        ylabel(c,'Kelvin','FontWeight','bold','FontSize',12);
        
        xlabel('x-position (mm)','FontWeight','bold','FontSize',12);
        ylabel('y-position (mm)','FontWeight','bold','FontSize',12);
        title({['Time = ',num2str(T),' s'];
            %['Time-step = ',num2str(cont1)]
            %['Re = ',num2str(Re)]
            ['r = ',num2str(r0),' mm']
            ['m = ',num2str(m0*1000),' g']
            },'FontWeight','bold','FontSize',12);
        
        imgname = strcat('Burn',num2str(cont1));
        imgname = strcat(imgname,'_Re=');
        imgname = strcat(imgname,num2str(0));
        set(gcf,'PaperPositionMode','auto');
    
        print(FigHandle,imgname,'-dpng','-r0');
        %print(FigHandle,imgname,'-dpng','..\CFD_burningDrop\diffOnly\');
        %print(FigHandle,imgname,'-depsc','-r0');
    end
end
disp('Time stepping finished');


%% Post-processing
%  The imagesc routine sets the origin of the coordinate system at the
%  top-left corner, x pointing downwards, y to the right.
%imagesc(RhoU); pause();
%imagesc(RhoV); pause();
%imagesc(RhoPhi'/Rho); pause();