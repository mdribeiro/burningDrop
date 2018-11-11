%% Pressure Example M
%  Program for Lecture 4 of Numerics and Flow Simulation
%  Andreas Kempf, 16.05.2013

%% Initialisation
clc; clear; hold off;

%% Create arrays X and Y
[X,Y]=meshgrid(-5.25:0.5:5.25,-5.25:0.5:5.25);

%% Calculate pressure field P
%P=exp(-X.^2-Y.^2);
%P=exp(-(X+1).^2-(Y+1).^2) - exp(-(X-1).^2-(Y-1).^2);
P=1/(2*pi)*log(sqrt(X.^2+Y.^2));

%% Calculate velocities, assume Dt/(2 rho Delta) = 1
% Initialisation
Ima=size(X,1); Jma=size(X,2);
U=zeros(Ima,Jma); V=zeros(Ima,Jma);
% Calculation
U(:,2:Jma-1)=(P(:,1:Jma-2)-P(:,3:Jma));
V(2:Ima-1,:)=(P(1:Ima-2,:)-P(3:Ima,:));

%% Calculate divergence
% Initialisation
DivUV=zeros(Ima,Jma);
% Calculation
DivUV(2:Ima-1,2:Jma-1) = (U(2:Ima-1,3:Jma)-U(2:Ima-1,1:Jma-2))...
                       + (V(3:Ima,2:Jma-1)-V(1:Ima-2,2:Jma-1));

%% Visualisation
surface(X,Y,P); pause();
quiver(X,Y,U,V); pause();
surface(X,Y,DivUV);