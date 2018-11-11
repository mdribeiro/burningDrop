%% Evaporation model

function [m0,r0,RhoO2,RhoN2,RhoFuel]=evapModel(Jma,m0,d0,r0,rho0,RhoO2,RhoN2,RhoFuel,det,B,Rho,Dx,Re,Sc,Mu)

    dA = 4*pi*((r0*0.001)^2);   %[m2]
    kc = (0.6*Re^0.5+Sc^0.333+2)*(Mu/(Rho*Sc))/(2*r0*0.001);
    mevap = kc*Rho*log(1+B);    %[kg/(s*m2)]
    
    % Recalculating droplet properties after evaporation
    mdif = m0 - mevap*det*dA;    %[kg]
    r0 = ( ( (3/4)*(mdif)/(pi*rho0)  )^(1/3)   )*1000;
    V0 = (4/3)*pi*(r0*0.001)^3;
    m0 = rho0*V0;
    
    mt = Rho*Dx^2;    % Talvez o mais correto fisicamente seria Rho = rho0      e mt = rho*Dx^3
    x0 = 31;
    y0 = Jma/2+1;
    %rfo2 = 10;
    %rfo2 = 1;  %Fazer esse numero extremamente baixo pra gerar uma condicao rica na superficie.... depois a Eq. de transporte se encarrega do resto
    for i=21:41
        for j=(Jma/2+1-d0*1000/2):(Jma/2+1+d0*1000/2)
            rij = floor(  sqrt(  (i-x0)^2 + (j-y0)^2  )  );
            if (rij <= (r0-1))
                RhoO2(i,j) = 0;
                RhoN2(i,j) = 0;
                RhoFuel(i,j) = 0;
            elseif (  (rij > (r0-1))&&(rij <= (r0+1))  )
                RhoFuel(i,j) = Rho*mevap*det*dA/mt;
                RhoO2(i,j) = RhoO2(i,j) - RhoO2(i,j)*RhoFuel(i,j)/Rho;
                RhoN2(i,j) = 1-RhoO2(i,j)-RhoFuel(i,j);
            end

        end
    end  