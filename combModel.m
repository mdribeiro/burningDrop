%% Evaporation model

function [RhoPhi,RhoFuel,RhoN2]=combModel(Ifi,Ila,Jfi,Jla,RhoPhi,RhoFuel,RhoN2,RhoO2,Jma,d0,x0,y0,B,cp,Rho,r0)

    Ts = 400*Rho;
    
    Q = 280;        % [kJ/kg] Latent Heat of Iso-octane
    H = -47900;     % [kJ/kg] Iso-octane Heat of Combustion
    %mO2 = RhoO2/Rho;    % [-] Mixture Fraction of O2
    r = 3.51;       % [kg/kg] Weight of Oxygen Required for Combustion of Unit Weight of Fuel
        
    for i=21:41
        for j=(Jma/2+1-d0*1000/2):(Jma/2+1+d0*1000/2)
            rij = floor(  sqrt(  (i-x0)^2 + (j-y0)^2  )  );
            if (rij <= (r0-1))
                RhoPhi(i,j) = Ts;
            end

        end
    end  

    resteq = 1/14;
    %% Complete flame zone combustion
    for i=Ifi:Ila
        for j=Jfi:Jla
            
            rel = RhoFuel(i,j)/RhoO2(i,j);
            
            if (rel > resteq)
                RhoPhi(i,j) = Rho*(Ts + Q*B/cp - (H/cp)*((RhoO2(i,j)/Rho)/r));

                RhoFuel(i,j) = 0;   % All fuel is consumed
                RhoO2(i,j) = RhoO2(i,j) - RhoO2(i,j)*resteq;  
                RhoN2(i,j) = 1-RhoO2(i,j)-RhoFuel(i,j);                
            end
        end
    end      


    %% Flame extintion at the base of the channel 
    for i=21:22
       for j=(y0-3):(y0+3)
            RhoPhi(i,j) = Ts;
       end
    end
