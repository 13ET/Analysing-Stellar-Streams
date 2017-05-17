function dwdt=dwHern(w, G, Mhalo, ahalo)
    
    %   Hernquist (halo) Potential:  phi=-G*Mhalo/(r+ahalo)
    %   (Specific) Acceleration from the halo:

    dxdt_2=-G*Mhalo*w(1)/(((w(1)^2+w(3)^2+w(5)^2)^0.5+ahalo)^2)*(w(1)^2+w(3)^2+w(5)^2)^(-0.5);
    
    dydt_2=-G*Mhalo*w(3)/(((w(1)^2+w(3)^2+w(5)^2)^0.5+ahalo)^2)*(w(1)^2+w(3)^2+w(5)^2)^(-0.5);
    
    dzdt_2=-G*Mhalo*w(5)/(((w(1)^2+w(3)^2+w(5)^2)^0.5+ahalo)^2)*(w(1)^2+w(3)^2+w(5)^2)^(-0.5);
   
    
    dwdt=[dxdt_2; dydt_2; dzdt_2];
end