function dwdt=dwMN(w, G, Mdisk, a, b)
    %   Miyamoto-Nagai (disc) Potential:   phi=-G*Mdisk/(R^2+(a+(z^2+b^2)^0.5)^2)^0.5
    %   (Specific) Acceleration from the disc:

    dxdt_2=-G*Mdisk*(w(1)/(w(1)^2+w(3)^2)^0.5)*((w(1)^2+w(3)^2)^0.5)/(w(1)^2+w(3)^2+(a+(w(5)^2+b^2)^0.5)^2)^1.5;
    
    dydt_2=-G*Mdisk*(w(3)/(w(1)^2+w(3)^2)^0.5)*((w(1)^2+w(3)^2)^0.5)/(w(1)^2+w(3)^2+(a+(w(5)^2+b^2)^0.5)^2)^1.5;
    
    dzdt_2=-(G*Mdisk/((w(1)^2+w(3)^2)+(a+(w(5)^2+b^2)^0.5)^2)^1.5)*(w(5)*(a+(w(5)^2+b^2)^0.5)*(w(5)^2+b^2)^(-0.5));
   
    
    dwdt=[dxdt_2; dydt_2; dzdt_2];
end