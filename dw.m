function dwdt=dw(t, w)
    load('constants'); %loads G, Mdisk, Mhalo, adisk, bdisk, ahalo
    
    %  total potential in all 3 directions, dwdt2(1) is for x, 2 for y, 3 for z
    dwdt2=dwMN(w, G, Mdisk, adisk, bdisk)+dwHern(w, G, Mhalo, ahalo);
    
    %  total potential in seperate directions dwdt2(x/y/z)
    x=w(2);
    dxdt1=dwdt2(1);
    
    y=w(4);
    dydt1=dwdt2(2);
    
    z=w(6);
    dzdt1=dwdt2(3);
    
    dwdt=[x; dxdt1; y; dydt1; z; dzdt1];
end