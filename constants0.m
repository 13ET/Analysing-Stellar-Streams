function constants0()   % I use this file to change the potential,
                        % edit numbers then run and it's saved to a file

G=0.004301;         %[(km/s)^2*pc*Msun^-1]

Mdisk=8*10^10;   %mass of MN disk [Msun]              

Mhalo=20*10^10;     %mass of Hernquist halo [Msun]   

adisk=3.7*10^3;     %scale length [pc]   

bdisk=0.2*10^3;     %scale height [pc]  

ahalo=18.5*10^3;      %scale length [pc] 


save('constants', 'G', 'Mdisk', 'Mhalo', 'adisk', 'bdisk', 'ahalo')
end