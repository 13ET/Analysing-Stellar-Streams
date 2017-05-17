function circular_orbit=circular_v(r0)
load('constants'); %loads G, Mdisk, Mhalo, adisk, bdisk, ahalo

%   for testing that orbits are circular: vc^2/r=-dphi/dr

%                 for disk
%vy0=(G*Mdisk/(x0^2+(adisk+bdisk)^2)^1.5*x0^2)^0.5;

%                 for disk + halo   
vc=abs((G.*Mdisk./(r0.^2+(adisk+bdisk)^2).^1.5.*r0.^2+(G*Mhalo.*r0)./(ahalo+r0).^2).^0.5);
   
%velocity when dphi(x, y,z)/dr=dphi/dx dx/dr + dphi/dy dy/dr + dphi/dz dz/dr - which is wrong     
% vc=(G*Mdisk/(r0^2+(adisk+bdisk)^2)^1.5*r0^2*(2+adisk/bdisk))^0.5;

circular_orbit=vc;

end