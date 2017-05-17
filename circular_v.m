function circular_orbit=circular_v(r0)
load('constants'); 

%            the equation for the circular velocity
vc=abs((G.*Mdisk./(r0.^2+(adisk+bdisk)^2).^1.5.*r0.^2+(G*Mhalo.*r0)./(ahalo+r0).^2).^0.5);

circular_orbit=vc;

end