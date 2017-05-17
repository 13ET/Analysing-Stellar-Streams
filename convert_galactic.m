function convert=convert_galactic(x, vx, y, vy, z, vz)

xsun=-(x-8000);    %x is directed towards GC, not outwards, therefore '-'
ysun=(y-0);        %heliocentric galactic coordinates - at the solar system ~(8, 0, 0)kpc
zsun=(z-0);        %the angles l & b are based from the solar system

vxsun=-(vx-0);    %x is directed towards GC, not outwards, therefore '-'
vysun=(vy-220);
vzsun=(vz-0);

R=(xsun^2+ysun^2)^0.5;
r=(xsun^2+ysun^2+zsun^2)^0.5;


% b value, galactic latitude
if R==0
    if zsun>0           %have to define what happens when R=0 - matlab can't handle zsun/0 (=NaN)
        b=pi/2;         %arctan(x) as x->infty = pi/2, x=zsun/R
    else                
        b=-pi/2;
    end
else
    b=atan(zsun/R);
end

% l value, galactic longitude - 4 quadrants & "boundary conditions"

if xsun~=0
    l=atan(abs(ysun/xsun));
    if xsun<0 && ysun>0         % visualise quadrants by making a drawing
        l=pi-l;                 % and consider the arctan function
    end                         
    if xsun<0 && ysun<0
        l=pi+l;
    end
    if xsun>0 && ysun<0
        l=2*pi-l;
    end
end
if xsun==0                             %at xsun=0 - boundaries
    if ysun>0                    
        l=pi/2; 
    else
        l=3*pi/2;
    end
    if ysun==0
        l=0; 
    end
end
if ysun==0                             %at ysun=0 - boundaries
   if xsun>0
      l=pi; 
   end
   if xsun<0
       l=0;
   end
end

%proper motion
K=4.7405;              %factor to turn [mas/yr*kpc] to [km/s]   -  milliarcsec

vvec=[vxsun; vysun; vzsun];
uvec=[cos(b)*cos(l); cos(b)*sin(l); sin(b)];             %vectors (the normal triad) shown in astm13 lecture notes, pg10
lvec=[-sin(l); cos(l); 0];                               %and equations for proper motion pg 10-11, eq. 2.6/11    
bvec=[-sin(b)*cos(l); -sin(b)*sin(l); cos(b)]; 

mul=(lvec.')*vvec*(1000/(K*r));                        %p(arallax)[mas]=1000/r[pc]
mub=(bvec.')*vvec*(1000/(K*r));
vr=(uvec.')*vvec;

convert=[l*180/pi, b*180/pi, mul, mub, vr];
end