clear all;
% close all;

%%%%% used to create the streams by integrating test particles

%%%%   time
dt=0.1;            %timestep [dt*0.9778 Myr]
t_span=0:dt:4500;    %time [0.9778 Myr]
%%%%   resets constants, just in case
constants0();       

%%%%   variable initial conditions 

n=100;               %number of test particles
r0=12*10^3;         %start position for the central particle of the satellite [pc]
rrange=0.1*10^3;    %the "radius" of the satellite, which the particles lie in
vc=circular_v(r0);        %circular velocity
v0=1.2*vc;     %0.9 - 1.3              %the initial speed for the centre particle
vrange=0.05*vc;                             %range of velocity
theta=45;                   %velocity's angle from disc plane, in degrees

dispersion=randn(1,n);     %randn gives a normal distribution between neg. pos. infty - unknown highest values         

%%%%   starting position from GC [pc]
x0=r0+rrange*(dispersion/max(abs(dispersion)));    %randn unpredictable max and min values, dividing with max gives values ~[-1, 1]
x0(1)=r0;                                          %this sets one particle at the actual centre
dispersion=randn(1,n);
y0=rrange*(dispersion/max(abs(dispersion)));      
y0(1)=0;                                          %this sets one particle at the actual centre
dispersion=randn(1,n);
z0=rrange*(dispersion/max(abs(dispersion)));
z0(1)=0;                                          %this sets one particle at the actual centre

dispersion=randn(1,n);              %new random numbers, so position & velocity dispersions are different

%%%%   starting velocities in x, y, z  [km/s]
vx00=0;                          %start-velocity in x
vy00=cos(theta*pi/180)*v0;          %start-velocity in y, angle made into radians
vz00=sin(theta*pi/180)*v0;          %start velocity in z, angle made into radians

vx0(1:n)=vx00+vrange*(dispersion/max(abs(dispersion)))/3^0.5;
vx0(1)=vx00;                     %particle at centre has no velocity dispersion
dispersion=randn(1,n);
vy0(1:n)=vy00+vrange*(dispersion/max(abs(dispersion)))/3^0.5;
vy0(1)=vy00;                     %particle at centre has no velocity dispersion
dispersion=randn(1,n);
vz0(1:n)=vz00+vrange*(dispersion/max(abs(dispersion)))/3^0.5;
vz0(1)=vz00;                     %particle at centre has no velocity dispersion

%%%%%   ode45 integrator for each particle 1:n

W(n, length(t_span), 6)=zeros();  %allocating memory for W, used below
for i=1:n
    w0=[x0(i)   vx0(i)   y0(i)    vy0(i)  z0(i)   vz0(i)];   %collection of all start values
    [t,w]=ode45(@dw, t_span, w0);           %ode45 integrates for each t in t_span, using start values w0
    W(i, :, :)=w(:, :);     %W(particle number(n), time, x/vx/y/vy/z/vz) = w(time, x/vx/y/vy/z/vz)
end



%%%%%   plots the integrated orbits

figure();
for i=1:n
    %plot(W(i, :,1)/1000, W(i, :, 3)/1000);                           %  2D   -  x&y    
    plot3(W(i, :,1)/1000, W(i, :, 3)/1000, W(i, :, 5)/1000);          %  3D   -  x&y&z
    hold on;
end

ylabel('y-value [kpc]');
xlabel('x-value [kpc]');
zlabel('z-value [kpc]');

%%%%  plots the values of the stream at some time 2D/3D

Wend(:,:)=zeros(n, 6);
figure();
for i=1:n
    Wend(i, :)=[W(i, length(t_span)-2000, 1) W(i, length(t_span)-2000, 2) W(i,length(t_span)-2000, 3) ...
    W(i, length(t_span)-2000, 4) W(i,length(t_span), 5) W(i, length(t_span)-2000, 6)];
    %plot(Wend(i, 1)/1000, Wend(i, 3)/1000, '*');                      %2D
    plot3(Wend(i, 1)/1000, Wend(i, 3)/1000, Wend(i, 5)/1000, '*');   %3D
    hold on;
end
% plot(Wend(1, 1)/1000, Wend(1, 3)/1000, '+');                      %2D
plot3(Wend(1, 1)/1000, Wend(1, 3)/1000, Wend(1, 5)/1000, '+');       %3D
hold on;
% plot(0,0, 'o');                                               %2D
plot3(0,0,0, 'o');                                              %3D
hold on;
% plot(W(1, :, 1)/1000, W(1, :, 3)/1000);                      %2D
plot3(W(1, :, 1)/1000, W(1, :, 3)/1000,  W(1, :, 5)/1000);             %3D
ylabel('y-value [kpc]');
xlabel('x-value [kpc]');
zlabel('z-value [kpc]');                                               %3D
    

save('observed_stream', 'Wend');        