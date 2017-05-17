clear all;  

%%%%% creates the confidence region by determining the log-likelihood of
%%%%% the orbit fit for a range of masses

%%%%   starting values
constants0();   %resets constants to "true" values, just in case
load('constants');  %saves constants
Mdisk0=Mdisk;
Mhalo0=Mhalo;
%%%%    range of mass values (from zero to the true values + this range)
Mhalorange=4*10^10;
Mdiskrange=6*10^10;
% both disc and halo use the same mass-step
dM=2*10^10;

j=1;
for Mdisk=(-0*10^10):dM:(Mdisk0+Mdiskrange)
    i=1;
    for Mhalo=(-0*10^10):dM:(Mhalo0+Mhalorange)
        save('constants', 'G', 'Mdisk', 'Mhalo', 'adisk', 'bdisk', 'ahalo');
        Ptot(i, j)=loglikelihood_of_orbit();
        Mtest(i, j, :)=[Mdisk Mhalo];   %necessary to determine the masses for max loglikelihood 
        i=i+1;
    end
    j=j+1;
end


%%%%            finds the max log-likelihood value
Ptotmax=-inf;
for i=1:length(Ptot(:,1))  
    for j=1:length(Ptot(1,:))
        if Ptot(i, j)>Ptotmax
          Ptotmax=Ptot(i,j);
          Ptotmaxij=[i j];
        end
    end
end
Ptotmax     % shows max log-likelihood value
Mtest(Ptotmaxij(1), Ptotmaxij(2), :)    %show max mass

%%%%    "calculate" the area of the 90% CR
A=0;            
for i=1:length(Ptot(:,1))
    for j=1:length(Ptot(1,:))
        if Ptot(i, j)>Ptotmax-2.3
          A=A+(dM/10^10)^2;     %area is dM^2
        end
    end
end
A  %show area

% plots CR
figure();
contourf(Ptot, (Ptotmax-2.3):0.5:(Ptotmax));   %2.3 is the confidence region with 90% confidence 

hold on;
scatter(Ptotmaxij(2), Ptotmaxij(1));            %plots maximal likelihood
hold on;
scatter(8/dM*10^10+1, 20/dM*10^10+1, 'x', 'linewidth', 1.5); %plots ``true'' value (works as long as dM is divisible by 8 & 20

% labels
colorbar();
xlabel('$M_{\mathrm{disc}} \hspace{1mm} [10^{10} M_{\odot}]$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$M_{\mathrm{halo}} \hspace{1mm} [10^{10} M_{\odot}]$', 'Interpreter', 'latex', 'FontSize', 20);
title('log-likelihood', 'FontSize', 20);

set(gca,'XTick', (0:2*dM:(Mdiskrange+Mdisk0+dM))/dM);                
set(gca,'XTickLabel', (0:2*dM:(Mdiskrange+Mdisk0+dM))/10^10-dM/10^10, 'FontSize', 16);
set(gca,'YTick', (0:2*dM:(Mhalorange+Mhalo0+dM))/dM);                
set(gca,'YTickLabel', (0:2*dM:(Mhalorange+Mhalo0+dM))/10^10-dM/10^10, 'FontSize', 16);


%%%% resets values to "true" potential
constants0();




