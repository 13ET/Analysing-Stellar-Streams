function probability=loglikelihood_of_orbit()
% clear all;

%%%%% makes the b and r values as a function of l and then determines the
%%%%% log-likelihood for that orbit to the stream

%%%%   time
dt=0.5;            %timestep [dt*0.9778 Myr]
t=50;

load('observed_stream');   %loads collection of all start values
k=0;
while k~=1      %loop used to change the time t according to the need of the stream used (prevents overlaps in l values from the model)
    clearvars U Ugal Wendgal t_span t_span2 u1 u2 bmodel lmodel Rmodel d2 b l
    
    t_span=0:dt:t;    %time [0.9778 Myr]  
    t_span2=0:dt:t/2;    %time [0.9778 Myr]

    %%%%   ode45 integrator for particle
    U=zeros(length(t_span), 6);
    [~,u1]=ode45(@dw, t_span2, [Wend(1,1), Wend(1,2), Wend(1,3), Wend(1,4), Wend(1,5), Wend(1,6)]);           %ode45 integrates for each t in t_span, using start values Wend
    [~,u2]=ode45(@dw, t_span2, [Wend(1,1), -Wend(1,2), Wend(1,3), -Wend(1,4), Wend(1,5), -Wend(1,6)]); 
    U(1:length(t_span2), :)=flip(u1(:, :));     %W(particle number(n), time, x/vx/y/vy/z/vz) = w(time, x/vx/y/vy/z/vz) - flip needed
    U(length(t_span2):(length(t_span)), :)=u2(:, :);

    %%%%   galactic coordinates
    Ugal(:,:)=zeros(length(t_span), 5);         %allocates memory
    Wendgal(:,:)=zeros(length(Wend(:,1)), 5);

    for i=1:length(t_span)                      %values for orbit of centre particle
        Ugal(i, :)=convert_galactic(U(i, 1), U(i, 2), U(i, 3), U(i, 4), U(i, 5), U(i, 6));
    end
    for i=1:length(Wend(:, 1))                  %observed stream
        Wendgal(i, :)=convert_galactic(Wend(i, 1), Wend(i, 2), Wend(i, 3), Wend(i, 4), Wend(i, 5), Wend(i, 6));
    end

    %%%% sorts the stream (galactic values) with respect to the galactic longitude
    [~, d1]=sort(Wendgal(:,1));     %d1 is the order they're placed
    Wendgal(:, :)=Wendgal(d1, :);   
    d3=d1;
    if (Wendgal(end,1)-Wendgal(1,1))>350    %this is to not make the stream split at the 360-0 transition
        Wendgal(1,1)=Wendgal(1,1)+360;
        for i=2:length(Wendgal(:,1))-1         
            if abs(Wendgal(i-1, 1)-Wendgal(i, 1))>=350    %this approach doesn't work if stream is bigger than 350°
                Wendgal(i, 1)=Wendgal(i, 1)+360; 
            else
                i=length(Wendgal(:,1)); %ends loop  if there are no more values
            end                                
        end
    end
    
    %%%% interpolation              
    [~, d2]=sort(Ugal(:,1));        %sorts with respect to the galactic longitude
    Ugal(:, :)=Ugal(d2, :);        
    
    if (Ugal(end,1)-Ugal(1,1))>350    %%%   this is to not make the stream split at the 360-0 transition
        Ugal(1,1)=Ugal(1,1)+360;
        for i=2:length(Ugal(:,1))-1         
            if abs(Ugal(i-1, 1)-Ugal(i, 1))>=350    %this approach doesn't work if stream is bigger than 350°
                Ugal(i, 1)=Ugal(i, 1)+360; 
            else
                i=length(Ugal(:,1)); %ends loop  if there are no more values
            end                                
        end
    end
    
    [~, d2]=sort(Ugal(:,1));        %sorts with respect to the galactic longitude, so that model is continuous
    Ugal(:, :)=Ugal(d2, :);
    %treats another possibility for error: if min of model is under 0°
    if min(Ugal(:,1))>min(Wendgal(:,1)) && max(Ugal(:,1))>360 && abs(max(Ugal(:,1))-max(Wendgal(:,1)))>180      
       Ugal(:,1)=Ugal(:,1)-360;                                     %then above code +360 even though stream is only close to 0
    end
    
    l=Wendgal(:,1);     %renaming for the sake of simplicity
    b=Wendgal(:,2);
    
    i=1;
    for j=1:length(Ugal(:,1))   %this sets model values in the same longitude range as the observed stream
        if Ugal(j,1)>=min(l)-5 && Ugal(j,1)<=max(l)+5           %within 5°, another safety precaution       
            lmodel(i)=Ugal(j, 1);
            bmodel(i)=Ugal(j, 2);
            Rmodel(i)=((U(j, 1)-8000)^2+U(j, 3)^2+U(j, 5)^2)^0.5;
            if (max(lmodel)<max(l) && max(Ugal(:,1))>max(l))        
                jend=j;                                             
            end
            if i==1 &&( (min(lmodel)>min(l) && min(Ugal(:,1))<min(l)) )
               lmodel(2)=lmodel(1);                 %these two if statements make sure values are found in the stream range
               lmodel(1)=Ugal(j-1, 1);              %by extending to one extra value before and after the value closest to 
               bmodel(2)=bmodel(1);                 %the stream l has been found
               bmodel(1)=Ugal(j-1, 2);              
               Rmodel(2)=Rmodel(1);
               Rmodel(1)=((U(j-1, 1)-8000)^2+U(j-1, 3)^2+U(j-1, 5)^2)^0.5;
               i=2;
            end
            i=i+1;
        end
    end
    if (max(lmodel)<max(l) && max(Ugal(:,1))>max(l))       %this is also meant to ensure l values in the range, adds the last value
        lmodel(i)=Ugal(jend+1, 1);
        bmodel(i)=Ugal(jend+1, 2);
        Rmodel(i)=((U(jend+1, 1)-8000)^2+U(jend+1, 3)^2+U(jend+1, 5)^2)^0.5;
    end
    
    k=1;
    if max(lmodel)<max(l) || min(lmodel)>min(l) %if we dont have values in the range of the stream,
        k=0;                                    %retry with longer integration time
        t=t+50;                                 %observed_stream4 shows a case where we never get l values corresponding to these 
    end
    if t>1000           %safety precaution, so there isn't an infinite loop in case values can't be found
       probability=-inf;
       return;
    end
end
%d3 is the order of the galactic longitudes for Wendgal, sort Wend with this 
%- problem before when trying to sort already sorted array
 Wend(:,:)=Wend(d3, :);  

r=((Wend(:, 1)-8000).^2+Wend(:, 3).^2+Wend(:, 5).^2).^0.5;     %radius for all the stream objects
ltot(1:length(lmodel))=lmodel(:);                       % I take all the lmodel values and add the longitude values of the stream lstream,
ltot((length(lmodel)+1):(length(lmodel)+length(l)))=l;  % this is done so interpolation can be done for the these longitudes
ltot(:)=sort(ltot(:));       %sorting for continuous

bmodeltot=zeros(length(ltot), 1);   %allocating memory
Rmodeltot=zeros(length(ltot), 1);
for j=1:(length(ltot))
    bmodeltot(j)=interp1(lmodel(:), bmodel(:), ltot(j)); %the total bmodel, when using l for bmodel and interpolating the values for
    Rmodeltot(j)=interp1(lmodel(:), Rmodel(:), ltot(j)); %stream longitudes Wendgal(:,1)
end

bmodelobjects=zeros(length(l), 1);    %allocating memory
Rmodelobjects=zeros(length(l), 1);
i=1;
for j=1:length(ltot)
    if i<=length(l)
        if ltot(j)==l(i)
           bmodelobjects(i)=bmodeltot(j);   %bmodelstream is the corresponding b value for the model at the l values of the stream
           Rmodelobjects(i)=Rmodeltot(j);
           i=i+1;
        end
    end
end
%%%% probability calculations
omegab=1;           %1° uncertainty
omegar=3*10^3;      %3kpc uncertainty
P=log(1);
for i=(1):(length(l))    
    P=P + log(1/(2*pi*omegab^2)^0.5*exp(-0.5*(b(i)-bmodelobjects(i))^2/omegab^2))...
        + log(1/(2*pi*omegar^2)^0.5*exp(-0.5*(r(i)-Rmodelobjects(i))^2/omegar^2));
end

%%%% confirmation that interpolation etc works
% figure();
% plot(lmodel, bmodel, '.');
% hold on;
% plot(l, b, '*');
% set(gca, 'FontSize', 16);
% ylabel('Galactic latitude, b [\circ]', 'FontSize', 18);
% xlabel('Galactic longitude, $\ell$ [$^\circ$]', 'interpreter', 'latex', 'FontSize', 18);
% if min(l)>15
%     set(gca,'XTick',15:30:720);                  %tick labels, just for accurate presentation when above 360
%     set(gca,'XTickLabel', (0:15:345)+15);        %(it's inaccurate betweeen 0-15) 
% end
% figure();
% plot(lmodel, Rmodel/1000, '.');
% hold on;
% plot(l, r/1000, '*');
% set(gca, 'FontSize', 16);
% ylabel('Distance from the Sun, r [kpc]', 'FontSize', 18);
% xlabel('Galactic longitude, $\ell$ [$^\circ$]', 'interpreter', 'latex', 'FontSize', 18);
% if min(l)>15
%     set(gca,'XTick',15:30:720);                  %tick labels, just for accurate presentation when above 360
%     set(gca,'XTickLabel', (0:15:345)+15);        %(it's inaccurate betweeen 0-15)
% end
% max(Wendgal(:,1))


probability=P;

end