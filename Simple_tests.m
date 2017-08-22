filename = 'C:\Users\Bradley\Documents\UNIVERSITY\Part D\Total Product Design\Newest Matlab\Route B\test 8 appended.csv';
data_raw = csvread(filename);
data_calibrated = Calibrate04008(data_raw);
sampling_rate = 0.02; %50Hz
%Angles
%starting_theta = [0, -PI./2 ,0]
l = length(data_calibrated);
Theta = zeros(l,3);
Theta(1,1) = 0.034113388;
Theta(1,2) = -1.34362605;
Theta(1,3) = 0.007888563786;
for i= 2:l;
    for j= 1:3;
        Theta(i,j) = Theta(i-1,j)+data_calibrated(i-1,j+3)*sampling_rate;
    end;
end;

%TIME matrix
T = zeros(l,1);
for i = 1:l;
    T(i) = (i-1)*sampling_rate;
end;
%figure;
% plot(Theta);
XYZ_temp = zeros(3,1);
XYZ = zeros(l,3);
for i=1:l;
        xyz = (data_calibrated(i, 1:3)).';
        X_rot = [1,0,0;0,cos(Theta(i,1)),-sin(Theta(i,1));0,sin(Theta(i,1)),cos(Theta(i,1))];
        Y_rot = [cos(Theta(i,2)),0,sin(Theta(i,2));0,1,0;-sin(Theta(i,2)),0,cos(Theta(i,2))];
        Z_rot = [cos(Theta(i,3)),-sin(Theta(i,3)),0;sin(Theta(i,3)),cos(Theta(i,3)),0;0,0,1];
        XYZ_temp = Z_rot*Y_rot*X_rot*xyz;
        XYZ(i,1) = XYZ_temp(1,1);
        XYZ(i,2) = XYZ_temp(2,1);
        XYZ(i,3) = XYZ_temp(3,1);
end


Xacc = XYZ(:,1);
Yacc = XYZ(:,2);
Zacc = XYZ(:,3);

% figure
% plot(Xacc,'r-')
% figure
% plot(Yacc,'g-')
% figure
% plot(Zacc,'b-')

Xfreq = fft(abs(Xacc));
X2 = abs(Xfreq/l);
X1 = X2(1:l/2+1);
X1(2:end-1) = 2*X1(2:end-1);

Yfreq = fft(abs(Yacc));
Y2 = abs(Yfreq/l);
Y1 = Y2(1:l/2+1);
Y1(2:end-1) = 2*Y1(2:end-1);

Zfreq = fft(abs(Zacc));
Z2 = abs(Zfreq/l);
Z1 = Z2(1:l/2+1);
Z1(2:end-1) = 2*Z1(2:end-1);

f = (0:(l/2))/(l*sampling_rate);
% figure
% plot(f,X1,'r-',f,Y1,'g-',f,Z1,'b-')

Xvel_start = 0;
Xdist_start = 0;
Xvel = zeros(l,1);
Xdist = zeros(l,1);
Xvel(1,1) = Xvel_start;
Xdist(1,1) = Xdist_start;
for i=2:l;
    Xvel(i,1) = Xvel(i-1,1)+sampling_rate*(Xacc(i,1)+Xacc(i-1))/2;
    Xdist(i,1) = Xdist(i-1,1)+sampling_rate*(Xvel(i,1)+Xvel(i-1))/2;
end

Yvel_start = 0;
Ydist_start = 0;
Yvel = zeros(l,1);
Ydist = zeros(l,1);
Yvel(1,1) = Yvel_start;
Ydist(1,1) = Ydist_start;
for i=2:l;
    Yvel(i,1) = Yvel(i-1,1)+sampling_rate*(Yacc(i,1)+Yacc(i-1))/2;
    Ydist(i,1) = Ydist(i-1,1)+sampling_rate*(Yvel(i,1)+Yvel(i-1))/2;
end

Zvel_start = 0;
Zdist_start = 0;
Zvel = zeros(l,1);
Zdist = zeros(l,1);
Zvel(1,1) = Zvel_start;
Zdist(1,1) = Zdist_start;
for i=2:l;
    Zvel(i,1) = Zvel(i-1,1)+sampling_rate*(Zacc(i,1)+1+Zacc(i-1)+1)/2;
    Zdist(i,1) = Zdist(i-1,1)+sampling_rate*(Zvel(i,1)+Zvel(i-1))/2;
end

% figure
% plot(T,Xdist,'r-')
% figure
% plot(T,Ydist,'g-')
% figure
% plot(T,Zdist,'b-')
avg_length1 = 40; %samples must even
avg_Zgiro = zeros(l,10);
avg_Zgiro(:,1) = T(:,1);
for i = (0.5*avg_length1+1):(l-0.5*avg_length1);
    for j = (i-0.5*avg_length1):(i+0.5*avg_length1);
        avg_Zgiro(i,2) = avg_Zgiro(i,2) + (1/avg_length1)*Theta(j,1);
    end;
end;
peak_length = 150;
for i = (0.5*peak_length+1):(l-0.5*peak_length);
    avg_Zgiro(i,3)= -(avg_Zgiro((i-0.5*peak_length),2))+(avg_Zgiro((i+0.5*peak_length),2));
end;

for i = 1:l;
    if avg_Zgiro(i,3) > 0.6 || avg_Zgiro(i,3) < -0.6 ;
        avg_Zgiro(i,4) = avg_Zgiro(i,3);
    else;
        avg_Zgiro(i,4) = 0;
    end;
    
end;

peak_centre_avg_width = 30;
x=0;
for i = (0.5*peak_centre_avg_width+1):(l-0.5*peak_centre_avg_width);
    if abs(avg_Zgiro(i,4)) > 0;
    for j = (i-0.5*peak_centre_avg_width):(i+0.5*peak_centre_avg_width);
        avg_Zgiro(i,5) = avg_Zgiro(i,5)+ (1/(peak_centre_avg_width+1))*avg_Zgiro(j,4);
    end;
    end;
end;

for i = 1:l;
    if avg_Zgiro(i,5)>0 && avg_Zgiro(i,5)<avg_Zgiro(i-1,5) && avg_Zgiro(i-1,5)>avg_Zgiro(i-2,5);
        avg_Zgiro(i,6) = avg_Zgiro(i,4);
    else if avg_Zgiro(i,5)<0 && avg_Zgiro(i,5)>avg_Zgiro(i-1,5) && avg_Zgiro(i-1,5)<avg_Zgiro(i-2,5);
            avg_Zgiro(i,6) = avg_Zgiro(i,4);
        else;
            avg_Zgiro(i,6) = 0;
        end;
    end;
end;
a = 1;
while a<=l;
while abs(avg_Zgiro(a,6))>0;
    c=a;
        while a<=l && abs(avg_Zgiro(a,5))>0
            if abs(avg_Zgiro(a,6))>0;
            avg_Zgiro(a,7) = avg_Zgiro(a,6);
            else;
            avg_Zgiro(a,7) = 0;
            end;
            a=a+1;
        end;
        d  = max(abs(avg_Zgiro(:,7)));
        for b = c:a;
            if abs(avg_Zgiro(b,7))==d;
                avg_Zgiro(b,8)=avg_Zgiro(b,7);
            else;
                avg_Zgiro(b,8)=0;
            end;
            end;
            avg_Zgiro(:,7)=0;
end;
a=a+1;
end;

% f = 1
% while f < l;
%     while abs(avg_Zgiro(f,6))>0;
%         if abs(avg_Zgiro(f,6))>abs(avg_Zgiro((f-1),6));
%             avg_Zgiro(f,6) = avg_Zgiro(f,6);
%         else;
%             avg_Zgiro(f,6) = 0;
%         end;
%         f = f+1;
%     end;
%     f=f+1;
% end; 

for i = 1:l;
    h = avg_Zgiro(i,8);
    if h>0.6 && h<1;
        avg_Zgiro(i,9) = 0.7854;
    else if h<-0.6 && h>-1;
            avg_Zgiro(i,9) = -0.7854;
    else if h>=1 && h<2.5;
        avg_Zgiro(i,9) = 1.5708;
    else if h<=-1 && h>-2.5 ;
        avg_Zgiro(i,9) = -1.5708;
    else if h>=2.5 && h<4;
        avg_Zgiro(i,9) = 3.1416;
    else if h<=-2.5 && h>-4;
        avg_Zgiro(i,9) = -3.1416;
    else if h>=4 && h<5.5;
        avg_Zgiro(i,9) = 4.7124;
    else if h<=-4 && h>-5.5;
        avg_Zgiro(i,9) = -4.7124;
    else if h>=5.5 && h<7;
        avg_Zgiro(i,9) = 6.2832;
    else if h<=-5.5 &&h>-7;
        avg_Zgiro(i,9) = -6.2832;
        else;
            avg_Zgiro(i,9)=0;
        end;
        end;
        end;
        end;
        end;
        end;
        end;
        end;
        end;
    end;
end;
for i = 2:l;
    avg_Zgiro(i,10)= avg_Zgiro(i-1,10)+avg_Zgiro(i,9);
end;


e=1;
for i = 1:l;
    if abs(avg_Zgiro(i,9))>0;
        e = e+1;
    else;
        e=e;
    end;
end;

dir_vec = zeros(e,1);
f = 1;
for i = 1:l;
    if abs(avg_Zgiro(i,9))>0;
        dir_vec(f,1) = avg_Zgiro(i,9);
        f=f+1;
    else;
        f=f;
    end;
end;
dir_vec
    
aMech_DSuite_stairs = [1.5708;1.5708;-0.7854;0.7854;3.1416;-1.5708;1.5708;1.5708;3.1416;3.1416;1.5708;1.5708;-1.5708;0;0;0;0;0;0;0];
bDSuite_Mech_stairs = [1.5708;-1.5708;-1.5708;-3.1416;-3.1416;-1.5708;-1.5708;1.5708;-3.1416;-0.7854;0.7854;-1.5708;-1.5708;0;0;0;0;0;0;0];
cMech_DSuite_nostairs = [-1.5708;-0.7854;0.7854;-1.5708;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
dDSuite_Mech_nostairs = [1.5708;0.7854;-0.7854;1.5708;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];

ma = zeros(length(aMech_DSuite_stairs),1);
na = zeros((length(aMech_DSuite_stairs)-length(dir_vec)+1),1);
for i = 1:(length(aMech_DSuite_stairs)-length(dir_vec)+1);
    for j = 1:length(dir_vec);
        if dir_vec(j) == aMech_DSuite_stairs(j+i-1);
            ma(j) = 1;
        else;
            ma(j) = 0;
        end;
        na(i) = mean(ma);
    end;
end;    
na        ;
    
mb = zeros(length(bDSuite_Mech_stairs),1);
nb = zeros((length(bDSuite_Mech_stairs)-length(dir_vec)+1),1);
for i = 1:(length(bDSuite_Mech_stairs)-length(dir_vec)+1);
    for j = 1:length(dir_vec);
        if dir_vec(j) == bDSuite_Mech_stairs(j+i-1);
            mb(j) = 1;
        else;
            mb(j) = 0;
        end;
        nb(i) = mean(mb);
    end;
end;    
nb            ;

mc = zeros(length(cMech_DSuite_nostairs),1);
nc = zeros((length(cMech_DSuite_nostairs)-length(dir_vec)+1),1);
for i = 1:(length(cMech_DSuite_nostairs)-length(dir_vec)+1);
    for j = 1:length(dir_vec);
        if dir_vec(j) == cMech_DSuite_nostairs(j+i-1);
            mc(j) = 1;
        else;
            mc(j) = 0;
        end;
        nc(i) = mean(mc);
    end;
end;    
nc      ;

md = zeros(length(dDSuite_Mech_nostairs),1);
nd = zeros((length(dDSuite_Mech_nostairs)-length(dir_vec)+1),1);
for i = 1:(length(dDSuite_Mech_nostairs)-length(dir_vec)+1);
    for j = 1:length(dir_vec);
        if dir_vec(j) == dDSuite_Mech_nostairs(j+i-1);
            md(j) = 1;
        else
            md(j) = 0;
        end;
        nd(i) = mean(md);
    end;
end;    
nd   ;

Route_factor = [max(na);max(nb);max(nc);max(nd)]

figure
plot(avg_Zgiro(:,1),avg_Zgiro(:,2),'r-', T,Theta(:,1),'g-',avg_Zgiro(:,1),avg_Zgiro(:,3),'b-',avg_Zgiro(:,1),avg_Zgiro(:,4),'cy-',avg_Zgiro(:,1),avg_Zgiro(:,8),'k-',avg_Zgiro(:,1),avg_Zgiro(:,6),'m-')
figure
plot(T,Theta(:,1),'r-',avg_Zgiro(:,1),avg_Zgiro(:,10),'k-')%,avg_Zgiro(:,1),avg_Zgiro(:,8),'m-',avg_Zgiro(:,1),avg_Zgiro(:,9),'cy-',




