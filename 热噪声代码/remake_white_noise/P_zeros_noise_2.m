function P=P_zeros_noise_2(Vt,V0,data,n)
P0=0;P1=0;P2=0;P3=0;P4=0;
P = [0;0;0;0;0];
for i=1:n+1
    P0=data(1,i)*(i-1)*(Vt-V0)^(i-2)+P0;
    P1=data(2,i)*(i-1)*(Vt-V0)^(i-2)+P1;
    P2=data(3,i)*(i-1)*(Vt-V0)^(i-2)+P2;
    P3=data(4,i)*(i-1)*(Vt-V0)^(i-2)+P3;
    P4=data(5,i)*(i-1)*(Vt-V0)^(i-2)+P4;
end
P = [P0;P1;P2;P3;P4];
end