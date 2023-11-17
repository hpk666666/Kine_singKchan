function [P0,P1,P2,P3,P4]=Pe_Vt(u,Vt,Ve,data,n)
P0=0;P1=0;P2=0;P3=0;P4=0;
for i=1:n+1
    P0=data(1,i)*(Vt-Ve).^(i+u-1+4)+P0;
    P1=data(2,i)*(Vt-Ve).^(i+u-1+3)+P1;
    P2=data(3,i)*(Vt-Ve).^(i+u-1+2)+P2;
    P3=data(4,i)*(Vt-Ve).^(i+u-1+1)+P3;
    P4=data(5,i)*(Vt-Ve).^(i+u-1)+P4;
end
end