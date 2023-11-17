function [P0,P1,P2,P3,P4]=leijia2_VL(u,Vt,VL,data,n)
P0=zeros(1,length(Vt));%初始化Pe(0,Vt),Pe(1,Vt)
P1=zeros(1,length(Vt));%初始化Pe(0,Vt),Pe(1,Vt)
P2=zeros(1,length(Vt));%初始化Pe(0,Vt),Pe(1,Vt)
P3=zeros(1,length(Vt));%初始化Pe(0,Vt),Pe(1,Vt)
P4=zeros(1,length(Vt));%初始化Pe(0,Vt),Pe(1,Vt)
for j=1:length(Vt)
    for i=1:n+1
        P0(j)=data(1,i)*(VL-Vt(j)).^(i+u-1)+P0(j);
        P1(j)=data(2,i)*(VL-Vt(j)).^(i+u-1)+P1(j);
        P2(j)=data(3,i)*(VL-Vt(j)).^(i+u-1)+P2(j);
        P3(j)=data(4,i)*(VL-Vt(j)).^(i+u-1)+P3(j);
        P4(j)=data(5,i)*(VL-Vt(j)).^(i+u-1+1)+P4(j);
    end
end
end