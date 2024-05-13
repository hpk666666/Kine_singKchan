function P=leijia2_V_noise(Vt,V0,data,n)
P0=zeros(1,length(Vt));%初始化Pe(0,Vt),Pe(1,Vt)
P1=zeros(1,length(Vt));%初始化Pe(0,Vt),Pe(1,Vt)
P2=zeros(1,length(Vt));%初始化Pe(0,Vt),Pe(1,Vt)
P3=zeros(1,length(Vt));%初始化Pe(0,Vt),Pe(1,Vt)
P4=zeros(1,length(Vt));%初始化Pe(0,Vt),Pe(1,Vt)
P=cell(1,5);
for j=1:length(Vt)
    for i=1:n+1
        P0(j)=data(1,i)*(Vt(j)-V0)^(i-1)+P0(j);
        P1(j)=data(2,i)*(Vt(j)-V0)^(i-1)+P1(j);
        P2(j)=data(3,i)*(Vt(j)-V0)^(i-1)+P2(j);
        P3(j)=data(4,i)*(Vt(j)-V0)^(i-1)+P3(j);
        P4(j)=data(5,i)*(Vt(j)-V0)^(i-1)+P4(j);
    end
end
P{1,1}=P0;
P{1,2}=P1;
P{1,3}=P2;
P{1,4}=P3;
P{1,5}=P4;
end