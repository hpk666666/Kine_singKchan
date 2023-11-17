function [DP1,DP2,DP3,DP4,DP5]=DnewP(dP0,dP1,dP2,dP3,dP4,Pi)
DP1=zeros(1,length(dP0));%初始化Pe(0,Vt),Pe(1,Vt)
DP2=zeros(1,length(dP0));%初始化Pe(0,Vt),Pe(1,Vt)
DP3=zeros(1,length(dP0));%初始化Pe(0,Vt),Pe(1,Vt)
DP4=zeros(1,length(dP0));%初始化Pe(0,Vt),Pe(1,Vt)
DP5=zeros(1,length(dP0));%初始化Pe(0,Vt),Pe(1,Vt)
PV=zeros(1,length(dP0));
for j=2:length(dP0)%1时pi为0在分母故去掉
        PV(j)=dP0(j)+dP1(j)+dP2(j)+dP3(j)+dP4(j);
        DP1(j)=dP0(j)*log(dP0(j)/(PV(j)*Pi(j)));
        DP2(j)=dP1(j)*log(dP1(j)/(PV(j)*Pi(j)));
        DP3(j)=dP2(j)*log(dP2(j)/(PV(j)*Pi(j)));
        DP4(j)=dP3(j)*log(dP3(j)/(PV(j)*Pi(j)));
        DP5(j)=dP4(j)*log(dP4(j)/(PV(j)*Pi(j)));
end
end