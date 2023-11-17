function [nP1,nP2,nP3,nP4,nP5]=newP(P0,P1,P2,P3,P4)
nP1=zeros(1,length(P0));%初始化Pe(0,Vt),Pe(1,Vt)
nP2=zeros(1,length(P0));%初始化Pe(0,Vt),Pe(1,Vt)
nP3=zeros(1,length(P0));%初始化Pe(0,Vt),Pe(1,Vt)
nP4=zeros(1,length(P0));%初始化Pe(0,Vt),Pe(1,Vt)
nP5=zeros(1,length(P0));%初始化Pe(0,Vt),Pe(1,Vt)
for j=1:length(P0)
        nP1(j)=P0(j)*log(P0(j));
        nP2(j)=P1(j)*log(P1(j));
        nP3(j)=P2(j)*log(P2(j));
        nP4(j)=P3(j)*log(P3(j));
        nP5(j)=P4(j)*log(P4(j));
end
end