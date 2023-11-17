function [FP01,FP12,FP23,FP34]=FnewP(P0,P1,P2,P3,P4,X,KT)
FP01=zeros(1,length(P0));%初始化Pe(0,Vt),Pe(1,Vt)
FP12=zeros(1,length(P0));%初始化Pe(0,Vt),Pe(1,Vt)
FP23=zeros(1,length(P0));%初始化Pe(0,Vt),Pe(1,Vt)
FP34=zeros(1,length(P0));%初始化Pe(0,Vt),Pe(1,Vt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alpha_V=0;
% Beta_V=0;
% fa=0.01*10^6*(x+55*10^-3)/(1-exp(-(x+55*10^-3)/0.01));%速度α(v)
% fb=0.125*10^3*exp(-(x+65*10^-3)/0.08);%速度β(v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
for j=1:length(P0)%1时pi为0在分母故去掉
          Alpha_V=0.01*10^6*(X(j)+55*10^-3)/(1-exp(-(X(j)+55*10^-3)/0.01));
          Beta_V=0.125*10^3*exp(-(X(j)+65*10^-3)/0.08);
          FP01(j)=KT*(P0(j)*4*Alpha_V-P1(j)*1*Beta_V)*log(4*Alpha_V/(1*Beta_V));
          FP12(j)=KT*(P1(j)*3*Alpha_V-P2(j)*2*Beta_V)*log(3*Alpha_V/(2*Beta_V));
          FP23(j)=KT*(P2(j)*2*Alpha_V-P3(j)*3*Beta_V)*log(2*Alpha_V/(3*Beta_V));
          FP34(j)=KT*(P3(j)*1*Alpha_V-P4(j)*4*Beta_V)*log(1*Alpha_V/(4*Beta_V));
end
end