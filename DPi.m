function [Pi]=DPi(dP0,dP1,dP2,dP3,dP4,dX)
Pi=zeros(1,length(dP0));
S0=0;%初始化Pe(0,Vt),Pe(1,Vt)
S1=0;%初始化Pe(0,Vt),Pe(1,Vt)
S2=0;%初始化Pe(0,Vt),Pe(1,Vt)
S3=0;%初始化Pe(0,Vt),Pe(1,Vt)
S4=0;%初始化Pe(0,Vt),Pe(1,Vt)
for i=2:length(dP0) %去掉最前面和最后面俩个异常点.第一个点积分为0
for j=2:i %去掉最前面和最后面俩个异常点
S0=(dP0(j-1)+dP0(j))*(dX(j)-dX(j-1))/2+S0;
S1=(dP1(j-1)+dP1(j))*(dX(j)-dX(j-1))/2+S1;
S2=(dP2(j-1)+dP2(j))*(dX(j)-dX(j-1))/2+S2;
S3=(dP3(j-1)+dP3(j))*(dX(j)-dX(j-1))/2+S3;
S4=(dP4(j-1)+dP4(j))*(dX(j)-dX(j-1))/2+S4;
end
Pi(i)=S1+S2+S3+S4+S0;
end