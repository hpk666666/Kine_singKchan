function [S0,S1,S2,S3,S4]=jifen(nP0,nP1,nP2,nP3,X)
S0=0;%��ʼ��Pe(0,Vt),Pe(1,Vt)
S1=0;%��ʼ��Pe(0,Vt),Pe(1,Vt)
S2=0;%��ʼ��Pe(0,Vt),Pe(1,Vt)
S3=0;%��ʼ��Pe(0,Vt),Pe(1,Vt)
for j=1:length(X)-1%ȥ����ǰ�������������쳣��
S0=(nP0(j)+nP0(j+1))*(X(j+1)-X(j))/2+S0;
S1=(nP1(j)+nP1(j+1))*(X(j+1)-X(j))/2+S1;
S2=(nP2(j)+nP2(j+1))*(X(j+1)-X(j))/2+S2;
S3=(nP3(j)+nP3(j+1))*(X(j+1)-X(j))/2+S3;
end
end