syms x m  
Iex=0.02*10^-12;%�ⲿ��������λ A
n=20;%��H_20
gK=20*10^-12;%��λ�ǵ絼�� S
VL1=-54.4*10^-3;%��λ�� V
VK=-77*10^-3;%��λ�� V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=linspace(0.05*10^-12,3*10^-12,n);%����Ĥ�������λ��m^2
ExJ=zeros(1,n);%��ֵ
Ex2J=zeros(1,n);%���׾�
DxJ=zeros(1,n);%����
PopenJ=zeros(1,n);%�򿪸���
for i=1:n
gL=3*M(i);%��λ�ǵ絼�� S/m ^2
ge=gL+gK;%��λ�ǵ絼�� S
Cm=0.01*M(i);%��λ�ǵ��� F
E=ge/Cm;%��λ�� pS/pF
L=gL/Cm;%��λ�� pS/pF
VL_star=(gL*VL1+Iex)/gL;%��λ�� V
Ve=(gL*VL1+gK*VK+Iex)/ge;%��λ�� V   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
faL=0.01*10^6*(VL_star+55*10^-3)/(1-exp(-(VL_star+55*10^-3)/0.01));%�ٶȦ�(v)��VL������ֵ
fbL=0.125*10^3*exp(-(VL_star+65*10^-3)/0.08);%�ٶȦ�(v)��VL������ֵ
fbe=0.125*10^3*exp(-(Ve+65*10^-3)/0.08);%�ٶȦ�(v)��Ve������ֵ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pz=1+4*faL/fbL+6*faL^2/fbL^2+4*faL^3/fbL^3+faL^4/(fbL^3*fbe);
PopenJ(i)=faL^4/(fbL^3*fbe)/Pz;
ExJ(i)=VL_star*(1-PopenJ(i))+Ve*PopenJ(i);
Ex2J(i)=VL_star^2*(1-PopenJ(i))+Ve^2*PopenJ(i);
DxJ(i)=Ex2J(i)-ExJ(i)^2;
end
fig2=figure(2);
subplot(1,3,1);
plot(M,ExJ)
xlabel('M');ylabel('E(V)');
subplot(1,3,2);
plot(M,DxJ)
xlabel('M');ylabel('D(V)');
subplot(1,3,3);
plot(M,PopenJ)
xlabel('M');ylabel('Popen');