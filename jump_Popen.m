syms x m  
Iex=0.02*10^-12;%外部电流，单位 A
n=20;%求H_20
gK=20*10^-12;%单位是电导率 S
VL1=-54.4*10^-3;%单位是 V
VK=-77*10^-3;%单位是 V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=linspace(0.05*10^-12,3*10^-12,n);%定义膜面积，单位是m^2
ExJ=zeros(1,n);%均值
Ex2J=zeros(1,n);%二阶矩
DxJ=zeros(1,n);%方差
PopenJ=zeros(1,n);%打开概率
for i=1:n
gL=3*M(i);%单位是电导率 S/m ^2
ge=gL+gK;%单位是电导率 S
Cm=0.01*M(i);%单位是电容 F
E=ge/Cm;%单位是 pS/pF
L=gL/Cm;%单位是 pS/pF
VL_star=(gL*VL1+Iex)/gL;%单位是 V
Ve=(gL*VL1+gK*VK+Iex)/ge;%单位是 V   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
faL=0.01*10^6*(VL_star+55*10^-3)/(1-exp(-(VL_star+55*10^-3)/0.01));%速度α(v)在VL这个点的值
fbL=0.125*10^3*exp(-(VL_star+65*10^-3)/0.08);%速度β(v)在VL这个点的值
fbe=0.125*10^3*exp(-(Ve+65*10^-3)/0.08);%速度β(v)在Ve这个点的值
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