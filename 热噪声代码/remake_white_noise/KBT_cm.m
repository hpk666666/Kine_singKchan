syms x m  
% Iex=1;%外部电流，单位 fA 法安   1 A = 10^15 fA ，
n=50;
gK=20;%单位是电导率 S
K_b=1.380649*10^-23;%玻尔兹曼常数,单位为J/K，J为焦耳，K是热力学温度
T=294.15;%绝对温度，单位K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=linspace(0.02*10^-12,0.08*10^-12,n);%定义膜面积，单位是um^2
%iex=1,M=linspace(0.02,0.08,n)
WHITENOSIE=zeros(1,n);
for i=1:n
Cm=0.01*M(i);%单位是电容 F
WHITENOSIE(i)=K_b*T/Cm;
end
plot(M,WHITENOSIE);