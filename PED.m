syms x m  
Iex=0.02*10^-12;%外部电流，单位 A
n=20;%求H_20//100
gK=20*10^-12;%单位是电导率 S
VL1=-54.4*10^-3;%单位是 V
VK=-77*10^-3;%单位是 V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=linspace(0.1*10^-12,100*10^-12,n);%定义膜面积，单位是m^2//0.05-3
ExC=zeros(1,n);%均值
Ex2C=zeros(1,n);%二阶矩
ExO=zeros(1,n);%均值
Ex2O=zeros(1,n);%二阶矩
Dx=zeros(1,n);%方差
Popen=zeros(1,n);%打开概率
Pclose=zeros(1,n);%关闭概率
PEDC=zeros(1,n);%通道关闭时候的泄露电流功耗
PEDO=zeros(1,n);%通道打开时候的泄露电流功耗
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
gL=3*M(i);%单位是电导率 S/m ^2
ge=gL+gK;%单位是电导率 S
Cm=0.01*M(i);%单位是电容 F
E=ge/Cm;%单位是 pS/pF
L=gL/Cm;%单位是 pS/pF
VL_star=(gL*VL1+Iex)/gL;%单位是 V
Ve=(gL*VL1+gK*VK+Iex)/ge;%单位是 V
faL=taylor(0.01*10^6*(x+55*10^-3)/(1-exp(-(x+55*10^-3)/0.01)),x,VL_star,'order',10);%速度α(v)在Ve这个点的泰勒展式
fbL=taylor(0.125*10^3*exp(-(x+65*10^-3)/0.08),x,VL_star,'order',10);%速度β(v)在Ve这个点的泰勒展式
newfaL=subs(faL,(x-VL_star),m);%将(x-VL)替换成m
newfbL=subs(fbL,(x-VL_star),m);%将(x-VL)替换成m
KaL=vpa(coeffs(newfaL));%按升幂排列泰勒展式的系数
LbL=vpa(coeffs(newfbL));%按升幂排列泰勒展式的系数
for j=1:10
    if mod(j,2)==0
        KaL(j)=-KaL(j);
        LbL(j)=-LbL(j);
    end
end
Ka_L=zeros(1,n+1);
Lb_L=zeros(1,n+1);
Ka_L(1:10)=KaL;
Lb_L(1:10)=LbL;
syms u
S=[L*(u+1)-4*Ka_L(1),Lb_L(1),                  0,                      0                       0;
   4*Ka_L(1),        L*(u+1)-(3*Ka_L(1)+Lb_L(1)),2*Lb_L(1),                0,                      0;
   0,              3*Ka_L(1),                L*(u+1)-2*(Ka_L(1)+Lb_L(1)),3*Lb_L(1),                0;
   0,              0,                      2*(Ka_L(1)),              L*(u+1)-(Ka_L(1)+3*Lb_L(1)),0;
   0,              0,                      0,                      Ka_L(1),                  (u+1)];
%E*(Ve-VL)*(u+1)->(u+1)E*(Ve-VL)对结果没有影响为节省计算时间，故抹去。
s=real(vpa(solve(det(S),'u')));
s1=double(vpa(s));
u_VL=[s1(2),s1(3),s1(4),s1(5)];%除去第一个-1的u值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fae=taylor(0.01*10^6*(x+55*10^-3)/(1-exp(-(x+55*10^-3)/0.01)),x,Ve,'order',10);%速度α(v)在Ve这个点的泰勒展式
fbe=taylor(0.125*10^3*exp(-(x+65*10^-3)/0.08),x,Ve,'order',10);%速度β(v)在Ve这个点的泰勒展式
newfa=subs(fae,(x-Ve),m);%将(x-Ve)替换成m
newfb=subs(fbe,(x-Ve),m);%将(x-Ve)替换成m
Kae=vpa(coeffs(newfa));%按升幂排列泰勒展式的系数
Lbe=vpa(coeffs(newfb));%按升幂排列泰勒展式的系数
Ka_e=zeros(1,n+1);
Lb_e=zeros(1,n+1);
Ka_e(1:10)=Kae;
Lb_e(1:10)=Lbe;
u_Ve=4*Lb_e(1)/E-1;%μ值,会随膜面积变化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short e
data_1=VL_diedai(1,u_VL,n,Ve,VL_star,E,Ka_L,Lb_L,L);%1,2,3,4分别代表四个不同的u值的pi
data_2=VL_diedai(2,u_VL,n,Ve,VL_star,E,Ka_L,Lb_L,L);
data_3=VL_diedai(3,u_VL,n,Ve,VL_star,E,Ka_L,Lb_L,L);
data_4=VL_diedai(4,u_VL,n,Ve,VL_star,E,Ka_L,Lb_L,L);
data=Ve_diedai(u_Ve,n,Ve,VL_star,E,Ka_e,Lb_e,L);
%%%%%%%%%%%%%%%%%%系数K%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vt=Ve+(VL_star-Ve)*3/4;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Pe_0,Pe_1,Pe_2,Pe_3,Pe_4]=Pe_Vt(u_Ve,Vt,Ve,data,n);
[PL1_0,PL1_1,PL1_2,PL1_3,PL1_4]=PL_Vt(u_VL(1),Vt,VL_star,data_1,n);
[PL2_0,PL2_1,PL2_2,PL2_3,PL2_4]=PL_Vt(u_VL(2),Vt,VL_star,data_2,n);
[PL3_0,PL3_1,PL3_2,PL3_3,PL3_4]=PL_Vt(u_VL(3),Vt,VL_star,data_3,n);
[PL4_0,PL4_1,PL4_2,PL4_3,PL4_4]=PL_Vt(u_VL(4),Vt,VL_star,data_4,n);
Pe_int=Ve_int(u_Ve,Vt,Ve,data,n);%这里利用牛顿莱布尼茨公式算(Pe(0,Vt)+...Pe(4,Vt))在【Ve,Vt】的积分
PL_int_1=VL_int(u_VL(1),Vt,VL_star,data_1,n);%这里利用牛顿莱布尼茨公式算(PL1(0,Vt)+...PL1(4,Vt))在【Vt,VL】的积分
PL_int_2=VL_int(u_VL(2),Vt,VL_star,data_2,n);
PL_int_3=VL_int(u_VL(3),Vt,VL_star,data_3,n);
PL_int_4=VL_int(u_VL(4),Vt,VL_star,data_4,n);
 A=[ -Pe_0, PL1_0,   PL2_0,   PL3_0,   PL4_0;
     -Pe_1, PL1_1,   PL2_1,   PL3_1,   PL4_1;
     -Pe_2, PL1_2,   PL2_2,   PL3_2,   PL4_2;
     -Pe_3, PL1_3,   PL2_3,   PL3_3,   PL4_3;
     -Pe_4, PL1_4,   PL2_4,   PL3_4,   PL4_4;
     Pe_int,PL_int_1,PL_int_2,PL_int_3,PL_int_4];
 B=[0;0;0;0;0;1];
 format rat
 K=A\B;
 %%%%%%%%%%E(X)%%%%%%%%%%%%%%%%%%%%%
[E_eT_0,E_eT_1,E_eT_2,E_eT_3,E_eT_4]=Ee(u_Ve,Vt,Ve,data,n);
[E_TL1_0,E_TL1_1,E_TL1_2,E_TL1_3,E_TL1_4]=EL(u_VL(1),Vt,VL_star,data_1,n);
[E_TL2_0,E_TL2_1,E_TL2_2,E_TL2_3,E_TL2_4]=EL(u_VL(2),Vt,VL_star,data_2,n);
[E_TL3_0,E_TL3_1,E_TL3_2,E_TL3_3,E_TL3_4]=EL(u_VL(3),Vt,VL_star,data_3,n);
[E_TL4_0,E_TL4_1,E_TL4_2,E_TL4_3,E_TL4_4]=EL(u_VL(4),Vt,VL_star,data_4,n);
%E_TL2_1代表第二组数据1状态，VL到Vt的均值
ExC(i)=K(1)*(E_eT_0+E_eT_1+E_eT_2+E_eT_3)+...
   K(2)*(E_TL1_0+E_TL1_1+E_TL1_2+E_TL1_3)+...
   K(3)*(E_TL2_0+E_TL2_1+E_TL2_2+E_TL2_3)+...
   K(4)*(E_TL3_0+E_TL3_1+E_TL3_2+E_TL3_3)+...
   K(5)*(E_TL4_0+E_TL4_1+E_TL4_2+E_TL4_3);
ExO(i)=K(1)*(E_eT_4)+...
   K(2)*(E_TL1_4)+...
   K(3)*(E_TL2_4)+...
   K(4)*(E_TL3_4)+...
   K(5)*(E_TL4_4);
 %%%%%%%%%%D(X)%%%%%%%%%%%%%%%%%%%%%
[E2_eT_0,E2_eT_1,E2_eT_2,E2_eT_3,E2_eT_4]=Ee2(u_Ve,Vt,Ve,data,n);
[E2_TL1_0,E2_TL1_1,E2_TL1_2,E2_TL1_3,E2_TL1_4]=EL2(u_VL(1),Vt,VL_star,data_1,n);
[E2_TL2_0,E2_TL2_1,E2_TL2_2,E2_TL2_3,E2_TL2_4]=EL2(u_VL(2),Vt,VL_star,data_2,n);
[E2_TL3_0,E2_TL3_1,E2_TL3_2,E2_TL3_3,E2_TL3_4]=EL2(u_VL(3),Vt,VL_star,data_3,n);
[E2_TL4_0,E2_TL4_1,E2_TL4_2,E2_TL4_3,E2_TL4_4]=EL2(u_VL(4),Vt,VL_star,data_4,n);
%E2_TL2_1代表第二组数据1状态，VL到Vt的二阶矩
Ex2C(i)=K(1)*(E2_eT_0+E2_eT_1+E2_eT_2+E2_eT_3)+...
   K(2)*(E2_TL1_0+E2_TL1_1+E2_TL1_2+E2_TL1_3)+...
   K(3)*(E2_TL2_0+E2_TL2_1+E2_TL2_2+E2_TL2_3)+...
   K(4)*(E2_TL3_0+E2_TL3_1+E2_TL3_2+E2_TL3_3)+...
   K(5)*(E2_TL4_0+E2_TL4_1+E2_TL4_2+E2_TL4_3);
Ex2O(i)=K(1)*(E2_eT_4)+...
   K(2)*(E2_TL1_4)+...
   K(3)*(E2_TL2_4)+...
   K(4)*(E2_TL3_4)+...
   K(5)*(E2_TL4_4);
u_all=[u_Ve,u_VL];
% Pclose(i)=P_close(u_all,Vt,Ve,VL,data,data_1,data_2,data_3,data_4,n,K);
Popen(i)=P_open(u_all,Vt,Ve,VL_star,data,data_1,data_2,data_3,data_4,n,K);
Pclose(i)=1-Popen(i);
PEDC(i)=gL*(Ex2C(i)-2*VL_star*ExC(i)+VL_star^2*Pclose(i));
PEDO(i)=ge*(Ex2O(i)-2*Ve*ExO(i)+Ve^2*Popen(i));
end
% subplot(1,2,1);
% plot(M,Popen)
% xlabel('M');ylabel('Popen');
% subplot(1,2,2);
% plot(M,Pclose)
% xlabel('M');ylabel('Pclose');
subplot(1,3,1);
plot(M,PEDC)
xlabel('M');ylabel('PEDC');
subplot(1,3,2);
plot(M,PEDO)
xlabel('M');ylabel('PEDO');
subplot(1,3,3);
plot(M,PEDC+PEDO)
xlabel('M');ylabel('PED');