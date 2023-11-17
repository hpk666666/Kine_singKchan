syms x m  
Iex=0.02*10^-12;%�ⲿ��������λ A
n=20;%��H_20//100
gK=20*10^-12;%��λ�ǵ絼�� S
VL1=-54.4*10^-3;%��λ�� V
VK=-77*10^-3;%��λ�� V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=linspace(0.1*10^-12,100*10^-12,n);%����Ĥ�������λ��m^2//0.05-3
ExC=zeros(1,n);%��ֵ
Ex2C=zeros(1,n);%���׾�
ExO=zeros(1,n);%��ֵ
Ex2O=zeros(1,n);%���׾�
Dx=zeros(1,n);%����
Popen=zeros(1,n);%�򿪸���
Pclose=zeros(1,n);%�رո���
PEDC=zeros(1,n);%ͨ���ر�ʱ���й¶��������
PEDO=zeros(1,n);%ͨ����ʱ���й¶��������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
gL=3*M(i);%��λ�ǵ絼�� S/m ^2
ge=gL+gK;%��λ�ǵ絼�� S
Cm=0.01*M(i);%��λ�ǵ��� F
E=ge/Cm;%��λ�� pS/pF
L=gL/Cm;%��λ�� pS/pF
VL_star=(gL*VL1+Iex)/gL;%��λ�� V
Ve=(gL*VL1+gK*VK+Iex)/ge;%��λ�� V
faL=taylor(0.01*10^6*(x+55*10^-3)/(1-exp(-(x+55*10^-3)/0.01)),x,VL_star,'order',10);%�ٶȦ�(v)��Ve������̩��չʽ
fbL=taylor(0.125*10^3*exp(-(x+65*10^-3)/0.08),x,VL_star,'order',10);%�ٶȦ�(v)��Ve������̩��չʽ
newfaL=subs(faL,(x-VL_star),m);%��(x-VL)�滻��m
newfbL=subs(fbL,(x-VL_star),m);%��(x-VL)�滻��m
KaL=vpa(coeffs(newfaL));%����������̩��չʽ��ϵ��
LbL=vpa(coeffs(newfbL));%����������̩��չʽ��ϵ��
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
%E*(Ve-VL)*(u+1)->(u+1)E*(Ve-VL)�Խ��û��Ӱ��Ϊ��ʡ����ʱ�䣬��Ĩȥ��
s=real(vpa(solve(det(S),'u')));
s1=double(vpa(s));
u_VL=[s1(2),s1(3),s1(4),s1(5)];%��ȥ��һ��-1��uֵ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fae=taylor(0.01*10^6*(x+55*10^-3)/(1-exp(-(x+55*10^-3)/0.01)),x,Ve,'order',10);%�ٶȦ�(v)��Ve������̩��չʽ
fbe=taylor(0.125*10^3*exp(-(x+65*10^-3)/0.08),x,Ve,'order',10);%�ٶȦ�(v)��Ve������̩��չʽ
newfa=subs(fae,(x-Ve),m);%��(x-Ve)�滻��m
newfb=subs(fbe,(x-Ve),m);%��(x-Ve)�滻��m
Kae=vpa(coeffs(newfa));%����������̩��չʽ��ϵ��
Lbe=vpa(coeffs(newfb));%����������̩��չʽ��ϵ��
Ka_e=zeros(1,n+1);
Lb_e=zeros(1,n+1);
Ka_e(1:10)=Kae;
Lb_e(1:10)=Lbe;
u_Ve=4*Lb_e(1)/E-1;%��ֵ,����Ĥ����仯
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short e
data_1=VL_diedai(1,u_VL,n,Ve,VL_star,E,Ka_L,Lb_L,L);%1,2,3,4�ֱ�����ĸ���ͬ��uֵ��pi
data_2=VL_diedai(2,u_VL,n,Ve,VL_star,E,Ka_L,Lb_L,L);
data_3=VL_diedai(3,u_VL,n,Ve,VL_star,E,Ka_L,Lb_L,L);
data_4=VL_diedai(4,u_VL,n,Ve,VL_star,E,Ka_L,Lb_L,L);
data=Ve_diedai(u_Ve,n,Ve,VL_star,E,Ka_e,Lb_e,L);
%%%%%%%%%%%%%%%%%%ϵ��K%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vt=Ve+(VL_star-Ve)*3/4;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Pe_0,Pe_1,Pe_2,Pe_3,Pe_4]=Pe_Vt(u_Ve,Vt,Ve,data,n);
[PL1_0,PL1_1,PL1_2,PL1_3,PL1_4]=PL_Vt(u_VL(1),Vt,VL_star,data_1,n);
[PL2_0,PL2_1,PL2_2,PL2_3,PL2_4]=PL_Vt(u_VL(2),Vt,VL_star,data_2,n);
[PL3_0,PL3_1,PL3_2,PL3_3,PL3_4]=PL_Vt(u_VL(3),Vt,VL_star,data_3,n);
[PL4_0,PL4_1,PL4_2,PL4_3,PL4_4]=PL_Vt(u_VL(4),Vt,VL_star,data_4,n);
Pe_int=Ve_int(u_Ve,Vt,Ve,data,n);%��������ţ��������Ĺ�ʽ��(Pe(0,Vt)+...Pe(4,Vt))�ڡ�Ve,Vt���Ļ���
PL_int_1=VL_int(u_VL(1),Vt,VL_star,data_1,n);%��������ţ��������Ĺ�ʽ��(PL1(0,Vt)+...PL1(4,Vt))�ڡ�Vt,VL���Ļ���
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
%E_TL2_1����ڶ�������1״̬��VL��Vt�ľ�ֵ
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
%E2_TL2_1����ڶ�������1״̬��VL��Vt�Ķ��׾�
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