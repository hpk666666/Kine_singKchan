syms x m  
M=8*10^-12;%����Ĥ�������λ��m^2%%% 2 4 6 8 
%�Ǵ����� M = 1*10-12 m^2 �����ֱ�Ϊ 0.05*10^-12 0.1*10^-12 0.15*10^-12 0.2*10^-12 A �����
Iex=0.6*10^-12;%�ⲿ��������λ A%%%-0.006-----0.1
%дһ��ǰ�ô����H��Iex�Ĳ������꣬�������ݴ�������������ex֮����ٵ������ݣ�
gL=3*M;%��λ�ǵ絼�� S/m ^2
gK=20*10^-12;%��λ�ǵ絼�� S
VL1=-54.4*10^-3;%��λ�� V
VK=-77*10^-3;%��λ�� V
ge=gL+gK;%��λ�ǵ絼�� S
Cm=0.01*M;%��λ�ǵ��� F
E=ge/Cm;%��λ�� pS/pF
L=gL/Cm;%��λ�� pS/pF
VL_star=(gL*VL1+Iex)/gL;%��λ�� V
Ve=(gL*VL1+gK*VK+Iex)/ge;%��λ�� V
% M=20*10^-12;%����Ĥ�������λ��m^2%%%0.5
% Iex=0.*10^-12;%�ⲿ��������λ A%%%-0.006-----0.1
% VL_star=(gL*VL1+Iex)/gL-10^-10;%��λ�� V
% Ve=(gL*VL1+gK*VK+Iex)/ge+10^-10;%��λ�� V
fa=taylor(0.01*10^6*(x+55*10^-3)/(1-exp(-(x+55*10^-3)/0.01)),x,VL_star,'order',10);%�ٶȦ�(v)��Ve������̩��չʽ
fb=taylor(0.125*10^3*exp(-(x+65*10^-3)/0.08),x,VL_star,'order',10);%�ٶȦ�(v)��Ve������̩��չʽ
newfa=subs(fa,(x-VL_star),m);%��(x-VL_star)�滻��m
newfb=subs(fb,(x-VL_star),m);%��(x-VL_star)�滻��m
Kaa=vpa(coeffs(newfa));%����������̩��չʽ��ϵ��
Lbb=vpa(coeffs(newfb));%����������̩��չʽ��ϵ��
Vt=Ve+(VL_star-Ve)*2.9/4;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:10
    if mod(i,2)==0
        Kaa(i)=-Kaa(i);
        Lbb(i)=-Lbb(i);
    end
end
n=200;%��H_20
Ka=zeros(1,n+1);
Lb=zeros(1,n+1);
Ka(1:10)=Kaa;
Lb(1:10)=Lbb;
VL_star=(gL*VL1+Iex)/gL-10^-11;%��λ�� V
Ve=(gL*VL1+gK*VK+Iex)/ge+10^-11;%��λ�� V
syms u
S=[L*(u+1)-4*Ka(1),Lb(1),                  0,                      0                       0;
   4*Ka(1),        L*(u+1)-(3*Ka(1)+Lb(1)),2*Lb(1),                0,                      0;
   0,              3*Ka(1),                L*(u+1)-2*(Ka(1)+Lb(1)),3*Lb(1),                0;
   0,              0,                      2*(Ka(1)),              L*(u+1)-(Ka(1)+3*Lb(1)),0;
   0,              0,                      0,                      Ka(1),                  E*(Ve-VL_star)*(u+1)];
s=real(vpa(solve(det(S),'u')));
%roots[a,b,c,d,e]
s1=double(vpa(s));
%f=a*x^4+b*x^3+c*x^2+d*x+e;
%p=sym2poly(f)
%X=roots(p)
u_VL=[s1(2),s1(3),s1(4),s1(5)];%��ȥ��һ��-1��uֵ
%�ĸ�uֵ�ý��ƽ������x�ǣ�1-u��
%L0=Lb(1),K0=Ka(1);
%1,2,3,4�ֱ�����ĸ���ͬ��uֵ
format short e
data_1=VL_diedai(1,u_VL,n,Ve,VL_star,E,Ka,Lb,L);%1,2,3,4�ֱ�����ĸ���ͬ��uֵ��pi
data_2=VL_diedai(2,u_VL,n,Ve,VL_star,E,Ka,Lb,L);
data_3=VL_diedai(3,u_VL,n,Ve,VL_star,E,Ka,Lb,L);
data_4=VL_diedai(4,u_VL,n,Ve,VL_star,E,Ka,Lb,L);

%data=zeros(2,5);%����һ��2��5�е�0����
%data(1,:)=[1,2,3,4,5];%�ڵ�1����������
%data(1,:)=H0
%��Ϊ�±겻��Ϊ0���±�ȫ����1
%��������inv(x)
%����˷���C = A*B ��ʾ���� A �� B �����Դ����˻���A ������������ B ��������ȡ�
%����ת�ã�A.' ��ʾ A ������ת�á����ڸ������ⲻ�漰���

fa=taylor(0.01*10^6*(x+55*10^-3)/(1-exp(-(x+55*10^-3)/0.01)),x,Ve,'order',10);%�ٶȦ�(v)��Ve������̩��չʽ
fb=taylor(0.125*10^3*exp(-(x+65*10^-3)/0.08),x,Ve,'order',10);%�ٶȦ�(v)��Ve������̩��չʽ
newfa=subs(fa,(x-Ve),m);%��(x-Ve)�滻��m
newfb=subs(fb,(x-Ve),m);%��(x-Ve)�滻��m
Kaa=vpa(coeffs(newfa));%����������̩��չʽ��ϵ��
Lbb=vpa(coeffs(newfb));%����������̩��չʽ��ϵ��
n=20;%��H_20
Ka=zeros(1,n+1);
Lb=zeros(1,n+1);
Ka(1:10)=Kaa;
Lb(1:10)=Lbb;
u_Ve=4*Lb(1)/E-1;%��ֵ
%L0=Lb(1),K0=Ka(1);
data=Ve_diedai(u_Ve,n,Ve,VL_star,E,Ka,Lb,L);
format short e
%vt�仯-0.006 3.7 0.1 2.9
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
%vpa(A,5);
X1=linspace(Ve,Vt,110);
X2=linspace(Vt,VL_star,330);
X=[X1,X2];
[Ve_P0,Ve_P1,Ve_P2,Ve_P3,Ve_P4]=leijia2_Ve(u_Ve,X1,Ve,data,n);
[VL1_P0,VL1_P1,VL1_P2,VL1_P3,VL1_P4]=leijia2_VL(u_VL(1),X2,VL_star,data_1,n);
[VL2_P0,VL2_P1,VL2_P2,VL2_P3,VL2_P4]=leijia2_VL(u_VL(2),X2,VL_star,data_2,n);
[VL3_P0,VL3_P1,VL3_P2,VL3_P3,VL3_P4]=leijia2_VL(u_VL(3),X2,VL_star,data_3,n);
[VL4_P0,VL4_P1,VL4_P2,VL4_P3,VL4_P4]=leijia2_VL(u_VL(4),X2,VL_star,data_4,n);
P0=[K(1)*Ve_P0,K(2)*VL1_P0+K(3)*VL2_P0+K(4)*VL3_P0+K(5)*VL4_P0];
P1=[K(1)*Ve_P1,K(2)*VL1_P1+K(3)*VL2_P1+K(4)*VL3_P1+K(5)*VL4_P1];
P2=[K(1)*Ve_P2,K(2)*VL1_P2+K(3)*VL2_P2+K(4)*VL3_P2+K(5)*VL4_P2];
P3=[K(1)*Ve_P3,K(2)*VL1_P3+K(3)*VL2_P3+K(4)*VL3_P3+K(5)*VL4_P3];
P4=[K(1)*Ve_P4,K(2)*VL1_P4+K(3)*VL2_P4+K(4)*VL3_P4+K(5)*VL4_P4];
u_all=[u_Ve,u_VL];
Popen=P_open(u_all,Vt,Ve,VL_star,data,data_1,data_2,data_3,data_4,n,K);
Pq0=P_q(u_all,Vt,Ve,VL_star,data,data_1,data_2,data_3,data_4,n,K,1);
Pq1=P_q(u_all,Vt,Ve,VL_star,data,data_1,data_2,data_3,data_4,n,K,2);
Pq2=P_q(u_all,Vt,Ve,VL_star,data,data_1,data_2,data_3,data_4,n,K,3);
Pq3=P_q(u_all,Vt,Ve,VL_star,data,data_1,data_2,data_3,data_4,n,K,4);
%�򿪸���
P=P0+P1+P2+P3+P4;
pq=1-Popen;
subplot(2,3,1);
plot(X,P0)
text(Ve,0.2,'\downarrowVe');text(VL_star,0.2,'\downarrowVL')
xlabel('V');ylabel('P(0,V)');
subplot(2,3,2);
plot(X,P1)
text(Ve,0.2,'\downarrowVe');text(VL_star,0.2,'\downarrowVL')
xlabel('V');ylabel('P(1,V)');
subplot(2,3,3);
plot(X,P2)
text(Ve,0.2,'\downarrowVe');text(VL_star,0.2,'\downarrowVL')
xlabel('V');ylabel('P(2,V)');
subplot(2,3,4);
plot(X,P3)
text(Ve,0.2,'\downarrowVe');text(VL_star,0.2,'\downarrowVL')
xlabel('V');ylabel('P(3,V)');
subplot(2,3,5);
plot(X,P4)
text(Ve,0.0002,'\downarrowVe');text(VL_star,0.0002,'\downarrowVL')
xlabel('V');ylabel('P(4,V)');
subplot(2,3,6);
plot(X,P0+P1+P2+P3+P4)
text(Ve,0.2,'\downarrowVe');text(VL_star,0.2,'\downarrowVL')
xlabel('V');ylabel('P(V)');