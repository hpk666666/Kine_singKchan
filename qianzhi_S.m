syms x m  
M=1*10^-12;%����Ĥ�������λ��m^2%%%0.5
Iex=0.006*10^-12;%�ⲿ��������λ A%%%-0.006-----0.1
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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