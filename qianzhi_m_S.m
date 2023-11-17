syms x m  
Iex=linspace(-0.006*10^-12,0.006*10^-12,5);
n=100;%求H_20
gK=20*10^-12;%单位是电导率 S
VL1=-54.4*10^-3;%单位是 V
VK=-77*10^-3;%单位是 V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=linspace(0.05*10^-12,3*10^-12,n);%定义膜面积，单位是m^2
data_1=cell(1,n);
data_2=cell(1,n);
data_3=cell(1,n);
data_4=cell(1,n);
data=cell(1,n);
ALLdata_1=cell(1,5);
ALLdata_2=cell(1,5);
ALLdata_3=cell(1,5);
ALLdata_4=cell(1,5);
ALLdata=cell(1,5);
u_VL=cell(1,n);
u_Ve=zeros(1,n);
Ve=zeros(1,n);
Vt=zeros(1,n);
ALLVt=cell(1,5);
ALLu_VL=cell(1,5);
ALLu_Ve=cell(1,5);
VL_star=zeros(1,n);
ALLVe=cell(1,5);
ALLVL_star=cell(1,5);
for k=1:5
for i=1:n
gL=3*M(i);%单位是电导率 S/m ^2
ge=gL+gK;%单位是电导率 S
Cm=0.01*M(i);%单位是电容 F
E=ge/Cm;%单位是 pS/pF
L=gL/Cm;%单位是 pS/pF
VL_star(i)=(gL*VL1+Iex(k))/gL;%单位是 V
Ve(i)=(gL*VL1+gK*VK+Iex(k))/ge;%单位是 V
Vt(i)=Ve(i)+(VL_star(i)-Ve(i))*2.9/4;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
faL=taylor(0.01*10^6*(x+55*10^-3)/(1-exp(-(x+55*10^-3)/0.01)),x,VL_star(i),'order',10);%速度α(v)在Ve这个点的泰勒展式
fbL=taylor(0.125*10^3*exp(-(x+65*10^-3)/0.08),x,VL_star(i),'order',10);%速度β(v)在Ve这个点的泰勒展式
newfaL=subs(faL,(x-VL_star(i)),m);%将(x-VL)替换成m
newfbL=subs(fbL,(x-VL_star(i)),m);%将(x-VL)替换成m
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
u_VL{1,i}=[s1(2),s1(3),s1(4),s1(5)];%除去第一个-1的u值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fae=taylor(0.01*10^6*(x+55*10^-3)/(1-exp(-(x+55*10^-3)/0.01)),x,Ve(i),'order',10);%速度α(v)在Ve这个点的泰勒展式
fbe=taylor(0.125*10^3*exp(-(x+65*10^-3)/0.08),x,Ve(i),'order',10);%速度β(v)在Ve这个点的泰勒展式
newfa=subs(fae,(x-Ve(i)),m);%将(x-Ve)替换成m
newfb=subs(fbe,(x-Ve(i)),m);%将(x-Ve)替换成m
Kae=vpa(coeffs(newfa));%按升幂排列泰勒展式的系数
Lbe=vpa(coeffs(newfb));%按升幂排列泰勒展式的系数
Ka_e=zeros(1,n+1);
Lb_e=zeros(1,n+1);
Ka_e(1:10)=Kae;
Lb_e(1:10)=Lbe;
u_Ve(i)=4*Lb_e(1)/E-1;%μ值,会随膜面积变化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short e
data_1{1,i}=VL_diedai(1,u_VL{1,i},n,Ve(i),VL_star(i),E,Ka_L,Lb_L,L);%1,2,3,4分别代表四个不同的u值的pi
data_2{1,i}=VL_diedai(2,u_VL{1,i},n,Ve(i),VL_star(i),E,Ka_L,Lb_L,L);
data_3{1,i}=VL_diedai(3,u_VL{1,i},n,Ve(i),VL_star(i),E,Ka_L,Lb_L,L);
data_4{1,i}=VL_diedai(4,u_VL{1,i},n,Ve(i),VL_star(i),E,Ka_L,Lb_L,L);
data{1,i}=Ve_diedai(u_Ve(i),n,Ve(i),VL_star(i),E,Ka_e,Lb_e,L);
end
ALLdata_1{1,k}=data_1;
ALLdata_2{1,k}=data_2;
ALLdata_3{1,k}=data_3;
ALLdata_4{1,k}=data_4;
ALLdata{1,k}=data;
ALLVt{1,k}=Vt;
ALLVe{1,k}=Ve;
ALLVL_star{1,k}=VL_star;
ALLu_VL{1,k}=u_VL;
ALLu_Ve{1,k}=u_Ve;
end