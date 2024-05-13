syms x m  
Iex=[-0.2,-0.1,0,0.1,0.2];
%linspace(-0.4,0.4,5);%1 A = 10^15 fA
% Iex=1;%外部电流，单位 fA 法安   1 A = 10^15 fA ，
n=100;
ddai=100;%求H_20迭代次数
gK=20;%单位是电导率 S
VL1=-54.4;%单位是 mV
VK=-77;%单位是 mV
K_b=1.380649*10^-23;%玻尔兹曼常数,单位为J/K，J为焦耳，K是热力学温度
T=294.15;%绝对温度，单位K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=linspace(0.02,0.08,n);%定义膜面积，单位是um^2
%iex=1,M=linspace(0.02,0.08,n)
ALLdata=cell(1,5);
data1=cell(1,n);
ALLP=cell(1,n);
ALLX=cell(1,n);
ALLK=cell(1,5);
Ve=zeros(1,n);
V_mid=zeros(1,n);
K=zeros(1,10);
K1=cell(1,n);
ALLV_mid=cell(1,5);
VL_star=zeros(1,n);
ALLVe=cell(1,5);
ALLVL_star=cell(1,5);
ALL_P=cell(1,5);
ALL_X=cell(1,5);
data=cell(1,10);%定义一行十列的元组一个里面是一个基解的幂级数系数
intP=cell(1,10);%归一化用到十个右边到左边零点基解的积分
P1_V0=cell(1,10);%元组的元胞是[P0,P1,P2,P3,P4]右边边界点的值
P2_V0=cell(1,10);%元组的元胞是[P0,P1,P2,P3,P4]左边边界点的值
P1_V0_2=cell(1,10);%元组的元胞是[P0,P1,P2,P3,P4]右边边界点导数的值
P2_V0_2=cell(1,10);%元组的元胞是[P0,P1,P2,P3,P4]左边边界点导数的值
ZY=eye(10);%自由未知量
o=10;%泰勒展开式展开次数
for k=1:5
for i=1:n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gL=3*M(i);%单位是电导率 S/m ^2
ge=gL+gK;%单位是电导率 S
Cm=0.01*M(i);%单位是电容 F
E=0.001*ge/Cm;%单位是 pS/pF
L=0.001*gL/Cm;%单位是 pS/pF
VL_star(i)=(gL*VL1+Iex(k))/gL;%单位是 mV %边界点
Ve(i)=(gL*VL1+gK*VK+Iex(k))/ge;%单位是 mV %边界点
V_mid(i)=Ve(i)+(VL_star(i)-Ve(i))*1/2;%展开点
delta_1=sqrt(2*K_b*T*gL*10^12*10^3)/Cm;%表示扩散的强度
delta_2=sqrt(2*K_b*T*ge*10^12*10^3)/Cm;%噪声误差1
fa=taylor(0.01*(x+55)/(1-exp(-(x+55)/10)),x,V_mid(i),'order',o);%速度α(v)在V0这个点的泰勒展式
fb=taylor(0.125*exp(-(x+65)/80),x,V_mid(i),'order',o);%速度β(v)在V0这个点的泰勒展式
newfa=subs(fa,(x-V_mid(i)),m);%将(x-V0)替换成m
newfb=subs(fb,(x-V_mid(i)),m);%将(x-V0)替换成m
Kaa=vpa(coeffs(newfa));%按升幂排列泰勒展式的系数
Lbb=vpa(coeffs(newfb));%按升幂排列泰勒展式的系数
Ka=zeros(1,ddai+1);
Lb=zeros(1,ddai+1);
Ka(1:o)=Kaa;
Lb(1:o)=Lbb;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:10
data{1,j}=V_diedai_noise(ddai,Ve(i),VL_star(i),V_mid(i),E,Ka,Lb,L,ZY(1:5,j),ZY(6:10,j),delta_1,delta_2);
%十个基础解              
intP{1,j}=V_int_noise(VL_star(i),Ve(i),V_mid(i),data{1,j},ddai);
%归一化用到的
P1_V0{1,j} = P_zeros_noise(VL_star(i),V_mid(i),data{1,j},ddai);%右边边界点的值
P2_V0{1,j} = P_zeros_noise(Ve(i),V_mid(i),data{1,j},ddai);%左边边界点的值
P1_V0_2{1,j} = P_zeros_noise_2(VL_star(i),V_mid(i),data{1,j},ddai);%右边导数的零点
P2_V0_2{1,j} = P_zeros_noise_2(Ve(i),V_mid(i),data{1,j},ddai);%左边导数的零点
end
data1{1,i}=data;
A=zeros(11,10);
for ii=1:11
    for jj=1:10
        if(ii==1)
            A(ii,jj)=P2_V0_2{1,jj}(5,1);%在Ve处n4状态的导数为0
        elseif(ii<6)
             A(ii,jj)=P2_V0_2{1,jj}(ii-1,1)-L*(VL_star(i)-Ve(i))*P2_V0{1,jj}(ii-1,1)*2/delta_1^2;
             %在Ve处n0 1 2 3状态的Pi导数等于L*(VL-Ve)*(2/delta_1^2)*Pi
         elseif(ii<10)
            A(ii,jj)=P1_V0_2{1,jj}(ii-5,1);  %在VL处n0 1 2 3 状态的导数为0
          elseif(ii==10)
             A(ii,jj)=P1_V0_2{1,jj}(5,1)-E*(Ve(i)-VL_star(i))*(2/delta_2^2)*P1_V0{1,jj}(5,1);
             %在VL处n4状态的Pi导数等于E*(Ve-VL)*(2/delta_2^2)*P4
              else
               A(ii,jj)=intP{1,jj}; %归一化
        end
    end
end
B=zeros(11,1);
B(11,1)=1;
K=A\B;
K1{1,i}=K;
X=linspace(Ve(i),VL_star(i),1000);
X_P=cell(1,10);
P=cell(1,5);
for jjj=1:10
    X_P{1,jjj}=leijia2_V_noise(X,V_mid(i),data{1,jjj},ddai);
end
for iii=1:5    
P{1,iii}=K(1)*X_P{1,1}{1,iii}+K(2)*X_P{1,2}{1,iii}+K(3)*X_P{1,3}{1,iii}+K(4)*X_P{1,4}{1,iii}...
        +K(5)*X_P{1,5}{1,iii}+K(6)*X_P{1,6}{1,iii}+K(7)*X_P{1,7}{1,iii}+K(8)*X_P{1,8}{1,iii}...
        +K(9)*X_P{1,9}{1,iii}+K(10)*X_P{1,10}{1,iii};
end
ALLP{1,i}=P;
ALLX{1,i}=X;
end
ALLK{1,k}=K1;%K=zeros(10,n);ALLK{1,k}(1到10,n)
ALLdata{1,k}=data1;%ALLdata{1,k}{1,1到n}{1,1到10}
ALLV_mid{1,k}=V_mid;
ALLVe{1,k}=Ve;
ALLVL_star{1,k}=VL_star;
ALL_P{1,k}=ALLP;
ALL_X{1,k}=ALLX;
end