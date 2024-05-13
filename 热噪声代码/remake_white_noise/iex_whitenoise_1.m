syms x m
M=0.02;%定义膜面积，单位是um^2%其实范围 0.15 
Iex=0 ;%外部电流，单位 fA 法安   1 A = 10^15 fA ，%其实范围 1
%M0.02到0.08|Iex-1到1fA 
%求H_20
ddai=210;
%iex=0,10,-10,确定膜面积
%10^-3是关于时间的换算系数与后面L和E处的系数一致
%电流等于0的时候，M面积的范围
gL=3*M;%单位是电导率 pS
gK=20;%单位是电导率 pS
VL1=-54.4;%单位是 mV
VK=-77;%单位是 mV
ge=gL+gK;%单位是电导率 pS
Cm=0.01*M;%单位是电容 pF
E=0.001*ge/Cm;%单位是 pS/pF
L=0.001*gL/Cm;%单位是 pS/pF
VL=(gL*VL1+Iex)/gL;%单位是 mV
Ve=(gL*VL1+gK*VK+Iex)/ge;%单位是 mV
V_mid=(VL-Ve)/2+Ve;%展开点
Ve_0=Ve;%边界点
VL_0=VL;%边界点
K_b=1.380649*10^-23;%玻尔兹曼常数,单位为J/K，J为焦耳，K是热力学温度
T=294.15;%绝对温度，单位K
delta_1=sqrt(2*K_b*T*gL*10^12*10^3)/Cm;%表示扩散的强度
delta_2=sqrt(2*K_b*T*ge*10^12*10^3)/Cm;%噪声误差1
o=10;
fa=taylor(0.01*(x+55)/(1-exp(-(x+55)/10)),x,V_mid,'order',o);%速度α(v)在V0这个点的泰勒展式
fb=taylor(0.125*exp(-(x+65)/80),x,V_mid,'order',o);%速度β(v)在V0这个点的泰勒展式
newfa=subs(fa,(x-V_mid),m);%将(x-V0)替换成m
newfb=subs(fb,(x-V_mid),m);%将(x-V0)替换成m
Kaa=vpa(coeffs(newfa));%按升幂排列泰勒展式的系数
Lbb=vpa(coeffs(newfb));%按升幂排列泰勒展式的系数
Ka=zeros(1,ddai+1);
Lb=zeros(1,ddai+1);
Ka(1:o)=Kaa;
Lb(1:o)=Lbb;
data=cell(1,10);%定义一行十列的元组一个里面是一个基解的幂级数系数
intP=cell(1,10);%归一化用到十个右边到左边零点基解的积分
P1_V0=cell(1,10);%元组的元胞是[P0,P1,P2,P3,P4]右边边界点的值
P2_V0=cell(1,10);%元组的元胞是[P0,P1,P2,P3,P4]左边边界点的值
P1_V0_2=cell(1,10);%元组的元胞是[P0,P1,P2,P3,P4]右边边界点导数的值
P2_V0_2=cell(1,10);%元组的元胞是[P0,P1,P2,P3,P4]左边边界点导数的值
ZY=eye(10);%自由未知量
for i=1:10
data{1,i}=V_diedai_noise(ddai,Ve,VL,V_mid,E,Ka,Lb,L,ZY(1:5,i),ZY(6:10,i),delta_1,delta_2);
%十个基础解              
intP{1,i}=V_int_noise(VL_0,Ve_0,V_mid,data{1,i},ddai);
%归一化用到的
P1_V0{1,i} = P_zeros_noise(VL_0,V_mid,data{1,i},ddai);%右边边界点的值
P2_V0{1,i} = P_zeros_noise(Ve_0,V_mid,data{1,i},ddai);%左边边界点的值
P1_V0_2{1,i} = P_zeros_noise_2(VL_0,V_mid,data{1,i},ddai);%右边导数的零点
P2_V0_2{1,i} = P_zeros_noise_2(Ve_0,V_mid,data{1,i},ddai);%左边导数的零点
end
%取-50mv和-70mv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=zeros(11,10);
for i=1:11
    for j=1:10
        if(i==1)
            A(i,j)=P2_V0_2{1,j}(5,1);%在Ve处n4状态的导数为0
        elseif(i<6)
             A(i,j)=P2_V0_2{1,j}(i-1,1)-L*(VL-Ve)*P2_V0{1,j}(i-1,1)*2/delta_1^2;
             %在Ve处n0 1 2 3状态的Pi导数等于L*(VL-Ve)*(2/delta_1^2)*Pi
         elseif(i<10)
            A(i,j)=P1_V0_2{1,j}(i-5,1);  %在VL处n0 1 2 3 状态的导数为0
          elseif(i==10)
             A(i,j)=P1_V0_2{1,j}(5,1)-E*(Ve-VL)*(2/delta_2^2)*P1_V0{1,j}(5,1);
             %在VL处n4状态的Pi导数等于E*(Ve-VL)*(2/delta_2^2)*P4
              else
               A(i,j)=intP{1,j}; %归一化
        end
    end
end
% Aa=size(A,1);
% while rank(A)~=Aa %当构造出的H不是行满秩时，找出它行向量组的极大无关组组成newH代替H
% [R,jb]=rref(A'); %结果中R为rref函数中目标矩阵的行阶梯最简矩阵，jb是一个向量，为目标矩阵的列极大无关组所在的列数。注意这里需要先将H转置，因为我们的目标是要求H的行向量组的极大无关组
% newH=A(jb,:);
% A=newH;
% Aa=size(A,1);
% end
B=zeros(11,1);
B(11,1)=1;
K=A\B;
X=linspace(Ve_0,VL_0,1000);
X_P=cell(1,10);
P=cell(1,5);
for j=1:10
    X_P{1,j}=leijia2_V_noise(X,V_mid,data{1,j},ddai);
end
for i=1:5    
P{1,i}=K(1)*X_P{1,1}{1,i}+K(2)*X_P{1,2}{1,i}+K(3)*X_P{1,3}{1,i}+K(4)*X_P{1,4}{1,i}...
        +K(5)*X_P{1,5}{1,i}+K(6)*X_P{1,6}{1,i}+K(7)*X_P{1,7}{1,i}+K(8)*X_P{1,8}{1,i}...
        +K(9)*X_P{1,9}{1,i}+K(10)*X_P{1,10}{1,i};
end
Pp=P{1,5}+P{1,4}+P{1,3}+P{1,2}+P{1,1};
%打开概率
P0=trapz(X,P{1,1});
P1=trapz(X,P{1,2});
P2=trapz(X,P{1,3});
P3=trapz(X,P{1,4});
P4=trapz(X,P{1,5});
Pall=P0+P1+P2+P3+P4;
P0=P0/Pall;
P1=P1/Pall;
P2=P2/Pall;
P3=P3/Pall;
P4=P4/Pall;
figure(1)
subplot(2,3,1);
plot(X,P{1,1})
xlabel('V');ylabel('P(0,V)');
set(gca,'xlim',[Ve_0,VL_0]);%x坐标轴范围
subplot(2,3,2);
plot(X,P{1,2})
xlabel('V');ylabel('P(1,V)');
set(gca,'xlim',[Ve_0,VL_0]);%x坐标轴范围
% set(gca,'ylim',[0,0.025]);%y坐标轴范围
subplot(2,3,3);
plot(X,P{1,3})
xlabel('V');ylabel('P(2,V)');
set(gca,'xlim',[Ve_0,VL_0]);%x坐标轴范围
% set(gca,'ylim',[0,0.035]);%y坐标轴范围
subplot(2,3,4);
plot(X,P{1,4})
xlabel('V');ylabel('P(3,V)');
set(gca,'xlim',[Ve_0,VL_0]);%x坐标轴范围
% set(gca,'ylim',[0,0.016]);%y坐标轴范围
subplot(2,3,5);
plot(X,P{1,5})
xlabel('V');ylabel('P(4,V)');
set(gca,'xlim',[Ve_0,VL_0]);%x坐标轴范围
subplot(2,3,6);
plot(X,P{1,5}+P{1,4}+P{1,3}+P{1,2}+P{1,1})
%text(Vt,0.2,'\downarrowVt');text(Ve,0.2,'\downarrowVe');text(VL,0.2,'\downarrowVL')
xlabel('V');ylabel('P(V)');
set(gca,'xlim',[Ve_0,VL_0]);%x坐标轴范围
% set(gca,'ylim',[0,0.08]);%y坐标轴范围