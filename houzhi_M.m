Popen=zeros(1,n);%打开概率
Pclose=zeros(1,n);%关闭概率
PEDC=zeros(1,n);%通道关闭时候的泄露电流功耗
PEDO=zeros(1,n);%通道打开时候的泄露电流功耗
PED=zeros(1,n);%通道打开时候的泄露电流功耗
ExC=zeros(1,n);%均值
Ex2C=zeros(1,n);%二阶矩
ExO=zeros(1,n);%均值
Ex2O=zeros(1,n);%二阶矩
PC=cell(1,6);
PO=cell(1,6);
PD=cell(1,6);
EXY=zeros(1,n);%IEX V
COV=zeros(1,n);%IEX V方差
Ex=zeros(1,n);%均值
Ex2=zeros(1,n);%V二阶矩
Dx=zeros(1,n);%V方差
g=zeros(1,n);%
Var_g=zeros(1,n);%
ALLCOV=cell(1,5);
ALLVar_g=cell(1,5);
ALLDx=cell(1,5);
ALLEx=cell(1,5);
ALLEx_g=cell(1,5);
ALLPEDC=cell(1,5);
ALLPEDO=cell(1,5);
ALLPED=cell(1,5);
ALLPopen=cell(1,5);
ALLNED=cell(1,5);
NED=zeros(1,n);
for k=1:5
for i=1:n
gL=3*M(i);%单位是电导率 S/m ^2
ge=gL+gK;%单位是电导率 S
Cm=0.01*M(i);%单位是电容 F
E=ge/Cm;%单位是 pS/pF
L=gL/Cm;%单位是 pS/pF
[Pe_0,Pe_1,Pe_2,Pe_3,Pe_4]=Pe_Vt(ALLu_Ve{1,k}(i),ALLVt{1,k}(i),ALLVe{1,k}(i),ALLdata{1,k}{1,i},ddai);
[PL1_0,PL1_1,PL1_2,PL1_3,PL1_4]=PL_Vt(ALLu_VL{1,k}{1,i}(1),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_1{1,k}{1,i},ddai);
[PL2_0,PL2_1,PL2_2,PL2_3,PL2_4]=PL_Vt(ALLu_VL{1,k}{1,i}(2),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_2{1,k}{1,i},ddai);
[PL3_0,PL3_1,PL3_2,PL3_3,PL3_4]=PL_Vt(ALLu_VL{1,k}{1,i}(3),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_3{1,k}{1,i},ddai);
[PL4_0,PL4_1,PL4_2,PL4_3,PL4_4]=PL_Vt(ALLu_VL{1,k}{1,i}(4),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_4{1,k}{1,i},ddai);
Pe_int=Ve_int(ALLu_Ve{1,k}(i),ALLVt{1,k}(i),ALLVe{1,k}(i),ALLdata{1,k}{1,i},ddai);%这里利用牛顿莱布尼茨公式算(Pe(0,Vt)+...Pe(4,Vt))在【Ve,Vt】的积分
PL_int_1=VL_int(ALLu_VL{1,k}{1,i}(1),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_1{1,k}{1,i},ddai);%这里利用牛顿莱布尼茨公式算(PL1(0,Vt)+...PL1(4,Vt))在【Vt,VL】的积分
PL_int_2=VL_int(ALLu_VL{1,k}{1,i}(2),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_2{1,k}{1,i},ddai);
PL_int_3=VL_int(ALLu_VL{1,k}{1,i}(3),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_3{1,k}{1,i},ddai);
PL_int_4=VL_int(ALLu_VL{1,k}{1,i}(4),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_4{1,k}{1,i},ddai);
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
[E_eT_0,E_eT_1,E_eT_2,E_eT_3,E_eT_4]=Ee(ALLu_Ve{1,k}(i),ALLVt{1,k}(i),ALLVe{1,k}(i),ALLdata{1,k}{1,i},ddai);
[E_TL1_0,E_TL1_1,E_TL1_2,E_TL1_3,E_TL1_4]=EL(ALLu_VL{1,k}{1,i}(1),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_1{1,k}{1,i},ddai);
[E_TL2_0,E_TL2_1,E_TL2_2,E_TL2_3,E_TL2_4]=EL(ALLu_VL{1,k}{1,i}(2),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_2{1,k}{1,i},ddai);
[E_TL3_0,E_TL3_1,E_TL3_2,E_TL3_3,E_TL3_4]=EL(ALLu_VL{1,k}{1,i}(3),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_3{1,k}{1,i},ddai);
[E_TL4_0,E_TL4_1,E_TL4_2,E_TL4_3,E_TL4_4]=EL(ALLu_VL{1,k}{1,i}(4),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_4{1,k}{1,i},ddai);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%EX
Ex(i)=K(1)*(E_eT_0+E_eT_1+E_eT_2+E_eT_3+E_eT_4)+...
   K(2)*(E_TL1_0+E_TL1_1+E_TL1_2+E_TL1_3+E_TL1_4)+...
   K(3)*(E_TL2_0+E_TL2_1+E_TL2_2+E_TL2_3+E_TL2_4)+...
   K(4)*(E_TL3_0+E_TL3_1+E_TL3_2+E_TL3_3+E_TL3_4)+...
   K(5)*(E_TL4_0+E_TL4_1+E_TL4_2+E_TL4_3+E_TL4_4);
%E_TL2_1代表第二组数据1状态，VL到Vt的均值
%Ec(X)%%%%%%%%%%%%%%%%%%%
ExC(i)=K(1)*(E_eT_0+E_eT_1+E_eT_2+E_eT_3)+...
   K(2)*(E_TL1_0+E_TL1_1+E_TL1_2+E_TL1_3)+...
   K(3)*(E_TL2_0+E_TL2_1+E_TL2_2+E_TL2_3)+...
   K(4)*(E_TL3_0+E_TL3_1+E_TL3_2+E_TL3_3)+...
   K(5)*(E_TL4_0+E_TL4_1+E_TL4_2+E_TL4_3);
%%%%%%%%%%%%%%%%%
ExO(i)=K(1)*(E_eT_4)+...
   K(2)*(E_TL1_4)+...
   K(3)*(E_TL2_4)+...
   K(4)*(E_TL3_4)+...
   K(5)*(E_TL4_4);
%%%%%%%%%%D(X)%%%%%%%%%%%%%%%%%%%%%
[E2_eT_0,E2_eT_1,E2_eT_2,E2_eT_3,E2_eT_4]=Ee2(ALLu_Ve{1,k}(i),ALLVt{1,k}(i),ALLVe{1,k}(i),ALLdata{1,k}{1,i},ddai);
[E2_TL1_0,E2_TL1_1,E2_TL1_2,E2_TL1_3,E2_TL1_4]=EL2(ALLu_VL{1,k}{1,i}(1),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_1{1,k}{1,i},ddai);
[E2_TL2_0,E2_TL2_1,E2_TL2_2,E2_TL2_3,E2_TL2_4]=EL2(ALLu_VL{1,k}{1,i}(2),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_2{1,k}{1,i},ddai);
[E2_TL3_0,E2_TL3_1,E2_TL3_2,E2_TL3_3,E2_TL3_4]=EL2(ALLu_VL{1,k}{1,i}(3),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_3{1,k}{1,i},ddai);
[E2_TL4_0,E2_TL4_1,E2_TL4_2,E2_TL4_3,E2_TL4_4]=EL2(ALLu_VL{1,k}{1,i}(4),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_4{1,k}{1,i},ddai);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
Ex2(i)=Ex2C(i)+Ex2O(i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_all=[ALLu_Ve{1,k}(i),ALLu_VL{1,k}{1,i}];
Popen(i)=P_open(u_all,ALLVt{1,k}(i),ALLVe{1,k}(i),ALLVL_star{1,k}(i),ALLdata{1,k}{1,i},ALLdata_1{1,k}{1,i},ALLdata_2{1,k}{1,i},ALLdata_3{1,k}{1,i},ALLdata_4{1,k}{1,i},ddai,K);
Pclose(i)=1-Popen(i);
%%%%%%%%%%%%%%%%%%%%PED
PEDC(i)=gL*(Ex2C(i)-2*ALLVL_star{1,k}(i)*ExC(i)+ALLVL_star{1,k}(i)^2*Pclose(i));
PEDO(i)=ge*(Ex2O(i)-2*ALLVe{1,k}(i)*ExO(i)+ALLVe{1,k}(i)^2*Popen(i));
NED(i)=Popen(i)*(VK-ALLVL_star{1,k}(i))^2/(1/gK+1/gL);
PED(i)=PEDO(i)+PEDC(i)+NED(i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%COV%%%V AND Ni
EXY(i)=K(1)*E_eT_4+K(2)*E_TL1_4+K(3)*E_TL2_4+K(4)*E_TL3_4+K(5)*E_TL4_4;
COV(i)=EXY(i)-Popen(i)*Ex(i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%EX AND DX
Dx(i)=Ex2(i)-Ex(i)^2;
Ex;
%%%%%%%%%%%%%%%%%%%%%%%%%
g(i)=gL*Pclose(i)+ge*Popen(i);
Var_g(i)=ge^2*Popen(i)+gL^2*Pclose(i)-g(i)^2;
end
ALLCOV{1,k}=COV;%协方差
ALLVar_g{1,k}=Var_g;%电导的方差
ALLEx_g{1,k}=g;
ALLDx{1,k}=Dx;%电压的方差
ALLEx{1,k}=Ex;
ALLPEDC{1,k}=PEDC;%关闭时的扩散功耗
ALLPEDO{1,k}=PEDO;
ALLPED{1,k}=PED;
ALLPopen{1,k}=Popen;
ALLNED{1,k}=NED;
end
%plot(M,ALLPopen{1,k})