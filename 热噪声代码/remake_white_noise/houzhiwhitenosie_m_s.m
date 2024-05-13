Popen=zeros(1,n);%打开概率
Pclose=zeros(1,n);%关闭概率
PEDC=zeros(1,n);%通道关闭时候的泄露电流功耗
PEDO=zeros(1,n);%通道打开时候的泄露电流功耗
PED=zeros(1,n);%通道打开时候的泄露电流功耗
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
ALL_D=cell(1,5);
ALL_S=cell(1,5);
D=zeros(1,n);
S=zeros(1,n);
for k=1:5
Ex_1=cell(1,10);
Ex_2=cell(1,10);
ExC=zeros(1,n);%均值
Ex2C=zeros(1,n);%二阶矩
ExO=zeros(1,n);%均值
Ex2O=zeros(1,n);%二阶矩
EXY=zeros(1,n);%IEX V
for i=1:n
gL=3*M(i);%单位是电导率 S/m ^2
ge=gL+gK;%单位是电导率 S
Cm=0.01*M(i);%单位是电容 F
for ii=1:10 %10个基础解
%%%%%%%%%%E(X)%%%%%%%%%%%%%%%%%%%%%
Ex_1{1,ii}=EX_1(ALLVL_star{1,k}(i),ALLV_mid{1,k}(i),ALLVe{1,k}(i),ALLdata{1,k}{1,i}{1,ii},ddai);%ALLdata{1,k}{1,1到n}{1,1到10}
%%%%%%%%%%D(X)%%%%%%%%%%%%%%%%%%%%%ALLK{1,k}(1到10,n)ALLK{1,k}(1到10,n)
Ex_2{1,ii}=EX_2(ALLVL_star{1,k}(i),ALLV_mid{1,k}(i),ALLVe{1,k}(i),ALLdata{1,k}{1,i}{1,ii},ddai);
ExC(i)=ALLK{1,k}{1,i}(ii)*(Ex_1{1,ii}(1,1)+Ex_1{1,ii}(1,2)+Ex_1{1,ii}(1,3)+Ex_1{1,ii}(1,4))+ExC(i);
ExO(i)=ALLK{1,k}{1,i}(ii)*(Ex_1{1,ii}(1,5))+ExO(i);
%%%%%%%%%%%
Ex2C(i)=ALLK{1,k}{1,i}(ii)*(Ex_2{1,ii}(1,1)+Ex_2{1,ii}(1,2)+Ex_2{1,ii}(1,3)+Ex_2{1,ii}(1,4))+Ex2C(i);
Ex2O(i)=ALLK{1,k}{1,i}(ii)*(Ex_2{1,ii}(1,5))+Ex2O(i);
%%%%%%%%%%%
EXY(i)=ALLK{1,k}{1,i}(ii)*(Ex_1{1,ii}(1,5))+EXY(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%EX AND EX2
Ex(i)=ExC(i)+ExO(i);
Ex2(i)=Ex2C(i)+Ex2O(i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Popen(i)=P_open_noise(ALLVL_star{1,k}(i),ALLVe{1,k}(i),ALLV_mid{1,k}(i),ALLdata{1,k}{1,i},ddai,ALLK{1,k}{1,i});
Pclose(i)=1-Popen(i);
%%%%%%%%%%%%%%%%%%%%PED
PEDC(i)=gL*(Ex2C(i)-2*ALLVL_star{1,k}(i)*ExC(i)+ALLVL_star{1,k}(i)^2*Pclose(i));
PEDO(i)=ge*(Ex2O(i)-2*ALLVe{1,k}(i)*ExO(i)+ALLVe{1,k}(i)^2*Popen(i));
NED(i)=Popen(i)*(VK-VL1)^2/(1/gK+1/gL);
PED(i)=PEDO(i)+PEDC(i)+NED(i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%COV%%%V AND Ni
Dx(i)=Ex2(i)-Ex(i)^2;
%%%%%%%%%%%%%%%%%%%%%%%%%
COV(i)=EXY(i)-Popen(i)*Ex(i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g(i)=gL*Pclose(i)+ge*Popen(i);
Var_g(i)=ge^2*Popen(i)+gL^2*Pclose(i)-g(i)^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%熵，相对熵

P0=ALL_P{1,k}{1,i}{1,1};
P1=ALL_P{1,k}{1,i}{1,2};
P2=ALL_P{1,k}{1,i}{1,3};
P3=ALL_P{1,k}{1,i}{1,4};
P4=ALL_P{1,k}{1,i}{1,5};
X=ALL_X{1,k}{1,i};
%%%%%%%%%%%%%%%%%%%%%%%%%%%
dP0=P0(1:length(P0));
dP1=P1(1:length(P0));
dP2=P2(1:length(P0));
dP3=P3(1:length(P0));
dP4=P4(1:length(P0));
dX=X(1:length(P0));
[nP0,nP1,nP2,nP3,nP4]=newP(dP0,dP1,dP2,dP3,dP4);
S0=trapz(dX,nP0);
S1=trapz(dX,nP1);
S2=trapz(dX,nP2);
S3=trapz(dX,nP3);
S4=trapz(dX,nP4);
%%%%%%%%%%%%%%%%%%%%%%%%熵
S(i)=-(S0+S1+S2+S3+S4);
%%%%%%%%%%%%%%%%%%%%%%%%散度
[Pi]=DPi(dP0,dP1,dP2,dP3,dP4,dX);%?
[DP0,DP1,DP2,DP3,DP4]=DnewP(dP0,dP1,dP2,dP3,dP4,Pi);%?
D0=trapz(dX,DP0);
D1=trapz(dX,DP1);
D2=trapz(dX,DP2);
D3=trapz(dX,DP3);
D4=trapz(dX,DP4);
% [D0,D1,D2,D3,D4]=jifen(DP1,DP2,DP3,DP4,DP5,dX);
%%%%%%%%%%%%%%%%%%%%%%%%相对熵
D(i)=(D0+D1+D2+D3+D4);
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
ALL_D{1,k}=D;
ALL_S{1,k}=S;
end

% for k=1:5
% for i=1:n
% P0=ALL_P{1,k}{1,i}{1,1};
% P1=ALL_P{1,k}{1,i}{1,2};
% P2=ALL_P{1,k}{1,i}{1,3};
% P3=ALL_P{1,k}{1,i}{1,4};
% P4=ALL_P{1,k}{1,i}{1,5};
% X=ALL_X{1,k}{1,i};
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% dP0=P0(2:length(P0)-1);
% dP1=P1(2:length(P0)-1);
% dP2=P2(2:length(P0)-1);
% dP3=P3(2:length(P0)-1);
% dP4=P4(2:length(P0)-1);
% dX=X(2:length(P0)-1);
% [nP0,nP1,nP2,nP3,nP4]=newP(dP0,dP1,dP2,dP3,dP4);
% S0=trapz(dX,nP0);
% S1=trapz(dX,nP1);
% S2=trapz(dX,nP2);
% S3=trapz(dX,nP3);
% S4=trapz(dX,nP4);
% S(i)=-(S0+S1+S2+S3+S4);
% [Pi]=DPi(dP0,dP1,dP2,dP3,dP4,dX);%?
% [DP0,DP1,DP2,DP3,DP4]=DnewP(dP0,dP1,dP2,dP3,dP4,Pi);%?
% D0=trapz(dX,DP0);
% D1=trapz(dX,DP1);
% D2=trapz(dX,DP2);
% D3=trapz(dX,DP3);
% D4=trapz(dX,DP4);
% % [D0,D1,D2,D3,D4]=jifen(DP1,DP2,DP3,DP4,DP5,dX);
% D(i)=-(D0+D1+D2+D3+D4);
% end
% ALL_D{1,k}=D;
% ALL_S{1,k}=S;
% end