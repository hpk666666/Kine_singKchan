syms x m  
Iex=[-0.2,-0.1,0,0.1,0.2];
%linspace(-0.4,0.4,5);%1 A = 10^15 fA
% Iex=1;%�ⲿ��������λ fA ����   1 A = 10^15 fA ��
n=100;
ddai=100;%��H_20��������
gK=20;%��λ�ǵ絼�� S
VL1=-54.4;%��λ�� mV
VK=-77;%��λ�� mV
K_b=1.380649*10^-23;%������������,��λΪJ/K��JΪ������K������ѧ�¶�
T=294.15;%�����¶ȣ���λK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=linspace(0.02,0.08,n);%����Ĥ�������λ��um^2
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
data=cell(1,10);%����һ��ʮ�е�Ԫ��һ��������һ��������ݼ���ϵ��
intP=cell(1,10);%��һ���õ�ʮ���ұߵ����������Ļ���
P1_V0=cell(1,10);%Ԫ���Ԫ����[P0,P1,P2,P3,P4]�ұ߽߱���ֵ
P2_V0=cell(1,10);%Ԫ���Ԫ����[P0,P1,P2,P3,P4]��߽߱���ֵ
P1_V0_2=cell(1,10);%Ԫ���Ԫ����[P0,P1,P2,P3,P4]�ұ߽߱�㵼����ֵ
P2_V0_2=cell(1,10);%Ԫ���Ԫ����[P0,P1,P2,P3,P4]��߽߱�㵼����ֵ
ZY=eye(10);%����δ֪��
o=10;%̩��չ��ʽչ������
for k=1:5
for i=1:n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gL=3*M(i);%��λ�ǵ絼�� S/m ^2
ge=gL+gK;%��λ�ǵ絼�� S
Cm=0.01*M(i);%��λ�ǵ��� F
E=0.001*ge/Cm;%��λ�� pS/pF
L=0.001*gL/Cm;%��λ�� pS/pF
VL_star(i)=(gL*VL1+Iex(k))/gL;%��λ�� mV %�߽��
Ve(i)=(gL*VL1+gK*VK+Iex(k))/ge;%��λ�� mV %�߽��
V_mid(i)=Ve(i)+(VL_star(i)-Ve(i))*1/2;%չ����
delta_1=sqrt(2*K_b*T*gL*10^12*10^3)/Cm;%��ʾ��ɢ��ǿ��
delta_2=sqrt(2*K_b*T*ge*10^12*10^3)/Cm;%�������1
fa=taylor(0.01*(x+55)/(1-exp(-(x+55)/10)),x,V_mid(i),'order',o);%�ٶȦ�(v)��V0������̩��չʽ
fb=taylor(0.125*exp(-(x+65)/80),x,V_mid(i),'order',o);%�ٶȦ�(v)��V0������̩��չʽ
newfa=subs(fa,(x-V_mid(i)),m);%��(x-V0)�滻��m
newfb=subs(fb,(x-V_mid(i)),m);%��(x-V0)�滻��m
Kaa=vpa(coeffs(newfa));%����������̩��չʽ��ϵ��
Lbb=vpa(coeffs(newfb));%����������̩��չʽ��ϵ��
Ka=zeros(1,ddai+1);
Lb=zeros(1,ddai+1);
Ka(1:o)=Kaa;
Lb(1:o)=Lbb;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:10
data{1,j}=V_diedai_noise(ddai,Ve(i),VL_star(i),V_mid(i),E,Ka,Lb,L,ZY(1:5,j),ZY(6:10,j),delta_1,delta_2);
%ʮ��������              
intP{1,j}=V_int_noise(VL_star(i),Ve(i),V_mid(i),data{1,j},ddai);
%��һ���õ���
P1_V0{1,j} = P_zeros_noise(VL_star(i),V_mid(i),data{1,j},ddai);%�ұ߽߱���ֵ
P2_V0{1,j} = P_zeros_noise(Ve(i),V_mid(i),data{1,j},ddai);%��߽߱���ֵ
P1_V0_2{1,j} = P_zeros_noise_2(VL_star(i),V_mid(i),data{1,j},ddai);%�ұߵ��������
P2_V0_2{1,j} = P_zeros_noise_2(Ve(i),V_mid(i),data{1,j},ddai);%��ߵ��������
end
data1{1,i}=data;
A=zeros(11,10);
for ii=1:11
    for jj=1:10
        if(ii==1)
            A(ii,jj)=P2_V0_2{1,jj}(5,1);%��Ve��n4״̬�ĵ���Ϊ0
        elseif(ii<6)
             A(ii,jj)=P2_V0_2{1,jj}(ii-1,1)-L*(VL_star(i)-Ve(i))*P2_V0{1,jj}(ii-1,1)*2/delta_1^2;
             %��Ve��n0 1 2 3״̬��Pi��������L*(VL-Ve)*(2/delta_1^2)*Pi
         elseif(ii<10)
            A(ii,jj)=P1_V0_2{1,jj}(ii-5,1);  %��VL��n0 1 2 3 ״̬�ĵ���Ϊ0
          elseif(ii==10)
             A(ii,jj)=P1_V0_2{1,jj}(5,1)-E*(Ve(i)-VL_star(i))*(2/delta_2^2)*P1_V0{1,jj}(5,1);
             %��VL��n4״̬��Pi��������E*(Ve-VL)*(2/delta_2^2)*P4
              else
               A(ii,jj)=intP{1,jj}; %��һ��
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
ALLK{1,k}=K1;%K=zeros(10,n);ALLK{1,k}(1��10,n)
ALLdata{1,k}=data1;%ALLdata{1,k}{1,1��n}{1,1��10}
ALLV_mid{1,k}=V_mid;
ALLVe{1,k}=Ve;
ALLVL_star{1,k}=VL_star;
ALL_P{1,k}=ALLP;
ALL_X{1,k}=ALLX;
end