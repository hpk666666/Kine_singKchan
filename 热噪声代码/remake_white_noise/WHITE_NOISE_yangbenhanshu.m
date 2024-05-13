tic;
clear;
M_area=0.05;%定义膜面积，单位是um^2
Iex=0;%外部电流，单位 fA 法安   1 A = 10^15 fA ，
Cm=0.01*M_area;%单位是电容 pF
gL=3*M_area;%单位是电导率 pS
gK=20;%单位是电导率 pS
V_L=-54.4;%单位是 mV
V_K=-77;%单位是 mV
ge=gL+gK;%单位是电导率 pS
V_e=(gL*V_L+gK*V_K+Iex)/ge;%单位是 V
V_L=(gL*V_L+Iex)/gL;%单位是 V
Cond_L=0.001*gL/Cm;
Cond_k=0.001*ge/Cm;
%%%%%%%%%%%%%%%%%%%%%%%%%%
K_b=1.380649*10^-23;%玻尔兹曼常数,单位为J/K，J为焦耳，K是热力学温度
T=294.15;%绝对温度，单位K
delta_1=sqrt(2*K_b*T*gL*10^12*10^3)/Cm;%表示扩散的强度
delta_2=sqrt(2*K_b*T*ge*10^12*10^3)/Cm;%噪声误差1
%Alpha_V=0.01*(V+55)/(1-exp(-(V+55)/10));
%Beta_V=0.125*exp(-(V+65)/80);
%K01=4*Alpha_V;
%K10=Beta_V;K12=3*Alpha_V;
%K21=2*Beta_V;K23=2*Alpha_V;
%K32=3*Beta_V;K34=Alpha_V;
%K43=4*Beta_V;
%转移速率
%lambda_t=[K01,K10+K12,K21+K23,K32+K34,K43];
%lambda_tt=0;
%lambda_tt=[K10/(K10+K12),K21/(K21+K23),K32/(K32+K34)];

delta_t=5*10^-6;%单位ms，这里是deltaB随机力的方差deltaB~（0，delta_t）
T_step=2*10^7;
%由于随机力的原因，这里用的莫得卡洛模拟。但是因为随机性的原因。随机出来的数个别可能
%会特别大，可能加一个限制条件比较好。设置膜电压范围为[-80，-30]
V=zeros(1,T_step);
S_k=zeros(1,T_step);

%这里的参数值要注意一下.  这里定的是:C=1

V(1)=-55;
n_i=1;

while n_i<T_step
    
   switch S_k(n_i)
      %-------------------------  
              case 0
                  %calculate time that the channel stay at "0"   lambda=4*alfa_n; 
                   U1=rand;
                   F=1;
                   integ_lambda=0;        %the integrel of lambda(s)
                   while F>U1  && n_i<T_step
                       U3=normrnd(0,sqrt(delta_t),1);
                       V(n_i+1)=V(n_i)+(-Cond_L*(V(n_i)-V_L))*delta_t+delta_1*U3;  
                           if V(n_i)>V_L
                             V(n_i)=(V_L-V(n_i))+V_L;            
                           end
      %-------------------------  
                           if V(n_i)<V_e
                                  V(n_i)=(V_e-V(n_i))+V_e;              
                           end
                       S_k(n_i+1)=S_k(n_i);
                       alfa_n=0.01*(V(n_i)+55)/(1-exp(-(V(n_i)+55)/10));
                       lambda=4*alfa_n;
                       integ_lambda=integ_lambda+lambda*delta_t;
                       F=exp(-integ_lambda);
                       n_i=n_i+1;
                   end
            
           
                 %at the jump time , system jump to "1"; definitely
            
                 S_k(n_i)=1;
         
            
              %-------------------------  
              case 1
                  %calculate time that the channel stay at "1"   lambda=3*alfa_n+1*bata_n; 
                   U1=rand;
                   F=1;
                   integ_lambda=0;        %the integrel of lambda(s)
                   while F>U1   &&   n_i<T_step
                       U3=normrnd(0,sqrt(delta_t),1);
                       V(n_i+1)=V(n_i)+(-Cond_L*(V(n_i)-V_L))*delta_t+delta_1*U3;
                           if V(n_i)>V_L
                             V(n_i)=(V_L-V(n_i))+V_L;            
                           end
      %-------------------------  
                           if V(n_i)<V_e
                                  V(n_i)=(V_e-V(n_i))+V_e;              
                           end
                       S_k(n_i+1)=S_k(n_i);
                       alfa_n=0.01*(V(n_i)+55)/(1-exp(-(V(n_i)+55)/10));
                       bata_n=0.125*exp(-(V(n_i)+65)/80); 
                       lambda=3*alfa_n+bata_n;
                       integ_lambda=integ_lambda+lambda*delta_t;
                       F=exp(-integ_lambda);
                       n_i=n_i+1;
                   end
           
                 
                 %at the jump time , system jump to "0" or "2"; with
                 %probability ...rate            
            
                U2=rand;
                if U2<=3*alfa_n/lambda; 
                    S_k(n_i)=2;   %jump to state "2"
                else
                   S_k(n_i)=0;   %jump to state "0"
                 end 
            
            
          %-------------------------  
              case 2
                   %calculate time that the channel stay at "2"   lambda=2*alfa_n+2*bata_n; 
                   U1=rand;
                   F=1;
                   integ_lambda=0;        %the integrel of lambda(s)
                   while F>U1  && n_i<T_step
                       U3=normrnd(0,sqrt(delta_t),1);
                       V(n_i+1)=V(n_i)+(-Cond_L*(V(n_i)-V_L))*delta_t+delta_1*U3;
                        if V(n_i)>V_L
                             V(n_i)=(V_L-V(n_i))+V_L;            
                        end
      %-------------------------  
                           if V(n_i)<V_e
                                  V(n_i)=(V_e-V(n_i))+V_e;              
                           end
                       S_k(n_i+1)=S_k(n_i);
                       alfa_n=0.01*(V(n_i)+55)/(1-exp(-(V(n_i)+55)/10));
                       bata_n=0.125*exp(-(V(n_i)+65)/80); 
                       lambda=2*alfa_n+2*bata_n;
                       integ_lambda=integ_lambda+lambda*delta_t;
                       F=exp(-integ_lambda);
                       n_i=n_i+1;
                   end
           
                  %at the jump time , system jump to "1" or "3"; with
                 %probability ...rate 
            
            
            
                   U2=rand;
                   if U2<=2*alfa_n/lambda;
                       S_k(n_i)=3;   %jump to state "3"
                   else
                       S_k(n_i)=1;   %jump to state "1"
                   end 
            
            
            
             %-------------------------  3
              case 3
                   %calculate time that the channel stay at "3"   lambda=1*alfa_n+3*bata_n; 
                   U1=rand;
                   F=1;
                   integ_lambda=0;        %the integrel of lambda(s)
                   while F>U1  &&   n_i<T_step
                       U3=normrnd(0,sqrt(delta_t),1);
                       V(n_i+1)=V(n_i)+(-Cond_L*(V(n_i)-V_L))*delta_t+delta_1*U3;
                        if V(n_i)>V_L
                             V(n_i)=(V_L-V(n_i))+V_L;            
                        end
      %-------------------------  
                           if V(n_i)<V_e
                                  V(n_i)=(V_e-V(n_i))+V_e;              
                           end
                       S_k(n_i+1)=S_k(n_i);
                       alfa_n=0.01*(V(n_i)+55)/(1-exp(-(V(n_i)+55)/10));
                       bata_n=0.125*exp(-(V(n_i)+65)/80); 
                       lambda=alfa_n+3*bata_n;
                       integ_lambda=integ_lambda+lambda*delta_t;
                       F=exp(-integ_lambda);
                       n_i=n_i+1;
                   end
           
                  %at the jump time , system jump to "2" or "4"; with
                 %probability ...rate 
            
                   U2=rand;
                   if U2<=alfa_n/lambda;
                       S_k(n_i)=4;   %jump to state "4"
                   else
                       S_k(n_i)=2;   %jump to state "2"
                   end 
         
            
            
             %-------------------------  4
              case 4
                    %calculate time that the channel stay at "4"   lambda=4*bata_n; 
                   U1=rand;
                   F=1;
                   integ_lambda=0;        %the integrel of lambda(s)
                   while F>U1    &&   n_i<T_step
                       U3=normrnd(0,sqrt(delta_t),1);
                       V(n_i+1)=V(n_i)+(-Cond_k*(V(n_i)-V_e))*delta_t+delta_2*U3;
                        if V(n_i)>V_L
                             V(n_i)=(V_L-V(n_i))+V_L;            
                        end
      %-------------------------  
                           if V(n_i)<V_e
                                  V(n_i)=(V_e-V(n_i))+V_e;              
                           end
                       S_k(n_i+1)=S_k(n_i); 
                       bata_n=0.125*exp(-(V(n_i)+65)/80);
                       lambda=4*bata_n;
                       integ_lambda=integ_lambda+lambda*delta_t;
                       F=exp(-integ_lambda);
                       n_i=n_i+1;
                   end
             %   %time that the system stay in state "4"
            %at jump time , system jump to  "3"; definitely
                  S_k(n_i)=3;
         
            
   end

end
subplot(2,1,1);
plot(S_k);   
xlabel('Time(ms)');
ylabel('States');
set(gca,'ytick',[0 1 2 3 4]) %x坐标轴上刻度的数据点位置
set(gca,'yticklabel',{'N0','N1','N2','N3','N4'});
% set(gca,'xlim',[0,500000]);%x坐标轴范围
subplot(2,1,2);                                                                                 
plot(V);   
xlabel('Time(ms)');
ylabel('Voltage(mV)');
% set(gca,'xlim',[0,500000]);%x坐标轴范围
%最后更改坐标轴刻度*0.001,这里的原横坐标是T_step,而不是T_step*delta_t