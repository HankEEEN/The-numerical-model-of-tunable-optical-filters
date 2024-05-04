clear all;
fs = 0.90; fe = 1.1; 
Np = 10000;
Normalized_frequency = fs:(fe-fs)/Np:fe;
%Sets the refractive indices
n1 = 3.717; n2 = 3.5;     %kL = 3
% n1 = 3.7; n2 = 3.5907;    %kL = 1.5
%Calculate the interface constants S
S12=(1/(2*n1))*[n1+n2 n1-n2;n1-n2 n1+n2];
S21=(1/(2*n2))*[n2+n1 n2-n1;n2-n1 n2+n1];
%Set the gain coefficient g (g_frequency = 1/2 g_power)
gL = 0.6;
g = gL / 200;


for i = 1:length(Normalized_frequency)                                              
fn=Normalized_frequency(i);    
%Calculate the optical thickness
theta1 = (pi/2)*fn+1i*g;       %For slot           
theta2_1=(pi/2)*fn+1i*g;       %For lamda         
theta2_2=(pi)*fn+1i*2*g;       %For a quarter-wave phase-shift
%Calculate the medium constant P
Pn1 = [exp(1i*theta1) 0;0 exp(-1i*theta1)];
Pn2_1=[exp(1i*theta2_1) 0;0 exp(-1i*theta2_1)];
Pn2_2=[exp(1i*theta2_2) 0;0 exp(-1i*theta2_2)];
%Construct the M matrix of the whole strcture
M_left = ((Pn1*S12*Pn2_1*S21)^24);
M_mid = Pn1*S12*Pn2_2*S21;
M_right = ((Pn1*S12*Pn2_1*S21)^25);
M = M_left * M_mid * M_right;
%Calculate the reflection coefficient
Ramp = M(2,1) / M(1,1);
%Calculate the power reflection coefficient
Rpow = (abs(Ramp)).^2;
Reflectance(i)=Rpow;
end

%Set the gain coefficient g (g_frequency = 1/2 g_power)
gL = 0;
g = gL / 200;

for i = 1:length(Normalized_frequency)                                              
fn=Normalized_frequency(i);    
%Calculate the optical thickness
theta1 = (pi/2)*fn+1i*g;       %For slot           
theta2_1=(pi/2)*fn+1i*g;       %For lamda         
theta2_2=(pi)*fn+1i*g;         %For a quarter-wave phase-shift
%Calculate the medium constant P
Pn1 = [exp(1i*theta1) 0;0 exp(-1i*theta1)];
Pn2_1=[exp(1i*theta2_1) 0;0 exp(-1i*theta2_1)];
Pn2_2=[exp(1i*theta2_2) 0;0 exp(-1i*theta2_2)];
%Construct the M matrix of the whole strcture
M_left = (Pn1*S12*Pn2_1*S21)^24;
M_mid = Pn1*S12*Pn2_2*S21;
M_right = (Pn1*S12*Pn2_1*S21)^25;
M = M_left * M_mid * M_right;
%Calculate the reflection coefficient
Ramp = M(2,1) / M(1,1);
%Calculate the power reflection coefficient
Rpow = (abs(Ramp)).^2;
Reflectance_1(i)=Rpow;
end

figure;
% plot(Normalized_frequency,Reflectance,'r','LineWidth',1);
% hold on;
plot(Normalized_frequency,Reflectance_1,'k','LineWidth',1);
xlim([fs fe]);
% set(gca, 'YScale', 'log');
% ylim([10^(-4) 10^3]);
xlabel('Normalized_frequency (f/fB)');                                                   
ylabel('Power reflectivity');







