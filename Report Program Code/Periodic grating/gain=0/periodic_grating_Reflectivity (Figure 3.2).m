clear all;

fs = 0; fe = 2; 
Np = 1e6;
Normalized_frequency = fs:(fe-fs)/Np:fe;
%Sets the refractive indices
n1 = 3.717; n2 = 3.5;     %kL = 3
% n1 = 3.7; n2 = 3.5907;    %kL = 1.5
%Number of the grating layers
N = 50;
%Calculate the interface constants S
S12=(1/(2*n1))*[n1+n2 n1-n2;n1-n2 n1+n2];
S21=(1/(2*n2))*[n2+n1 n2-n1;n2-n1 n2+n1];
%Set the gain coefficient g (g_frequency = 1/2 g_power)
gL = 1.2;
g = gL / (4 * N);

for i = 1:length(Normalized_frequency)                                              
fn=Normalized_frequency(i);    
%Calculate the optical thickness
theta1 = (pi/2)*fn+1i*g;       %For slot     lamda/4
theta2 = (pi/2)*fn+1i*g;       %For lamda    lamda/4
%Calculate the medium constant P
Pn1 = [exp(1i*theta1) 0;0 exp(-1i*theta1)];
Pn2 = [exp(1i*theta2) 0;0 exp(-1i*theta2)];
%Construct the M matrix of the whole strcture
M = (Pn1*S12*Pn2*S21)^N;
%Calculate the reflection coefficient
Ramp = M(2,1) / M(1,1);
%Calculate the power reflection coefficient
Rpow = (abs(Ramp)).^2;
Reflectance(i)=Rpow;
end


%Set the gain coefficient g (g_frequency = 1/2 g_power)
gL = 0;
g = gL / (4 * N);

for i = 1:length(Normalized_frequency)                                              
fn=Normalized_frequency(i);    
%Calculate the optical thickness
theta1 = (pi/2)*fn+1i*g;       %For slot     lamda/4
theta2 = (pi/2)*fn+1i*g;       %For lamda    lamda/4
%Calculate the medium constant P
Pn1 = [exp(1i*theta1) 0;0 exp(-1i*theta1)];
Pn2 = [exp(1i*theta2) 0;0 exp(-1i*theta2)];
%Construct the M matrix of the whole strcture
M = (Pn1*S12*Pn2*S21)^N;

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
xlim([0.9 1.1]);
% set(gca, 'YScale', 'log');
% ylim([10^(-6) 10^4]);
xlabel('Normalized_frequency (f/fB)');                                                   
ylabel('Power reflectivity');






