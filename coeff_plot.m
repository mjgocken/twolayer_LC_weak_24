% plot A/B v thetaB

% close all, clear, clc

% alpha1_dim = -0.006*3;      % alpha_1 (Pa.s)
% alpha2_dim = -0.0812*3;     % alpha_2 (Pa.s)
% alpha3_dim = -0.0036*3;     % alpha_3 (Pa.s)
% alpha4_dim = 0.0652*3;      % alpha_4 (Pa.s)
% alpha5_dim = -0.064*3;      % alpha_5 (Pa.s)
% alpha6_dim = -0.0208*3;     % alpha_6 (Pa.s)

alpha1_dim = -0.006/2;      % alpha_1 (Pa.s)
alpha2_dim = -0.0812/2;     % alpha_2 (Pa.s)
alpha3_dim = -0.0036/2;     % alpha_3 (Pa.s)
alpha4_dim = 0.0652/2;      % alpha_4 (Pa.s)
alpha5_dim = -0.064/2;      % alpha_5 (Pa.s)
alpha6_dim = -0.0208/2;     % alpha_6 (Pa.s)

% nondimensionalized viscosities 
alpha1 = alpha1_dim/alpha4_dim;
alpha2 = alpha2_dim/alpha4_dim;
alpha3 = alpha3_dim/alpha4_dim;
alpha4 = alpha4_dim/alpha4_dim;
alpha5 = alpha5_dim/alpha4_dim;
alpha6 = alpha6_dim/alpha4_dim;

A =@(thetaB) 2*cos(2*thetaB)*(alpha1+alpha5+alpha6+4)*(alpha2+alpha3-alpha5+alpha6) + cos(4*thetaB)*(alpha1*(alpha2-alpha3)...
    + (alpha2+alpha3)*(alpha5-alpha6)) + alpha1*alpha2 - alpha1*alpha3...
    - 2*alpha1*alpha5 - 2*alpha1*alpha6 - 8*alpha1 + alpha2*alpha5 + 3*alpha2*alpha6 + 8*alpha2 - 3*alpha3*alpha5...
    - alpha3*alpha6 - 8*alpha3 - 2*alpha5^2 - 4*alpha5*alpha6 - 16*alpha5 - 2*alpha6^2 - 16*alpha6 - 32;           % coeff in u2 eqn
B =@(thetaB) -alpha1*cos(4*thetaB) + alpha1 - 2*cos(2*thetaB)*(alpha2+alpha3-alpha5+alpha6) + 2*(-alpha2+alpha3+alpha5+alpha6) + 8;            % coeff in u2 eqn

thetaB = 0:pi/200:pi/2;
coeff = A(thetaB)./B(thetaB);

%plot coefficient function
figure();
plot(thetaB,coeff,'ok','linewidth',2.5)

%function is not symmetric
sum(coeff(1:50)-flip(coeff(52:end)))

%plot the difference
figure();
plot(thetaB(1:50),coeff(1:50)-flip(coeff(52:end)))

%plot the first half and the mirror image of the second half of the
%function
figure();
plot(thetaB(1:50),coeff(1:50))
hold on
plot(thetaB(1:50),flip(coeff(52:end)))