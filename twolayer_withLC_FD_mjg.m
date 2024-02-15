%% Extensional flow simulation
% updated: 6/19/23
%
% Implementing equations for a two layer model of the tear film. The aqueous
% layer is Newtonian, and the lipid layer is liquid crystal (weak
% elasticity).
% 
% Unknowns:  h1,p1,h2,u2,g,m=c*h1 
%
% The PDE problem to solve is 
%
% $$ h1_t - v0 (x/s) h1_x + (1/s)(u1ave h1)_x = Jo - Je $$ 
%
% $$ h2_t - v0 (x/s) h2_x + (1/s)(u2 h2)_x = 0 $$ 
%
% $$ p1 = -Upsilon C2 (1/s^2) h_xx - C1 (1/s^2) h_1xx $$
%
% $$ -A/B d Upsilon (1/s^2) (u2x h2)_x - u2/h1 = (1/2s) h1 p1_x - C2 d Upsilon (1/s^3) h2 hxxx + (1/s) M g_x $$
%
% $$ m_t - v0 (x/s) m_x = (1/Pe1) (1/s^2) (m_x - m h1_x / h1)_x - (1/s) (m u1ave)_x $$
%
% $$ g_t - v0 (x/s) g_x + (1/s)(g u2)_x = (1/Pe2)(1/s^2)g_xx $$
%
% $$ 0 < y < 1 $$
% 
% bcs: $$  $$ 
% 
% ics: 
%
% We have the following definitions:
%
% h = h1 + d h2
%
% Jo = Op (m/h1 - 1)
% 
% Je = E / (1 + R h2)

%% setup
%
% clean up and set parameters

close all, clear, clc

% paper parameters
% k=1;j=1;
% for thetaB = 0:pi/100:pi/2
% for R=10.^(-2:.1:4)
% for M=10.^(-3:.1:3)
% for H2=(5:5:210)./1e9 
% for sigma = 0.03:0.1:2.9
% physical constants
v0 = 0;           % speed of moving end 

%% CASE 1: leslie viscosities for 5CB from Stewart
% alpha1_dim = -0.006;      % alpha_1 (Pa.s)
% alpha2_dim = -0.0812;     % alpha_2 (Pa.s)
% alpha3_dim = -0.0036;     % alpha_3 (Pa.s)
% alpha4_dim = 0.0652;      % alpha_4 (Pa.s)
% alpha5_dim = -0.064;      % alpha_5 (Pa.s)
% alpha6_dim = -0.0208;     % alpha_6 (Pa.s)

%% CASE 2: leslie viscosities for 5CB from Stewart but with a4 = 2*mu2 from Newtonian case
% alpha1_dim = -0.006;      % alpha_1 (Pa.s)
% alpha2_dim = -0.0812;     % alpha_2 (Pa.s)
% alpha3_dim = -0.0036;     % alpha_3 (Pa.s)
% alpha4_dim = 2*0.1;      % alpha_4 (Pa.s)
% alpha5_dim = -0.064;      % alpha_5 (Pa.s)
% alpha6_dim = -0.0208;     % alpha_6 (Pa.s)

%% CASE 3: leslie viscosities for 3 times 5CB from Stewart
scl=3;
alpha1_dim = -0.006*scl;      % alpha_1 (Pa.s)
alpha2_dim = -0.0812*scl;     % alpha_2 (Pa.s)
alpha3_dim = -0.0036*scl;     % alpha_3 (Pa.s)
alpha4_dim = 0.0652*scl;      % alpha_4 (Pa.s)
alpha5_dim = -0.064*scl;      % alpha_5 (Pa.s)
alpha6_dim = -0.0208*scl;     % alpha_6 (Pa.s)

%% CASE 4: newtonian case
% alpha1_dim = 0;     % alpha_1 (Pa.s)
% alpha2_dim = 0;     % alpha_2 (Pa.s)
% alpha3_dim = 0;     % alpha_3 (Pa.s)
% alpha4_dim = 0.0652*3;     % alpha_4 (Pa.s)
% alpha5_dim = 0;     % alpha_5 (Pa.s)
% alpha6_dim = 0;     % alpha_6 (Pa.s)

% nondimensionalized viscosities 
alpha1 = alpha1_dim/alpha4_dim;
alpha2 = alpha2_dim/alpha4_dim;
alpha3 = alpha3_dim/alpha4_dim;
alpha4 = alpha4_dim/alpha4_dim;
alpha5 = alpha5_dim/alpha4_dim;
alpha6 = alpha6_dim/alpha4_dim;
thetaB = pi/2;       % angle of the director
% thetaB=0;

A = 2*cos(2*thetaB)*(alpha1+alpha5+alpha6+4)*(alpha2+alpha3-alpha5+alpha6) + cos(4*thetaB)*(alpha1*(alpha2-alpha3)...
    + (alpha2+alpha3)*(alpha5-alpha6)) + alpha1*alpha2 - alpha1*alpha3...
    - 2*alpha1*alpha5 - 2*alpha1*alpha6 - 8*alpha1 + alpha2*alpha5 + 3*alpha2*alpha6 + 8*alpha2 - 3*alpha3*alpha5...
    - alpha3*alpha6 - 8*alpha3 - 2*alpha5^2 - 4*alpha5*alpha6 - 16*alpha5 - 2*alpha6^2 - 16*alpha6 - 32;           % coeff in u2 eqn
B = -alpha1*cos(4*thetaB) + alpha1 - 2*cos(2*thetaB)*(alpha2+alpha3-alpha5+alpha6) + 2*(-alpha2+alpha3+alpha5+alpha6) + 8;            % coeff in u2 eqn

H1 = 3.5e-6;      % AL thickness (m)
H2 = 5e-8;        % LL thickness (m)
E0 = 5.093e-7;    % thinning rate (m/s)
V = 5.093e-7;     % thinning rate (m/s)
U = 4.624e-5;     % horizontal velocity (m/s)
C = 300;          % salt concentration (mOsM)
Gamma0 = 4e-7;    % surfactant concentration (mol/m^2)
mu1 = 1.3e-3;     % dynamic viscosity of water (Pa*s)
mu2 = alpha4_dim/2;       % dynamic viscosity of lipid (Pa*s) calculated from leslie viscosity
gam1 = 0.027;     % surface tension of water (N/m)
gam2 = 0.018;     % surface tension of lipid (N/m)
D1 = 1.6e-9;      % salt diffusion coefficient (m^2/s)
D2 = 3e-8;        % polar lipid diffusion constant (m^2/s)
Igas = 8.3145;    % ideal gas constant (J/(K*mol))
T0 = 273.15+35;   % eye surface temperature (K)
Pc = 2.3e-7;      % corneal osmosis coeff (kg/(mOsM*m^2*s))
rho1 = 1000;      % density of liquid water (kg/m^3)
rho2 = 900;       % density of lipid (kg/m^3)
km = 0.018157;    % mass transfer coeff (m/s)
Dk = 5.136e-11;   % lipid permeability (m^2/s)

% nondimensional parameters
L = H1*( (gam1+gam2)/(mu1*V) )^(1/4); % length scale (m)
epsilon = H1/L;                 % aspect ratio of the AL (0.011)
d = H2/H1;                      % lipid to aqueous thickness ratio (0.0143);
t = L/U;                        % time scale (s)
Upsilon = epsilon^2*mu2/mu1;    % viscosity ratio (0.0093)
C1 = epsilon^3*(gam1 + gam2)/(mu1*U); % aqueous capillary number (1)
C2 = epsilon*gam2/(mu2*U);      % lipid capillary number (42.871)
Pe1 = (U*L)/D1;                 % Peclet number for salt diff. (9.185)
Pe2 = (U*L)/D2;                 % Peclet number for surfactant diff. (0.4898)
M = epsilon*Igas*T0*Gamma0/(mu1*U);   % Marangoni number (187.7) 
Op = C*Pc/(rho1*epsilon*U);     % Osmosis parameter (0.135)
E = E0/(epsilon*U);             % Evaporation parameter (1)
R = km*H2/(Dk);                 % evaporative resistance of LL (17.68)

nu = .5;          % for robin boundary condition
P = 2;            % controls period of initial condition: h0 = aa+bb*cos(P*pi*y);

% Set grid size and solver tolerances; includes event detection
N = 512;         %  grid size: number of intervals
a = -1/(L*1000); % right endpoint; 1mm
% a = -3/(L*1000);   % 3 mm for Fig 12
b = -a;          % left endpoint
% a = 0/(L*1e3);
% b = 2/(L*1e3);
y = linspace(a,b,N+1)';

% Define constant mass matrix
% This is how Matlab knows that this is a system of DAEs.
% Ones on diagonal where equations have time derivatives.
% Zero rows where BCs are applied, and for all p1 and u2 equations.
% BCs at both ends, so rows for 1st and last of each variable are zero
% 6 PDE: h1, h2, p1, u2, m, g
Ms=speye(6*(N+1));
Ms(1,1) = 0;                    % for h1 bc at right end
Ms(N+1,N+1) = 0;                % for h1 bc at left end (origin)
Ms(N+2,N+2) = 0;                % for h2 bc at right end
Ms(2*(N+1),2*(N+1)) = 0;        % for h2 bc at left end
for i=2*N+3:4*(N+1)
    Ms(i,i)=0;    % all points for p1 and u2
end
Ms(4*N+5,4*N+5) = 0;            % for m bc at right end
Ms(5*(N+1),5*(N+1)) = 0;        % for m bc at left end
Ms(5*N+6,5*N+6) = 0;            % for g bc at right end
Ms(6*(N+1),6*(N+1)) = 0;        % for g bc at left end

%spy(Ms)
 
MyTols = odeset('RelTol',5e-4,'AbsTol',1e-4,'Mass',Ms,'MStateDependence','none',...
       'MassSingular','yes','MaxOrder',2,'BDF','on','Events',@(t,w)zerothickness(w,N,epsilon,L)); 
   
% Pick one for BC for h:
    setbc = 0;    % at left: h2x = 0, at right: h2x = 0
%     setbc = 1;    % at left: h2x != 0, at right: h2x != 0
%     setbc = 2;    % at left: h2x = 0, at right: h2 = 1
%     setbc = 5;    % at left: h2x = 0, at right: nu*h2x+(1-nu)*(h2-1) = 0
%     setbc = 6;    % at left: -nu*h2x+(1-nu)*(h2-1) = 0, at right: h2x = 0

% Time steps for output--pick one:
% Regular steps :
tend = 60*U/L; dtp = tend/6; tspan = 0:dtp:tend; % could make finer
% tend = 0.060*U/L; dtp = tend/6; tspan = 0:dtp:tend; % short ts
% tspan = [0:.1:60].*U/L;
% tspan=[1,10,60].*U/L;
% tspan = (2:10:52)*(U/L); % fig 7/8
% tspan = [0 0.01 0.02 0.05 0.1 0.25 0.5 1 2 3];  
% tspan = [0. 0.00001 0.0001 0.001 0.01 0.1 0.25 0.5 1 4 8];


%% spatial discretization. 
% variables are for fixed domain
% $$ h(y(t),t),\ u(y(t),y) $$ are on fixed domain.
% $$ y = x/s(t),\ s(t) = 1+v0*t $$ ($$ y $$ is $$ \xi $$ in paper)

% delta x
delx = (b-a)/N;

% first derivative matrix
% All rows treated as interior first, replace 1st and last below
D1 = (1/(2*delx))*full(spdiags([-ones(N+1,1),ones(N+1,1)],[-1,1],N+1,N+1));  % (N+1)x(N+1) first derivative matrix using second order center diff; 
D1(1,:) = (1/(2*delx))*[-3 4 -1 zeros(1,N-2)];    % one sided 1st derivative at left end
D1(end,:) = (1/(2*delx))*[zeros(1,N-2) 1 -4 3];   % one sided 1st derivative at right end

% second derivative matrix
a2 = [2 -5 4 -1 zeros(1, N-3)]; %one sided 2nd derivative at left end
b2 = full(spdiags([ones(N+1,1), -2*ones(N+1, 1), ones(N+1,1)],[0, 1, 2],N-1,N+1)); % 2nd derivative centered
c2 = [zeros(1,N-3) -1 4 -5 2]; %one sided 2nd derivative on right end
D2 = (1/delx^2)*[a2;b2;c2]; %(N+1)x(N+1)

% third derivative matrix
% a3 is a 2x(N+1) matrix: forward fd for end, then fd using one point from left, 3 from right for second from the end
a3 = [-5 18 -24 14 -3 zeros(1,N-4); -3 10 -12 6 -1 zeros(1,N-4)]; % interior points, centered
b3 = full(spdiags([-ones(N+1,1),2*ones(N+1,1),-2*ones(N+1,1),ones(N+1,1)],[0,1,3,4],N-3,N+1)); % (N-3)x(N+1) second order center diff
% c3 is a 2x(N+1) matrix: uses 3 points from left, one from right for second from end, and backward fd for right end
c3 = [zeros(1,N-4) 1 -6 12 -10 3; zeros(1,N-4) 3 -14 24 -18 5]; 
% assemble D3 matrix
D3 = (1/(2*delx^3))*[a3;b3;c3]; % third derivative matrix (N+1)x(N+1)

% fourth derivative
a4 = [3 -14 26 -24 11 -2 zeros(1,N-5); 2 -9 16 -14 6 -1 zeros(1,N-5)];
b4 = full(spdiags([ones(N+1,1),-4*ones(N+1,1),6*ones(N+1,1),-4*ones(N+1,1),ones(N+1,1)],[0,1,2,3,4],N-3,N+1));
c4 = [zeros(1,N-5) -1 6 -14 16 -9 2; zeros(1,N-5) -2 11 -24 26 -14 3];
D4 = (1/delx^4)*[a4;b4;c4];


%% ICs 

% for h1
AA = 1; BB = 0; P1 = 2;
h10 = AA+BB*cos(P1*pi*y);
h10prime = -P1*pi*BB*sin(P1*pi*y);
h10dprime = -(P1*pi)^2*BB*cos(P1*pi*y);
h10tprime = (P1*pi)^3*BB*sin(P1*pi*y);
h10fprime = (P1*pi)^4*BB*cos(P1*pi*y);

% for h2
if setbc == 1 % needs quadratic term to be consistent
    aa = .9; bb = 0.1; dd = 0.1;
    h20 = aa+bb*cos(P*pi*y)+dd*(y.*(y-1));
    h20prime = -P*pi*bb*sin(P*pi*y)+dd(2*y-1);
    h20dprime = -(P*pi)^2*bb*cos(P*pi*y)+2*dd;
    h20tprime = (P*pi)^3*bb*sin(P*pi*y);
    h20fprime = (P*pi)^4*bb*cos(P*pi*y);
else

    %***exact IC from Stapf paper figures!!***
    dd=0;
    mu = (b+a)/2;
    sigma = .2364; % This is used for stapf figures!
    %sigma = (max(y)-min(y))/2*(0.2/3); % correct format for stapf
    gaussian = 0.9*exp(-(y-mu).^2/(2*sigma^2));
    h20 = 1 - gaussian;
    h20prime = (y-mu)/sigma^2.*gaussian;
    h20dprime = ( 1/sigma^2 - (y-mu).^2/sigma^4 ).*gaussian;
    h20tprime = ( -(3*(y-mu))/sigma^4 + (y-mu).^3/sigma^6 ).*gaussian;

% %   ***flat IC***
%     h20 = ones(N+1,1);
%     h20prime = zeros(N+1,1);
%     h20dprime = zeros(N+1,1);
%     h20tprime = zeros(N+1,1);


end

h0 = h10 + d*h20; % define total film thickness
h0dprime = h10dprime + d*h20dprime; 
h0tprime = h10tprime + d*h20tprime; 

g0 = ones(N+1,1); % initial surfactant concentration
g0prime = zeros(N+1,1);

% %%*** tanh IC for surfactant ***
% dee = 0.9;
% w = 1;
% mu = 2.2;
% 
% [gamIC, min_gamIC, max_gamIC] = tanh_dist(y, w, mu, dee);
% g0 = 1 + 2 - 2*gamIC;
% g0prime = (-2*dee/max_gamIC)*(sech(6*(y-w)/mu).^2 - sech(6*(-y-w)/mu).^2)*(3/mu);
% plot(y,g0,y,g0prime)

% pressure
p10 = -Upsilon*C2*h0dprime - C1*h10dprime; % initial AL pressure
p10prime = -Upsilon*C2*h0tprime - C1*h10tprime;

% Solve for the initial u2
lhs = - (A/B)*d*Upsilon*(diag(h20)*D2 + diag(h20prime)*D1)-diag(1./h10);
lhs(1,:) = [1 zeros(1,N)];
lhs(end,:) = [zeros(1,N) 1];
rhs = p10prime.*h10./2 - C2*d*Upsilon*h20.*h0tprime + M*g0prime;
rhs(1,1) = 0;
rhs(end,1) = v0;
u20 = lhs\rhs;

c0 = 1;      % initial salt concentration
m0 = h10*c0; % initial mass m=h10*c0

w0 = [h10; h20; p10; u20; m0; g0];      % ic for solver

% plotic(y,h10,h20,u20,p10,m0,g0)

%% call solver

[T, W, te, ye, ie] = ode15s(@(t,w) efnewt(t,w,y,v0,Upsilon,C1,C2,Pe1,Pe2,A,B,d,M,dd,Op,E,R,D1,D2,D3,setbc,nu,thetaB),tspan,w0,MyTols);
    % te: column vector of the times at which events occurred.
    % ye: the solution value at each of the event times in te.
    % ie: indices into the vector returned by the event function. The values indicate which event the solver detected.

% readable variables 
hh1 = W(:,1:N+1);  hh2 = W(:,N+2:2*(N+1)); pp1 = W(:,2*N+3:3*(N+1)); 
uu2 = W(:,3*N+4:4*(N+1)); mm = W(:,4*N+5:5*(N+1)); gg = W(:,5*N+6:6*(N+1)); 
y2 = ( 1+v0*T' ).*y; % time on a moving domain

% calculate c
cc = mm./hh1;

%dimensionalize
% t_dim = (L/U)*T; % (s)
% h1_dim = 1e6*epsilon*L*hh1; % (micrometers)
% h2_dim = 1e9*epsilon*d*L*hh2; % (nanometers)
% p1_dim = (mu1*U)/(epsilon^2*L)*pp1; % (Pa)
% u2_dim = 6e7*U*uu2; % (micrometers/min)
% g_dim = Gamma0*gg; % (mol/m^2)
% c_dim = C*cc; % (mOsM)
% y_dim=L*1000*(1+v0*T)*y'; % (mm)
% 
% save LC_3times5CB_thetaDepEvap_thetaBpi18_dim_finer t_dim h1_dim h2_dim p1_dim u2_dim g_dim c_dim y_dim
% save LC_53timesCB_thetaDepEvap_thetaB0 T hh1 hh2 pp1 uu2 gg cc
% save newtonian_thetaDepEvap_thetaBpi2_dim t_dim h1_dim h2_dim p1_dim u2_dim g_dim c_dim y_dim


%% check volume at all timepoints
for i = 1:length(T)
    y3 =( 1+v0*T(i) )*y;
    fprintf('At time = %.1f, the volume of the AL  is %f, and of m is %f\n',T(i),trapz(y3,hh1(i,:)),trapz(y3,mm(i,:)))
end

%% Plot results
% % plots all times in tspan.  h first, then u, then p, and perhaps more.
% % 
% % %plot initial conditions
% plotic(y,h0,u0,p0,D1,D2,h0prime,p0prime,p0dprime)
%    
% dimensional plots
moving_dim(y,T,hh1,hh2,uu2,pp1,mm,gg,cc,v0,L,U,epsilon,d,mu1,Gamma0,C)
% 
% % yyplot dimensional
% consolidated(y,T,hh1,hh2,uu2,pp1,mm,gg,cc,v0,L,U,epsilon,d,mu1,Gamma0,C,t,H1,H2,V)


%% plot all components of equation for u2
% figure();
% plot(y,Upsilon*d*D1*((D1*uu2')'.*hh2)','linewidth',2)
% figure();
% plot(y,(hh2.*(D2*uu2')')','linewidth',2)
% figure();
% plot(y,((D1*hh2').*(D1*uu2'))','linewidth',2)
% figure();
% plot(y,(uu2.*hh1)','linewidth',2)
% figure();
% plot(y,(hh1.*(D1*pp1')')'./2,'linewidth',2)
% figure();
% plot(y,(hh2.*(D3*(hh1+d*hh2)')')','linewidth',2)
% figure();
% plot(y,D1*gg','linewidth',2)

%% plot thetaB vs h1min, m2min, max osmolarity

% if T(end) == 60*U/L
%     h1min(k)=min(hh1(end,:));
%     h2min(k)=min(hh2(end,:));
%     cmax(k)=max(cc(end,:));
%     thetaBval(k)=thetaB;
%     k=k+1;
% else
%     TBUT(j)=T(end);
%     TBUTthetaB(j)=thetaB;
%     j=j+1;
%     h1min(k)=min(hh1(end,:));
%     h2min(k)=min(hh2(end,:));
%     cmax(k)=max(cc(end,:));
%     thetaBval(k)=thetaB;
%     k=k+1;
% end
% % 
% % h1min(k)=min(hh1(end,:));
% % h2min(k)=min(hh2(end,:));
% % cmax(k)=max(cc(end,:));
% % thetaBval(k)=thetaB;
% % k=k+1;
% end
% 
% % AL min at final time
% figure();
% plot(thetaBval,h1min*H1*1e6,'linewidth',2)
% xlabel('$\theta_B$','Interpreter','LaTex')
% ylabel('Min AL thickness $(\mu m)$','Interpreter','LaTex')
% set(gca,'FontSize',20)
% set(gca,'XTick',[0,pi/6,pi/3,pi/2]) 
% set(gca,'XTickLabel',{'0','pi/6','pi/3','pi/2'})
% 
% % AL min at final time
% figure();
% yyaxis right
% plot(thetaBval,h1min*H1*1e6,'linewidth',2)
% xlabel('\theta_B')
% ylabel('Min AL thickness (\mu m)')
% ylim([0.45,2.1])
% yyaxis left
% plot(TBUTthetaB,TBUT*t,'linewidth',2)
% ylabel('TBUT (s)')
% ylim([0,60])
% set(gca,'FontSize',20)
% set(gca,'XTick',[0,pi/6,pi/3,pi/2]) 
% set(gca,'XTickLabel',{'0','\pi/6','\pi/3','\pi/2'})
% 
% % LL min at final time
% figure();
% plot(thetaBval,h2min*H2*1e9,'linewidth',2)
% xlabel('$\theta_B$','Interpreter','LaTex')
% ylabel('Min LL thickness $(nm)$','Interpreter','LaTex')
% set(gca,'FontSize',20)
% set(gca,'XTick',[0,pi/6,pi/3,pi/2]) 
% set(gca,'XTickLabel',{'0','pi/6','pi/3','pi/2'})
% 
% figure();
% yyaxis right
% plot(thetaBval,h2min*H2*1e9,'linewidth',2)
% xlabel('\theta_B')
% ylabel('Min LL thickness (nm)')
% ylim([6,22])
% yyaxis left
% plot(TBUTthetaB,TBUT*t,'linewidth',2)
% ylabel('TBUT (s)')
% ylim([0,60])
% set(gca,'FontSize',20)
% set(gca,'XTick',[0,pi/6,pi/3,pi/2]) 
% set(gca,'XTickLabel',{'0','\pi/6','\pi/3','\pi/2'})
% 
% 
% % max osmolarity at final time
% figure();
% plot(thetaBval,cmax*C,'linewidth',2)
% xlabel('$\theta_B$','Interpreter','LaTex')
% ylabel('Max osmolarity $(mOsM)$','Interpreter','LaTex')
% set(gca,'FontSize',20)
% set(gca,'XTick',[0,pi/6,pi/3,pi/2]) 
% set(gca,'XTickLabel',{'0','pi/6','pi/3','pi/2'})
% 
% figure();
% yyaxis right
% plot(thetaBval,cmax*C,'linewidth',2)
% xlabel('\theta_B')
% ylabel('Max Osmolarity (mOsM)')
% ylim([500,1700])
% yyaxis left
% plot(TBUTthetaB,TBUT*t,'linewidth',2)
% ylabel('TBUT (s)')
% ylim([0,60])
% set(gca,'FontSize',20)
% set(gca,'XTick',[0,pi/6,pi/3,pi/2]) 
% set(gca,'XTickLabel',{'0','\pi/6','\pi/3','\pi/2'})


%% Fig 2 from Stapf
% plot(hh2(end,end)*1e9*epsilon*d*L, (E./(1+R*hh2(end,end)))*epsilon*U*1e6*60,'ko','linewidth',2)
% hold on
% xlabel('$h_2\;(nm)$', 'Interpreter', 'LaTeX')
% ylabel('$J_e\;(micrometer/min)$', 'Interpreter', 'LaTeX')
% end

% semilogx(R,1e6*epsilon*L*hh1(end,end),'ko','linewidth',2)
% hold on
% xlabel('$\mathcal{R}$', 'Interpreter', 'LaTeX')
% ylabel('$h_1 \;(\mu m)$', 'Interpreter', 'LaTeX')
% end

%% Fig 5 from Stapf 

% if T(end) == 60*U/L
%     h1min(k)=min(hh1(end,:));
%     h1minR(k)=R;
%     k=k+1;
% else
%     TBUT(j)=T(end);
%     TBUTR(j)=R;
%     j=j+1;
% end
% end
% 
% TBUT_dim = (L/U)*TBUT; % (s)
% h1min_dim = 1e6*epsilon*L*h1min; % (micrometers)
% 
% save fig5_LC_3times5CB_thetaDepEvap_thetaBpi18_dim h1min_dim h1minR TBUT_dim TBUTR
% 
% figure();
% yyaxis left
% semilogx(TBUTR,TBUT*t,'linewidth',2)
% ylabel('TBUT (s)','Interpreter','LaTex')
% yyaxis right
% semilogx(h1minR,h1min*H1*1e6,'linewidth',2)
% ylabel('$h_1^{min}\;(\mu m)$','Interpreter','LaTex')
% xlabel('$\mathcal{R}$','Interpreter','LaTex')

%% Fig 10 from Stapf 

% if T(end) == 60*U/L
%     h1min(k)=min(hh1(end,:));
%     h1minM(k)=M;
%     k=k+1;
% else
%     TBUT(j)=T(end);
%     TBUTM(j)=M;
%     j=j+1;
% end
% end
% 
% figure();
% yyaxis left
% semilogx(h1minM,h1min*H1*1e6,'linewidth',2)
% ylabel('$h_1^{min}\;(\mu m)$','Interpreter','LaTex')
% yyaxis right
% semilogx(TBUTM,TBUT*t,'linewidth',2)
% ylabel('TBUT (s)','Interpreter','LaTex')
% xlabel('$\mathcal{M}$','Interpreter','LaTex')

%% Fig 12 from Stapf 

% if T(end) == 60*U/L
%     h1min(k)=min(hh1(end,:));
%     h1minsigma(k)=sigma;
%     k=k+1;
% else
%     TBUT(j)=T(end);
%     TBUTsigma(j)=sigma;
%     j=j+1;
% end
% end
% 
% figure();
% yyaxis left
% plot(h1minsigma,h1min*H1*1e6,'linewidth',2)
% ylabel('$h_1^{min}\;(\mu m)$','Interpreter','LaTex')
% yyaxis right
% plot(TBUTsigma,TBUT*t,'linewidth',2)
% ylabel('TBUT (s)','Interpreter','LaTex')
% xlabel('width\,(mm)','Interpreter','LaTex')

%% Fig 14 from Stapf 
%     TBUT_gaussian(j)=T(end);
%     TBUTH1_gaussian(j)=H1;
%     j=j+1;
% end
% 
% save TBUT_vs_H1_gaussian TBUT_gaussian TBUTH1_gaussian
% 
% figure();
% loglog(TBUTH1_gaussian*1e3,TBUT_gaussian*t,'linewidth',2)
% ylabel('TBUT (s)','Interpreter','LaTex')
% xlabel('$H_1\,(mm)$','Interpreter','LaTex')

%% experimenting with A/B values
% minAL(:,k)=min(hh1,[],2);
% AB(k)=A;
% k=k+1;
% % minLL(:,i)=min(hh2,[],2);
% % maxos(:,i)=max(cc,[],2);
% % maxlipid(:,i)=max(gg,[],2);
% end
% plot(AB,minAL,'o','linewidth',3)
% hold on
% % end
% xlabel('A/B')
% ylabel('min AL thickness')



%% Functions needed
%  Put 'em all in this file.
%% RHS function
function f = efnewt(t,w,y,v0,Upsilon,C1,C2,Pe1,Pe2,A,B,d,M,dd,Op,E,R,D1,D2,D3,setbc,nu,thetaB)
% [T,W] = ode15s(@(t,u) efnewt(t,u,),tspan,y0,MyTols);
% function for nonlinear evolution of two layer model with osmolarity
% and surfactant; local scaling
% Homogeneous Neumann conditions on both ends, all variables
% input:  t  (scalar, not always used)
%         w  row vector of h1, h2, p1, u2, m, g
%         y  independent variable (in equation of h after mapping)
%         v0, Upsilon, C1, C2, Pe1, Pe2, A, B, d, M, dd, Op, E, R  parameters
%         D1, D2  differentiation matrices
%         setbc   choose Neumann or Dirichlet bcs
%         nu      for mixed bcs 
% output: f rhs of f  (column vector)
N = length(D1);

h1 = w(1:N); h2 = w(N+1:2*N); p1 = w(2*N+1:3*N); % readable variables
u2 = w(3*N+1:4*N); m = w(4*N+1:5*N); g = w(5*N+1:6*N);
s = 1+v0*t;  % end location

Jo = Op*(m./h1 - 1);    % Osmosis
% Je = E./(1+R*h2); %no theta dependence
Je = E./(1+R*(.1+.9*sin(thetaB))*h2);       % Evaporation
% Je = E./(1+(10.1299*thetaB+1.768)*h2);        % linear function for R(thetaB)
h = h1 + d*h2;          % total film thickness

h1x = D1*h1/s;  h2x = D1*h2/s;  u2x = D1*u2/s;  % convenient variables
p1x = D1*p1/s;  mx = D1*m/s;    gx = D1*g/s;
hx = D1*h/s;

% depth average velocity in AL
%u1ave = -p1x.*h1.^2/3 + (Upsilon*C2*d*h2.*(D2*hx)/s^2 - Upsilon*(B/A)*d*(D1*(u2x.*h2))/s - M*gx).*(h1/2); 
u1ave = -p1x.*h1.^2/12 + u2/2; % correct form
u1ave(1) = 0;    % bc at left
u1ave(end) = v0; % bc at right

q1 = u1ave.*h1;    dq1 = D1*q1/s;  % flux and derivative
q2 = u2.*h2;       dq2 = D1*q2/s;

% h1, AL thickness: all points
f(1,1)           = h1x(1);
f(2:N-1,1)       = v0*y(2:end-1).*h1x(2:end-1) - dq1(2:end-1) + Jo(2:end-1) - Je(2:end-1);
f(N,1)           = h1x(end);

% h2, LL thickness: interior points
f(N+2:2*N-1,1)   = v0*y(2:end-1).*h2x(2:end-1) - dq2(2:end-1);

% p1, pressure in AL: all points for now
f(2*N+1:3*N,1)   = p1 + Upsilon*C2*(D1*hx)/s + C1*(D1*h1x)/s; 
% f(2*N+1,1)         = p1x(1);
% f(2*N+2:3*N-1,1)   = p1(2:end-1) + Upsilon*C2*(D1(2:end-1,:)*hx)/s + C1*(D1(2:end-1,:)*h1x)/s; 
% f(3*N,1)           = p1x(end);
            
% u2, axial velocity component in LL: all points
f(3*N+1,1)       = u2(1);          % fixed at origin
% f(3*N+2:4*N-1,1) = (A/B)*d*Upsilon*(D1(2:end-1,:)*(u2x.*h2))/s + u2(2:end-1)./h1(2:end-1) ...
%                     + h1(2:end-1).*p1x(2:end-1)./2 - C2*d*Upsilon*h2(2:end-1).*(D2(2:end-1,:)*hx)/s^2 ...
%                     + M*gx(2:end-1); 
f(3*N+2:4*N-1,1) = (A/B)*d*Upsilon*(D1(2:end-1,:)*(u2x.*h2))/s + u2(2:end-1)./h1(2:end-1) ...
                    + h1(2:end-1).*p1x(2:end-1)./2 - C2*d*Upsilon*h2(2:end-1).*(D3(2:end-1,:)*h)/s^3 ...
                    + M*gx(2:end-1); % using third derivative matrix
f(4*N,1)         = u2(end) - v0;     % right end moves

% m, mass=h1*c: all points for now
f(4*N+1,1)           = mx(1);
f(4*N+2:5*N-1,1)     = v0*y(2:end-1).*mx(2:end-1) + (1/Pe1)*(D1(2:end-1,:)*(mx - m.*h1x./h1))/s - D1(2:end-1,:)*(m.*u1ave)/s;
f(5*N,1)             = mx(end);

% g, surfactant: all points for now
f(5*N+1,1)           = gx(1);
f(5*N+2:6*N-1,1)     = v0*y(2:end-1).*gx(2:end-1) - D1(2:end-1,:)*(g.*u2)/s + (1/Pe2)*(D1(2:end-1,:)*gx)/s;
f(6*N,1)             = gx(end);

if setbc == 0
   f(N+1,1)       = h2x(1);         % h2x = 0 at left end
   f(2*N,1)       = h2x(end);       % h2x = 0 at right end
end
if setbc == 1 
   f(N+1,1)       = h2x(1) + dd;    % h2x nonzero at left end
   f(2*N,1)       = h2x(end) - dd;  % h2x nonzero at right end
end
if setbc == 2
   f(N+1,1)       = h2x(1);         % h2x = 0 at left end
   f(2*N,1)       = h2(end) - 1;    % h2 = 1 at right end
end
if setbc == 5   %robin at right
   f(N+1,1)       = h2x(1);         % h2x = 0 at left end
   f(2*N,1)       = nu*h2x(end)+(1-nu)*(h2(end) - 1);   % robin at right
end
if setbc == 6   %robin at left
   f(N+1,1)       = -nu*h2x(1)+(1-nu)*(h2(1) - 1);   % robin at left
   f(2*N,1)       = h2x(end);       % h2 = 0 at right end
end
end

%% event detection 
function [position, isterminal, direction] = zerothickness(w,N,epsilon,L)
    position = double(min(w(1:N+1))>=0.5/(epsilon*L*1e6)); %look at the minimum thickness of AL, and compare to 0.5 micrometers
%     position = w(35)-0.5;
    isterminal = 1; % halt integration if thickness is too small
    direction = 0; % zero can be approached from either side
end

%% Plot the ICs
function plotic(y,h10,h20,u20,p10,m0,g0)
  %plot u0 as a check
   figure
   plot(y,h10,'r--',y,h20,'b--',y,p10,'k--',y,u20,'b-',y,m0,'g:',y,g0,'b:','LineWidth',2)
   xlabel('x'), title('Initial Conditions')  
   hold on
   %plot(y,h0.*(D1*u0),'g-.','LineWidth',2)
   %legend('u_0','h_0','p_0','h_0 u_{0,x}','Location','East')
   legend('h1_0','h2_0','p1_0','u2_0','m_0','g_o','Location','East')
   title('Initial Conditions','FontSize',28,'Interpreter','LaTeX')
   hold off

end

%% Dimensional plots
% in order to compare to Stapf's results

function moving_dim(y,T,hh1,hh2,uu2,pp1,mm,gg,cc,v0,L,U,epsilon,d,mu1,Gamma0,C)
% first dimensionalize results

t_dim = (L/U)*T; % (s)
h1_dim = 1e6*epsilon*L*hh1; % (micrometers)
h2_dim = 1e9*epsilon*d*L*hh2; % (nanometers)
p1_dim = (mu1*U)/(epsilon^2*L)*pp1; % (Pa)
u2_dim = 6e7*U*uu2; % (micrometers/min)
g_dim = Gamma0*gg; % (mol/m^2)
c_dim = C*cc; % (mOsM)
y_dim=L*1000*(1+v0*T)*y'; % (mm)

% plot all dimensional quantities using subplot
h = figure('Name','results');
subplot(2,3,1)
plot(y_dim',h1_dim','LineWidth',2);
xlabel('$x \;(mm)$','FontSize',26,'Interpreter','LaTeX')
ylabel('$h_1 \; (micrometers)$','FontSize',26,'Interpreter','LaTeX')
ax=gca;
ax.FontSize = 18;
ylim([0 4])
xlim([-1 1])
title('$h1$ on moving domain','FontSize',18,'Interpreter','LaTeX')

subplot(2,3,2)
plot(y_dim',h2_dim','LineWidth',2);
xlabel('$x \;(mm)$','FontSize',26,'Interpreter','LaTeX')
ylabel('$h_2 \;(nm)$','FontSize',26,'Interpreter','LaTeX')
ax=gca;
ax.FontSize = 18;
ylim([0 60])
xlim([-1 1])
title('$h2$ on moving domain','FontSize',18,'Interpreter','LaTeX')


% subplot(2,3,3)
% plot(y_dim',u2_dim','LineWidth',2);
% xlabel('$x \;(mm)$','FontSize',26,'Interpreter','LaTeX')
% ylabel('$u_2 \;(m/s)$','FontSize',26,'Interpreter','LaTeX')
% ax=gca;
% ax.FontSize = 18;
% % ylim([0 60])
% xlim([-1 1])
% title('$u2$ on moving domain','FontSize',18,'Interpreter','LaTeX')

subplot(2,3,3)
plot(y_dim',uu2','LineWidth',2);
xlabel('$x \;(mm)$','FontSize',26,'Interpreter','LaTeX')
ylabel('$u_2$','FontSize',26,'Interpreter','LaTeX')
ax=gca;
ax.FontSize = 18;
ylim([-.06 0.06])
%xlim([-1 1])
title('$u2$ on moving domain','FontSize',18,'Interpreter','LaTeX')

subplot(2,3,4)
plot(y_dim',p1_dim','LineWidth',2);
xlabel('$x \;(mm)$','FontSize',26,'Interpreter','LaTeX')
ylabel('$p_1 \;(Pa)$','FontSize',26,'Interpreter','LaTeX')
ax=gca;
ax.FontSize = 18;
% ylim([0 60])
xlim([-1 1])
title('$p1$ on moving domain','FontSize',18,'Interpreter','LaTeX')

subplot(2,3,5)
plot(y_dim',c_dim','LineWidth',2);
xlabel('$x \;(mm)$','FontSize',26,'Interpreter','LaTeX')
ylabel('$c \;(mOsM)$','FontSize',26,'Interpreter','LaTeX')
ax=gca;
ax.FontSize = 18;
ylim([0 600])
xlim([-1 1])
title('$c$ on moving domain','FontSize',18,'Interpreter','LaTeX')
legendCell = cellstr(num2str(t_dim, 't=%.4f s\n'));
lgd=legend(legendCell,"location","south");
title(lgd,'Time','FontWeight','normal','FontSize',14)

% subplot(2,3,6)
% plot(y_dim',g_dim','LineWidth',2);
% xlabel('$x \;(mm)$','FontSize',26,'Interpreter','LaTeX')
% ylabel('$\Gamma \;(mol/m^2)$','FontSize',26,'Interpreter','LaTeX')
% ax=gca;
% ax.FontSize = 18;
% % ylim([.997 1.003])
% xlim([-1 1])
% title('$\Gamma$ on moving domain','FontSize',18,'Interpreter','LaTeX')

subplot(2,3,6)
plot(y_dim',gg','LineWidth',2);
xlabel('$x \;(mm)$','FontSize',26,'Interpreter','LaTeX')
ylabel('$\Gamma$','FontSize',26,'Interpreter','LaTeX')
ax=gca;
ax.FontSize = 18;
ylim([.997 1.003])
%xlim([-1 1])
title('$\Gamma$ on moving domain','FontSize',18,'Interpreter','LaTeX')

set(h, 'Position', [0 0 1500 1000])

% % individual plots for manuscript
fpath = '/Users/mjgocken/Desktop/sample figs nov 8/formatting test';

h = figure('Name','aqueous');
ax=gca;
set(ax,'FontSize',32,'Units','pixels','Position', [150 100 383 385])
mk = {'--','-','-.',':','--','-','-.'};
hold on
for j = 1:length(T) 
   plot(y_dim(j,:),h1_dim(j,:),'LineWidth',3.5,'LineStyle',mk{j})
end
xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
ylabel('$h_1$ ($\mu$m)','FontSize',32,'Interpreter','LaTeX')
ylim([0 4])
xlim([-1 1])
xticks([-1 -0.5 0 0.5 1])
set(h, 'Position', [0 0 550 500])
box on
hold off
fname = 'ALprofile_thetaPi12_wThetaDep';
% saveas(gcf,fullfile(fpath,fname),'epsc');
% saveas(gcf,fullfile(fpath,fname),'fig');

h = figure('Name','lipid');
ax=gca;
set(ax,'FontSize',32,'Units','pixels','Position', [150 100 383 385])
mk = {'--','-','-.',':','--','-','-.'};
hold on
for j = 1:length(T) 
   plot(y_dim(j,:),h2_dim(j,:),'LineWidth',3.5,'LineStyle',mk{j},'DisplayName',strcat('t = ',num2str(round(t_dim(j),2))))
end
xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
ylabel('$h_2$ (nm)','FontSize',32,'Interpreter','LaTeX')
ylim([0 60])
xlim([-1 1])
xticks([-1 -0.5 0 0.5 1])
set(h, 'Position', [0 0 550 500])
box on
% legendCell = cellstr(num2str(t_dim, 't=%.f s\n'));
% lgd=legend(legendCell,"location","southwest");
% lgd=legend('$t=0$','$t=10$','$t=20$','$t=30$','$t=40$','$t=50$','$t=60$','Interpreter','LaTex','FontSize',29,'Location','southwest')
lgd=legend('FontSize',29,'Interpreter','LaTeX',"location","southwest");
title(lgd,'Time (s)','FontWeight','normal','FontSize',29);
hold off
fname = 'LLprofile_thetaPi12_wThetaDep';
%saveas(gcf,fullfile(fpath,fname),'epsc');
% saveas(gcf,fullfile(fpath,fname),'fig');

h = figure('Name','velocity');
ax=gca;
set(ax,'FontSize',32,'Units','pixels','Position', [150 100 383 385])
mk = {'--','-','-.',':','--','-','-.'};
plot(y_dim(1,:)+4,u2_dim(j,:))
hold on
for j = 2:length(T) 
   plot(y_dim(j,:),u2_dim(j,:),'LineWidth',3.5,'LineStyle',mk{j})
end
xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
ylabel('$u_2$ ($\mu$m/min)','FontSize',32,'Interpreter','LaTeX')
% ylim([-0.1 0.1])
xlim([-1 1])
xticks([-1 -0.5 0 0.5 1])
set(h, 'Position', [0 0 550 500])
box on
hold off
fname = 'u2profile_thetaPi12';
% saveas(gcf,fullfile(fpath,fname),'epsc');
% saveas(gcf,fullfile(fpath,fname),'fig');

% h = figure('Name','velocity');
% ax=gca;
% ax.FontSize = 32;
% mk = {'--','-','-.',':','--','-','-.'};
% hold on
% for j = 2:length(T) 
%    plot(y_dim(j,:),uu2(j,:),'LineWidth',2.5,'LineStyle',mk{j})
% end
% xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
% ylabel('$u_2$','FontSize',32,'Interpreter','LaTeX')
% % ylim([-0.005 0.005])
% xlim([-1 1])
% xticks([-1 -0.5 0 0.5 1])
% set(h, 'Position', [0 0 500 500])
% box on
% hold off
% fname = 'u2profile_thetaPi12_wThetaDep';
%saveas(gcf,fullfile(fpath,fname),'epsc');
% saveas(gcf,fullfile(fpath,fname),'fig');

h = figure('Name','pressure');
ax=gca;
set(ax,'FontSize',32,'Units','pixels','Position', [150 100 383 385])
mk = {'--','-','-.',':','--','-','-.'};
hold on
for j = 1:length(T) 
   plot(y_dim(j,:),p1_dim(j,:),'LineWidth',3.5,'LineStyle',mk{j})
end
xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
ylabel('$p_1$ (Pa)','FontSize',32,'Interpreter','LaTeX')
%ylim([-1 0.2])
xlim([-1 1])
xticks([-1 -0.5 0 0.5 1])
% yticks([-1 -0.5 0 0.2])
set(h, 'Position', [0 0 550 500])
box on
hold off
fname = 'p1profile_thetaPi12_wThetaDep';
%saveas(gcf,fullfile(fpath,fname),'epsc');
% saveas(gcf,fullfile(fpath,fname),'fig');

h = figure('Name','osmolarity');
ax=gca;
set(ax,'FontSize',32,'Units','pixels','Position', [150 100 383 385])
mk = {'--','-','-.',':','--','-','-.'};
hold on
for j = 1:length(T) 
   plot(y_dim(j,:),c_dim(j,:),'LineWidth',3.5,'LineStyle',mk{j})
end
xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
ylabel('$c$ (mOsM)','FontSize',32,'Interpreter','LaTeX')
%ylim([0 600])
xlim([-1 1])
xticks([-1 -0.5 0 0.5 1])
set(h, 'Position', [0 0 550 500])
box on
hold off
fname = 'osmprofile_thetaPi12_wThetaDep';
%saveas(gcf,fullfile(fpath,fname),'epsc');
% saveas(gcf,fullfile(fpath,fname),'fig');

h = figure('Name','surfactant');
ax=gca;
set(ax,'FontSize',32,'Units','pixels','Position', [150 100 383 385])
mk = {'--','-','-.',':','--','-','-.'};
hold on
for j = 1:length(T) 
   plot(y_dim(j,:),gg(j,:),'LineWidth',3.5,'LineStyle',mk{j})
end
xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
ylabel('$\Gamma$','FontSize',32,'Interpreter','LaTeX')
%ylim([.998 1.002])
xlim([-1 1])
xticks([-1 -0.5 0 0.5 1])
set(h, 'Position', [0 0 550 500])
box on
hold off
fname = 'surfprofile_thetaPi12_wThetaDep';
%saveas(gcf,fullfile(fpath,fname),'epsc');
% saveas(gcf,fullfile(fpath,fname),'fig');


end

%% consolidated plots
function consolidated(y,T,hh1,hh2,uu2,pp1,mm,gg,cc,v0,L,U,epsilon,d,mu1,Gamma0,C,t,H1,H2,V)

% first dimensionalize results
t_dim = (L/U)*T; % (s)
h1_dim = 1e6*epsilon*L*hh1; % (micrometers)
h2_dim = 1e9*epsilon*d*L*hh2; % (nanometers)
p1_dim = (mu1*U)/(epsilon^2*L)*pp1; % (Pa)
u2_dim = 6e7*U*uu2; % (micrometers/min)
g_dim = Gamma0*gg; % (mol/m^2)
c_dim = C*cc; % (mOsM)
y_dim=L*1000*(1+v0*T)*y'; % (mm)
half=floor(length(y_dim)/2); % for plotting half on yyplot

% h1 and h2
figure();
ax=gca;
ax.FontSize = 32;
mk = {'--','-','-.',':','--','-','-.'};
yyaxis left
hold on
for j = 1:length(T) 
   plot(y_dim(j,1:half),h1_dim(j,1:half),'LineWidth',2.5,'LineStyle',mk{j},'Marker','none')
end
colororder('default');
xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
xticks([-1 -0.5 0 0.5 1])
ylabel('$h_1$ ($\mu$m)','FontSize',32,'Interpreter','LaTeX')
ylim([0 4])
yyaxis right
hold on
for j = 1:length(T) 
   plot(y_dim(j,half:end),h2_dim(j,half:end),'LineWidth',2.5,'LineStyle',mk{j},'Marker','none','DisplayName',strcat('t = ',num2str(round(t_dim(j),2))))
end
hold off
colororder('default');
xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
ylabel('$h_2$ (nm)','FontSize',32,'Interpreter','LaTeX')
ylim([0 60])
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

% osmolarity and surfactant
figure();
ax=gca;
ax.FontSize = 32;
mk = {'--','-','-.',':','--','-','-.'};
yyaxis left
hold on
for j = 1:length(T) 
   plot(y_dim(j,1:half),c_dim(j,1:half),'LineWidth',2.5,'LineStyle',mk{j},'Marker','none')
end
colororder('default');
xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
xticks([-1 -0.5 0 0.5 1])
ylabel('$c$ (mOsM)','FontSize',32,'Interpreter','LaTeX')
ylim([0 615])
yyaxis right
hold on
for j = 1:length(T) 
   plot(y_dim(j,half:end),gg(j,half:end),'LineWidth',2.5,'LineStyle',mk{j},'Marker','none','DisplayName',strcat('t = ',num2str(round(t_dim(j),2))))
end
hold off
colororder('default');
xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
ylabel('$\Gamma$','FontSize',32,'Interpreter','LaTeX')
% ylim([0.998 1.002])
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

% figure();
% hold on
% for j = 1:length(T)
%     plot(t_dim(j),min(h1_dim(j,:)),'o','LineWidth',2)
% end
% xlabel('$t$','FontSize',32,'Interpreter','LaTeX')
% ylabel('$h_{1min}$','FontSize',32,'Interpreter','LaTeX')
end

%% Surfactant IC
function [gamIC, min_gamIC, max_gamIC] = tanh_dist(y, w, mu, d)
	gamIC = .5+tanh(6*(y-w)/mu)/2 + 0.5+tanh(6*(-y-w)/mu)/2;
    min_gamIC = min(gamIC);
	gamIC = gamIC - min_gamIC;
    max_gamIC = max(gamIC);
	gamIC = gamIC / max_gamIC;
	gamIC = 1-d+d*gamIC;
	return
end
