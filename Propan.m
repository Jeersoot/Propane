clf
clc
clear

R = 8.31446; % Enhet: J/(K*mol)
Tc = 369.8; % Enhet: K
Pc = 4.249*10^6; % Enhet: Pa
w = 0.152; % Enhet: -
%% Deluppgift 1

% Skapar en cell med all data vi har fått
data = {(csvread('propane_JMV.dat'))
    (csvread('propane_PV222K.dat'))
    (csvread('propane_PV240K.dat'))
    (csvread('propane_PV259K.dat'))
    (csvread('propane_PV277K.dat'))
    (csvread('propane_PV296K.dat'))
    (csvread('propane_PV314K.dat'))
    (csvread('propane_PV333K.dat'))
    (csvread('propane_PV351K.dat'))
    (csvread('propane_PV362K.dat'))
    (csvread('propane_PV388K.dat'))
    (csvread('propane_PV407K.dat'))
    (csvread('propane_PV425K.dat'))
    (csvread('propane_PV444K.dat'))
    (csvread('propane_PV462K.dat'))
    (csvread('propane_PV481K.dat'))};

% Skapar en cell med lika många färger som element vi har i data-cellen
colours = {('rx'), (1/255*[255 165 0]), ('yellow'), ('green'), ('blue'), ...
    (1/255*[75 0 130]), ('magenta'), (1/255*[148 0 211]), ...
    (1/255*[255 0 255]), (1/255*[192 192 192]), (1/255*[255 215 0]), ...
    (1/255*[207 185 151]), (1/255*[128 0 0]), (1/255*[128 128 128]), ... 
    ('black'), (1/255*[64 224 208])};

figure(1)
% Skapar en plot och parar ihop rätt datasträng med rätt färg
for i = 1:length(data)
    if i == 1
        loglog(data{1}(:,3), data{1}(:,2), 'bx')
        hold on
        loglog(data{1}(:,4), data{1}(:,2), 'rx')
    else
        loglog(data{i}(:,2), data{i}(:,1), 'Color', colours{1,i})
    end
end

% Finslipar plottens utseende
ylim([3*10^4 10^7])
legend('Mättad vätska', 'Mättad gas', '222 K','240 K', '259 K', ...
    '277 K', '296 K', '314 K', '333 K', '351 K','362 K', '388 K', ...
    '407 K', '425 K', '444 K', '462 K', '481 K')
title('PV-diagram, Propan')
xlabel('V [m^{3}/mol]')
ylabel('P [Pa]')

%% Deluppgift 2
% De fyra temperaturerna som väljs
lowest_temp = data{2}; % 222 K
highest_temp = data{end}; % 481 K
under_critical_temp = data{8}; % 333 K
over_critical_temp = data{end-3}; % 425 K

% Sparar de fyra temperaturernas färger från figur 1
colour4temp = {(colours{1,2})
        (colours{1,8}) 
        (colours{1,end-3})
        (colours{1,end})};

% 4x1 cell med 40x1-matriser
Zaprx = {(lowest_temp(:,2).*lowest_temp(:,1)./(222*R))
    (under_critical_temp(:,2).*under_critical_temp(:,1)./(333*R))
    (over_critical_temp(:,2).*over_critical_temp(:,1)./(425*R))
    (highest_temp(:,2).*highest_temp(:,1)./(481*R))};

% Skapar en cell med respektive tryck
P = {(lowest_temp(:,1))
    (under_critical_temp(:,1))
    (over_critical_temp(:,1))
    (highest_temp(:,1))};

% Skapar en cell med resp. temperaturer
T = {(222)
    (333)
    (425)
    (481)};

% Skapar en cell med de olika B-värden
B = {(-0.00076)
    (-0.0003)
    (-0.00021)
    (-0.00013)};

for i = 1:4
    M = [];
    for j = 1:length(P{i})
        [Z,G] = PengRobinson(P{i}(j),T{i},Pc,Tc,w);
        M = [M; [Z, G']];
    end
    figure(2)
    plot(P{i},Zaprx{i}, 'Color', colour4temp{i})
    hold on
    plot(P{i},M(:,1),'o', 'Color', colour4temp{i})
    
    figure(3)
    plot(P{i},Zaprx{i}, 'Color', colour4temp{i})
    hold on
    plot(P{i}, 1 + B{i}.*P{i}./(R.*T{i}), '--', 'Color', colour4temp{i})
    plot(P{i},M(:,1),'o', 'Color', colour4temp{i})
end

% Finslipar 2:a och 3:e figuren
for n = 2:3
    figure(n)
    
    % Modifierar 3:e figuren
    if n == 3
        axis([0 4*10^6 0 1])
        legend('Lägsta temperatur (222 K)', 'Z-approx (222 K)', ...
        'Peng-Rob (222 K)','Under kritisk temp. (333 K)', ...
        'Z-approx (333 K)', 'Peng-Rob (333 K)', ...
        'Över kritisk temp. (425 K)', 'Z-approx (425 K)', ...
        'Peng-Rob (425 K)', 'Högsta temperatur (481 K)', ...
        'Z-approx (481 K)', 'Peng-Rob (481 K)')
    else
         legend('Lägsta temperatur (222 K)', 'Peng-Rob (222 K)', ... 
             'Under kritisk temp. (333 K)', 'Peng-Rob (333 K)', ...
            'Över kritisk temp. (425 K)', 'Peng-Rob (425 K)', ...
            'Högsta temperatur (481 K)', 'Peng-Rob (481 K)')
    end
    title('ZP-diagram, Propan')
    ylabel('Z [-]')
    xlabel('P [Pa]')
end

%% 4.c)

% 298 K
PV296 = data{6};
P296 = PV296(:,1);

M = [];
for i = 1:length(P296)
    [Z,G] = PengRobinson(P296(i),298.15,Pc,Tc,w);
    M = [M; [Z, G']];
end

% Plottar kolumnerna
figure(4)
for j = 2:4
    plot(P296,M(:,j))
    hold on
end
ylabel('(G-G^{ig})/(RT) [-]')
xlabel('P [Pa]')
legend('Z_{1}', 'Z_{2}', 'Z_{3}')

function [Z,G] = PengRobinson(P,T,Pc,Tc,w)
R = 8.31446; % Enhet: J/(K*mol)
Tr = T/Tc;

kappa = 0.37464 + 1.54226*w - 0.26992*(w^2); % Eq. (7.17)
alpha = (1 + kappa*(1 - sqrt(Tr)))^2; % Eq. (7.17)

ac = 0.45723553*((R*Tc)^2)/Pc; % Eq. (7.16)
a = ac*alpha; % Eq. (7.16)
A = a*P/(R^2*T^2); % Eq. (7.21)

b = 0.07779607*R*Tc/Pc; % Eq. (7.16)
B = b*P/(R*T); % Eq. (7.22)

p = B - 1;
q = A - 3*B^2 - 2*B;
r = B^3 + B^2 - A*B;
C = [1 p q r]; % C beskriver ekvationen: Z^3 + pZ^2 + qZ + r = 0
z = roots(C);

G = z-1-log(z-B)-(A/(B*sqrt(8)))*log((z+(1+sqrt(2))*B)./(z+(1-sqrt(2))*B));
if isreal(z)
    g = G;
else
    g = abs(imag(z));
end
[g, index] = sort(g);
Z = z(index(1));
end