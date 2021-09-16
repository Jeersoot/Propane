clf
clc
clear

R = 8.31446; % Enhet: J/(K*mol)

%% Deluppgift 1

% Skapar en lista med all data vi har fått
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

% Skapar en lista med lika många färger som element vi har i data-listan
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

% 4x1 cell med 40x1-matriser
Z = {(lowest_temp(:,2).*lowest_temp(:,1)./(222*R))
    (under_critical_temp(:,2).*under_critical_temp(:,1)./(333*R))
    (over_critical_temp(:,2).*over_critical_temp(:,1)./(425*R))
    (highest_temp(:,2).*highest_temp(:,1)./(481*R))};

% Skapar en cell med respektive tryck
P = {(lowest_temp(:,1))
    (under_critical_temp(:,1))
    (over_critical_temp(:,1))
    (highest_temp(:,1))};

% Skapar en cell med resp. temperaturer
T = {(222) (333) (425) (481)};
 
% Plottar figur 2
figure(2)
for i = 1:4
    plot(P{i},Z{i})
    hold on
end

% Plottar figur 3
figure(3)
for i = 1:4
    plot(P{i},Z{i})
    if i == 1
        
    end
    hold on
end

% Finslipar 2:a och 3:e figuren
for n = 2:3
    legend('Lägsta temperatur (222 K)', 'Under kritisk temp. (333 K)', ...
        'Över kritisk temp. (425 K)', 'Högsta temperatur (481 K)')
    % Modifierar 3:e figuren
    if n == 3
        xlim([0 4*10^6])
    end
    title('ZP-diagram, Propan')
    ylabel('Z [-]')
    xlabel('P [Pa]')
end

%function [Z G] = PengRobinson(P,T,Pc,Tc,w):
%
%end
