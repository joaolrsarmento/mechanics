% Instituto Tecnológico de Aeronáutica
% Engenharia Mecânica-Aeronáutica
% MPD-11 - Dinâmica de máquinas
% Aluno: João Sarmento
% 
% Laboratório 1 - Cinemática de mecanismos
% Dados experimentais - Quatro barras
% 
%% Inicialização

clear;
close all;
clc;

%% Medidas da Posição

% Comprimentos dos elos 1 a 4 (valores fixos, necessários para comparação
% de cálculos analíticos/numéricos com dados experimentais):

r_1 = 0.255;                        % [m]
r_2 = 0.06860 + 0.00790;            % [m]
r_3 = 0.280;                        % [m]
r_4 = 0.255;                        % [m]

% Valores de theta_2 considerados para medição do deslocamento angular 
% do elo 4:

theta_2 = 0:5:360;                  % [graus]

% O valor nulo corresponde à situação em que o elo 2 está posicionado na
% horizontal (para a direita). Os valores aumentam no sentido anti-horário

% Medidas do deslocamento angular do elo 4:

theta_4 = [
    14.6
    12.6
    10.6
    9.0
    7.5
    6.2
    5.3
    4.6
    4.2
    4.1
    4.1
    4.4
    4.8
    5.6
    6.4
    7.4
    8.4
    9.8
    11.2
    12.5
    14.0
    15.5
    17.0
    18.6
    20.2
    21.8
    23.4
    25.1
    26.5
    28.1
    29.6
    31.2
    32.5
    34.0
    35.2
    36.5
    37.7
    38.9
    39.9
    40.8
    41.7
    42.6
    43.4
    44.0
    44.5
    45.0
    45.3
    45.6
    45.8
    45.8
    45.8
    45.7
    45.5
    45.2
    44.7
    44.2
    43.4
    42.6
    41.5
    40.4
    39.1
    37.6
    36.1
    34.4
    32.5
    30.5
    28.3
    26.2
    24.0
    21.6
    19.2
    16.8
    14.6
    ];                              % [graus]

% O valor nulo corresponde à situação em que o elo 4 está posicionado na
% vertical (para cima). Os valores aumentam no sentido anti-horário


%% Cálculo númerico da posição
grau_para_rad = pi / 180;
vetor_theta_2 = theta_2 * grau_para_rad;
vetor_theta_3 = zeros(size(vetor_theta_2));
vetor_theta_4 = zeros(size(vetor_theta_2));
x0 = [45 * grau_para_rad; 60 * grau_para_rad];

for i = 1:length(vetor_theta_2)
    theta_2_it = vetor_theta_2(i);
    theta_1_it = 0;
    if i > 1
        x0 = [vetor_theta_3(i-1); vetor_theta_4(i-1)];
    end
    c = [r_1 r_2 r_3 r_4 theta_1_it theta_2_it]; % valor constante para cálculo da função f do mecanismo
    x = newton_raphson(@f_quatro_barras, @j_quatro_barras, x0, c); % método de newton-raphson para cada c
    vetor_theta_3(i) = x(1);
    vetor_theta_4(i) = x(2);
end

theta_4_numerico = vetor_theta_4 / grau_para_rad - 90;
theta_3_numerico = vetor_theta_3 / grau_para_rad;

%% Cálculo de erros

erro_posicao = abs(theta_4_numerico' - theta_4);
erro_posicao_relativo = abs(erro_posicao./theta_4_numerico');
%% Gráficos da posição

figure(1);
plot(theta_2,theta_4, '-o', 'Color', 'Black');
hold on
plot(theta_2,theta_4_numerico, 'LineWidth', 2, 'Color', 'Black');
hold off
legend('Experimental', 'Teórico');
xlabel('$\theta_{2}$ [$^{\circ}$]','Interpreter','latex');
ylabel('$\theta_{4}$ [$^{\circ}$]','Interpreter','latex');
title('\bf Quatro Barras','Interpreter','latex');
xlim([0,360]);
xticks(0:30:360);
xtickformat('$%g$');
% ylim([0,50]);
% yticks(0:5:50);
ytickformat('$%g$');
set(gca,'TickLabelInterpreter','latex');

%% Calculo númerico da velocidade angular

dt = 0.001; % assume-se um delta_t de 1ms
w_4_numerico = zeros(length(theta_2));
t = 0:dt:(length(theta_2)-1)*dt;
for i = 1:length(theta_2)
    if i == length(theta_2)
        w_4_numerico(i) = w_4_numerico(i-1);
    else
        if i == 1
            w_4_numerico(i) = (theta_4(i+2)-theta_4(i)) .* grau_para_rad ./ (2 * dt);
        else
            w_4_numerico(i) = (theta_4(i+1)-theta_4(i-1)) .* grau_para_rad ./ (2 * dt);
        end
    end
end
w_4_numerico = w_4_numerico(:, 1);
%% Cálculo analítico da velocidade angular
vetor_theta_2 = theta_2 * grau_para_rad;
w_2 = 5 * grau_para_rad / dt;
w_3_analitico = r_2 * w_2 .* sin(vetor_theta_4 - vetor_theta_2)./(r_3 .* sin(vetor_theta_3 - vetor_theta_4));
w_4_analitico = r_2 * w_2 .* sin(vetor_theta_2 - vetor_theta_3)./(r_4 .* sin(vetor_theta_4 - vetor_theta_3));
%% Cálculo erros
erro_velocidade = abs(w_4_analitico' - w_4_numerico);
erro_velocidade_relativo = abs(erro_velocidade./w_4_analitico');
%% Gráficos da velocidade angular
figure;
plot(t,w_4_analitico,'LineWidth', 2, 'Color', 'Black');
hold on
plot(t,w_4_numerico, '-o', 'Color', 'Black');
hold off
legend('Teórico', 'Númerico');
xlabel('t(s)','Interpreter','latex');
ylabel('$\omega_{4}$ [rad/s]','Interpreter','latex');
title('\bf Quatro Barras','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

%% Cálculo númerico da aceleração angular
a_4_numerico = zeros(length(theta_2));
t = 0:dt:(length(theta_2)-1)*dt;
for i = 1:length(theta_2)
    if i == length(theta_2)
        a_4_numerico(i) = a_4_numerico(i-1);
    else
        if i == 1
            a_4_numerico(i) = (w_4_numerico(i+2)-w_4_numerico(i+1)) ./ (2*dt);
        else
            a_4_numerico(i) = (w_4_numerico(i+1) - w_4_numerico(i-1)) ./ (2*dt);
        end
    end
end
a_4_numerico = a_4_numerico(:, 1);

%% Cálculo númerico da aceleração angular pelo método das dif finitas de 2 ordem
a_4_numerico_segunda = zeros(length(theta_2));
t = 0:dt:(length(theta_2)-1)*dt;
for i = 1:length(theta_2)
    if i == length(theta_2)
        a_4_numerico_segunda(i) = a_4_numerico_segunda(i-1);
    else
        if i == 1
            a_4_numerico_segunda(i) = (theta_4(i+2)-2*theta_4(i+1)+theta_4(i)) .* grau_para_rad./ dt ^ 2;
        else
            a_4_numerico_segunda(i) = (theta_4(i+1)-2*theta_4(i)+theta_4(i-1)) .* grau_para_rad ./ dt ^ 2;
        end
    end
end
a_4_numerico_segunda = a_4_numerico_segunda(:, 1);
%% Cálculo analítico da aceleração angular
A = r_4 .* sin(vetor_theta_4);
B = r_3 .* sin(vetor_theta_3);
C = r_2 * w_2 ^ 2 * cos(vetor_theta_2) + ...
    r_3 .* w_3_analitico .^ 2 .* cos(vetor_theta_3) - ...
    r_4 .* w_4_analitico .^ 2 .* cos(vetor_theta_4);
D = r_4 .* cos(vetor_theta_4);
E = r_3 .* cos(vetor_theta_3);
F = -r_2 .* w_2 ^ 2 .* sin(vetor_theta_2) - ...
    r_3 .* w_3_analitico .^ 2 .* sin(vetor_theta_3) + ...
    r_4 .* w_4_analitico .^ 2 .* sin(vetor_theta_4);

a_4_analitico = (C .* E - B .* F)' ./ (A .* E - B .* D);
a_4_analitico = a_4_analitico(:, 1);
%% Cálculo de erros
erro_aceleracao = abs(a_4_analitico - a_4_numerico);
erro_aceleracao_relativo = abs(erro_aceleracao./a_4_analitico);

erro_aceleracao_2 = abs(a_4_analitico - a_4_numerico_segunda);
erro_aceleracao_2_relativo = abs(erro_aceleracao_2./a_4_analitico);

%% Gráficos da aceleração angular
figure;
plot(t,a_4_analitico,'LineWidth', 2, 'Color', 'Black');
hold on
plot(t,a_4_numerico, '-o', 'Color', 'Black');
hold off
legend('Teórico', 'Númerico');
xlabel('t(s)','Interpreter','latex');
ylabel('$\alpha_{4}$ [rad/$s^2$]','Interpreter','latex');
title('\bf Quatro Barras','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

figure;
plot(t,a_4_analitico,'LineWidth', 2, 'Color', 'Black');
hold on
plot(t,a_4_numerico_segunda, '-o', 'Color', 'Black');
hold off
legend('Teórico', 'Númerico - 2 ordem');
xlabel('t(s)','Interpreter','latex');
ylabel('$\alpha_{4}$ [rad/$s^2$]','Interpreter','latex');
title('\bf Quatro Barras','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');


%% Gráficos erros

figure;
histogram(erro_posicao, 'FaceColor', 'Black');
title('Erro Posicao (graus): Frequencia x Valor');

figure;
histogram(erro_velocidade, 'FaceColor', 'Black');
title('Erro Velocidade (rad/s): Frequencia x Valor');

figure;
histogram(erro_aceleracao, 'FaceColor', 'Black');
title('Erro Aceleracao (rad/s2): Frequencia x Valor');

figure;
histogram(erro_aceleracao_2, 'FaceColor', 'Black');
title('Erro Aceleracao 2 (rad/s2): Frequencia x Valor');


media_erro_posicao = mean(erro_posicao) * grau_para_rad
media_erro_velocidade = mean(erro_velocidade)
media_erro_aceleracao = mean(erro_aceleracao)
media_erro_aceleracao_2 = mean(erro_aceleracao_2)

std_posicao = std(erro_posicao) * grau_para_rad
std_erro_velocidade = std(erro_velocidade)
std_erro_aceleracao = std(erro_aceleracao)
std_erro_aceleracao_2 = std(erro_aceleracao_2)