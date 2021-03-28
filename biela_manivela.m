% Instituto Tecnológico de Aeronáutica
% Engenharia Mecânica-Aeronáutica
% MPD-11 - Dinâmica de máquinas
% Aluno: João Sarmento
% 
% Laboratório 1 - Cinemática de mecanismos
% Dados experimentais - Biela-manivela
% 
%% Inicialização

clear;
close all;
clc;

%% Medidas

% Comprimentos dos elos 2 e 3 (valores fixos, necessários para comparação
% de cálculos analíticos/numéricos com dados experimentais):

r_2 = 0.09150 - 0.125*0.0254 - 0.375*0.0254;    % [m]
r_3 = 0.195;                                    % [m]

% Valores de theta_2 considerados para medição do deslocamento do 'pistão'
% (elo 4):

theta_2 = 0:5:360;                              % [graus]

% O valor nulo corresponde à situação em que o elo 2 está posicionado na
% horizontal (para a direita). Os valores aumentam no sentido anti-horário

% Medidas do deslocamento do pistão (elo 4):

r_4 = [
    158
    157.8
    157.2
    156.2
    155.0
    153.4
    151.2
    148.8
    146.0
    142.8
    139.2
    134.8
    130.4
    125.6
    120.4
    114.6
    108.4
    102.0
    95.4
    88.3
    81.2
    73.8
    66.1
    58.8
    51.6
    44.6
    37.4
    31.2
    25.2
    19.6
    14.4
    10.0
    6.4
    3.6
    1.4
    0.2
    -0.2
    0.0
    1.2
    3.4
    6.2
    9.8
    14.0
    19.2
    24.6
    30.8
    37.4
    44.0
    51.2
    58.8
    66.2
    73.6
    81.0
    88.2
    95.0
    102.0
    108.2
    114.6
    120.0
    125.6
    130.4
    135.0
    139.0
    142.8
    146.0
    148.8
    151.4
    153.2
    154.8
    156.2
    157.0
    157.6
    157.8
    ]/1000;                                     % [m]

% Os valores assumidos por r_4 aumentam para a direita (conforme fotos e 
% vídeo relacionados ao laboratório).

%% Cálculo analítico da posição
grau_para_rad = pi / 180;
vetor_theta_2 = theta_2 * grau_para_rad;
x = (r_2 .* (1 - cos(pi-vetor_theta_2)) + (r_2 ^ 2 / (2 * r_3 )) .* sin(pi-vetor_theta_2) .^ 2);


%% Cálculo do erro
erro_posicao = abs(x'-r_4);
erro_posicao_relativo = abs(erro_posicao./x');
%% Gráfico da posição
figure;
plot(theta_2,r_4*1e3,'-o', 'Color', 'Black');
hold on
plot(theta_2,x*1e3, 'LineWidth', 2, 'Color', 'Black');
hold off
legend('Experimental', 'Teórico');
xlabel('$\theta_{2}$ [$^{\circ}$]','Interpreter','latex');
ylabel('$r_{4}$ [mm]','Interpreter','latex');
title('\bf Biela-Manivela','Interpreter','latex');
xlim([0,360]);
xticks(0:30:360);
xtickformat('$%g$');
% ylim([-5,160]);
% yticks(0:10:160);
ytickformat('$%g$');
set(gca,'TickLabelInterpreter','latex');

%% Calculo númerico da velocidade

dt = 0.001; % assume-se um delta_t de 1ms
v_numerico = zeros(length(theta_2));
t = 0:dt:(length(theta_2)-1)*dt;
for i = 1:length(theta_2)
    if i == length(theta_2)
        v_numerico(i) = v_numerico(i-1);
    else
        if i == 1
            v_numerico(i) = (r_4(i+2)-r_4(i)) ./ (2 * dt);
        else
            v_numerico(i) = (r_4(i+1)-r_4(i-1)) ./ (2 * dt);
        end
    end
end
v_numerico = v_numerico(:, 1);
%% Cálculo analítico da velocidade 
vetor_theta_2 = theta_2 * grau_para_rad;
w_2 = 5 * grau_para_rad / dt;
v_analitico = -w_2 * r_2 .* (sin(vetor_theta_2) + r_2 / (2 * r_3) .* sin(2* pi - 2 .* vetor_theta_2));

%% Cálculo do erro
erro_velocidade = abs(v_analitico' - v_numerico);
erro_velocidade_relativo = abs(erro_velocidade./v_analitico');

%% Gráficos da velocidade angular
figure;
plot(t,v_analitico,'LineWidth', 2, 'Color', 'Black');
hold on
plot(t,v_numerico, '-o', 'Color', 'Black');
hold off
legend('Teórico', 'Númerico');
xlabel('t(s)','Interpreter','latex');
ylabel('v[m/s]','Interpreter','latex');
title('\bf Biela-Manivela','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

%% Cálculo númerico da aceleração angular
a_numerico = zeros(length(theta_2));
t = 0:dt:(length(theta_2)-1)*dt;
for i = 1:length(theta_2)
    if i == length(theta_2)
        a_numerico(i) = a_numerico(i-1);
    else
        if i == 1
            a_numerico(i) = (v_numerico(i+2) - v_numerico(i)) ./ (dt * 2);
        else
            a_numerico(i) = (v_numerico(i+1) - v_numerico(i-1)) ./ (dt * 2);
        end
    end
end
a_numerico = a_numerico(:, 1);

%% Cálculo analítico da aceleração angular
a_analitico = w_2 ^ 2 * r_2 .* (cos(pi - vetor_theta_2) + (r_2/r_3) * cos(2 * pi - 2 .* vetor_theta_2));

%% Cálculo do erro
erro_aceleracao = abs(a_analitico' - a_numerico);
erro_aceleracao_relativo = abs(erro_aceleracao./a_analitico');
%% Gráficos da aceleração angular
figure;
plot(t,a_analitico,'LineWidth', 2, 'Color', 'Black');
hold on
plot(t,a_numerico, '-o', 'Color', 'Black');
hold off
legend('Teórico', 'Númerico');
xlabel('t(s)','Interpreter','latex');
ylabel('$a$ [$m/s^2$]','Interpreter','latex');
title('\bf Biela-Manivela','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

%% Gráficos erros

figure;
histogram(erro_posicao, 'FaceColor', 'Black');
title('Erro Posicao (m): Frequencia x Valor');

figure;
histogram(erro_velocidade, 'FaceColor', 'Black');
title('Erro Velocidade (m/s): Frequencia x Valor');

figure;
histogram(erro_aceleracao, 'FaceColor', 'Black');
title('Erro Aceleracao (m/s2): Frequencia x Valor');

media_erro_posicao = mean(erro_posicao) 
media_erro_velocidade = mean(erro_velocidade)
media_erro_aceleracao = mean(erro_aceleracao)

std_posicao = std(erro_posicao) 
std_erro_velocidade = std(erro_velocidade)
std_erro_aceleracao = std(erro_aceleracao)