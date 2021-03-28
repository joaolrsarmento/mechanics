% Instituto Tecnológico de Aeronáutica
% Engenharia Mecânica-Aeronáutica
% MPD-11 - Dinâmica de máquinas
% Aluno: João Sarmento
% 
% Laboratório 1 - Cinemática de mecanismos
% Dados experimentais - Jugo escocês
% 
%% Inicialização

clear;
close all;
clc;

%% Medidas

% Comprimento do elo 2 (valor fixo, necessários para comparação
% de cálculos analíticos/numéricos com dados experimentais):

r_2 = 0.09150 - 0.125*0.0254 - 0.375*0.0254;    % [m]

% Valores de theta_2 considerados para medição do deslocamento do elo 4:

theta_2 = 0:5:360;                              % [graus]

% O valor nulo corresponde à situação em que o elo 2 está posicionado na
% horizontal (para a direita). Os valores aumentam no sentido anti-horário

% Medidas do deslocamento do elo 4:

r_4 = [
    157.2
    157.0
    156.0
    154.6
    152.4
    149.4
    146.4
    142.6
    138.4
    133.4
    128.4
    122.8
    117.0
    111.0
    104.6
    98.2
    92.4
    84.4
    77.4
    70.7
    64.0
    56.8
    50.4
    44.4
    38.2
    32.6
    27.0
    22.2
    17.4
    13.2
    9.4
    6.4
    4.0
    2.0
    0.7
    0.0
    -0.4
    0.1
    0.9
    2.6
    4.6
    7.4
    10.8
    14.6
    19.0
    23.6
    28.8
    34.5
    40.4
    46.4
    52.8
    59.3
    66.2
    72.8
    80.0
    86.6
    93.4
    100.3
    106.6
    113.2
    119.0
    124.6
    130.4
    135.2
    139.6
    144.2
    147.5
    150.5
    153.0
    155.0
    156.4
    157.0
    157.2
    ]/1000;                                     % [m]

% Os valores assumidos por r_4 aumentam para a direita (conforme fotos e 
% vídeo relacionados ao laboratório).

%% Cálculo analítico da posição
grau_para_rad = pi/180;
vetor_theta_2 = theta_2 * grau_para_rad;
r_4_analitico = 2*r_2 - r_2 .* (1 - cos(vetor_theta_2));

%% Cálculo do erro
erro_posicao = abs(r_4_analitico' - r_4);
erro_posicao_relativo = abs(erro_posicao ./ r_4_analitico');

%% Gráfico da posicao

figure;
plot(theta_2,r_4*1e3,'-o', 'Color', 'Black');
hold on
plot(theta_2, r_4_analitico*1e3, 'Color', 'Black', 'LineWidth', 2);
hold off
legend('Experimental', 'Teórico');
xlabel('$\theta_{2}$ [$^{\circ}$]','Interpreter','latex');
ylabel('$r_{4}$ [mm]','Interpreter','latex');
title('\bf Jugo Escoces','Interpreter','latex');
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
v_analitico = -w_2  * r_2 * sin(vetor_theta_2);

%% Cálculo do erro
erro_velocidade = abs(v_analitico' - v_numerico);
erro_velocidade_relativo = abs(erro_velocidade ./ v_analitico');

%% Gráficos da velocidade angular
figure;
plot(t,v_analitico,'LineWidth', 2, 'Color', 'Black');
hold on
plot(t,v_numerico, '-o', 'Color', 'Black');
hold off
legend('Teórico', 'Númerico');
xlabel('t(s)','Interpreter','latex');
ylabel('v[m/s]','Interpreter','latex');
title('\bf Jugo Escoces','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

%% Cálculo númerico da aceleração angular
a_numerico = zeros(length(theta_2));
t = 0:dt:(length(theta_2)-1)*dt;
for i = 1:length(theta_2)
    if i == length(theta_2)
        a_numerico(i) = a_numerico(i-1);
    else
        if i == 1
            a_numerico(i) = (v_numerico(i+2)-v_numerico(i)) ./ (2*dt) ;
        else
            a_numerico(i) = (v_numerico(i+1) - v_numerico(i-1)) ./ (2*dt);
        end
    end
end
a_numerico = a_numerico(:, 1);

%% Cálculo analítico da aceleração angular
a_analitico = w_2 ^ 2 * r_2 .* cos(pi-vetor_theta_2);

%% Cálculo do erro
erro_aceleracao = abs(a_analitico' - a_numerico);
erro_aceleracao_relativo = abs(erro_aceleracao ./ a_analitico');

%% Gráficos da aceleração angular
figure;
plot(t,a_analitico,'LineWidth', 2, 'Color', 'Black');
hold on
plot(t,a_numerico, '-o', 'Color', 'Black');
hold off
legend('Teórico', 'Númerico');
xlabel('t(s)','Interpreter','latex');
ylabel('$a$ [$m/s^2$]','Interpreter','latex');
title('\bf Jugo Escoces','Interpreter','latex');
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