clear; close all;
a0 = 0;                            
a1 = 0;                           
a2 = 0;                           
K = 0;                             
beta = 0;                          

mdl = 'model_2';
open_system(mdl);
a = 9; b = 19; c = 33;             
delta = 0.02;                      
P = [-a, -a, -a;                   
     -a, -b, -b;                   
     -a, -c, -c;                   
     -a, -b+1i*c, -b-1i*c;         
     -c, -a+1i*a, -a-1i*a;         
     -c, -a+1i*c, -a-1i*c];        
R = zeros(6, 6);                   % массив результатов: eta, mu, tp*, sigma*, tp, sigma для 6 экспериментов
opts = {'Solver', 'ode45', 'StopTime', '2', 'MaxStep', '1e-4', 'RelTol', '1e-6', ...  % параметры расчета: решатель ode45,
'ReturnWorkspaceOutputs', 'on'};   % время 2 с, макс. шаг 1e-4, точность 1e-6, результат возвращается в переменную out
for i = 1:6                        % цикл по шести наборам полюсов
    p = P(i, :);                   % текущий набор полюсов (строка матрицы P)
    cf = real(poly(p));            % коэффициенты полинома [1 a2 a1 a0]; real убирает нулевую мнимую часть
    a2 = cf(2); a1 = cf(3); a0 = cf(4);   % записать коэффициенты в переменные, которые читает модель
    K = abs(prod(p));              % числитель |l1*l2*l3| (равен a0, поэтому W(0) = 1)
    out = sim(mdl, opts{:});       % запуск моделирования с заданными опциями
    hs = out.get('h');             % получить сигнал переходной характеристики h из результатов
    t = hs.Time; h = hs.Data;      % вектор времени и вектор значений h(t)
    W = zpk([], p, K);             % передаточная функция: нулей нет, полюса p, усиление K
    hy = dcgain(W);                % установившееся значение переходной характеристики (= 1)
    eta = min(abs(real(p)));       % степень устойчивости: расстояние до мнимой оси ближайшего полюса
    mu = max(abs(imag(p) ./ real(p)));   % степень колебательности: максимум |Im/Re| по полюсам
    tp_est = log(1/delta) / eta;   % оценка времени переходного процесса по eta
    s_est = 0;                     % оценка перерегулирования по умолчанию 0
    if mu > 0, s_est = 100 * exp(-pi/mu); end   % при колебательности оценка sigma* = 100*exp(-pi/mu), %
    e = abs(h - hy) - delta*abs(hy);   % отклонение от установившегося значения минус ширина трубки; e > 0 - вне трубки
    k = find(e > 0, 1, 'last');    % индекс последней точки, лежащей вне трубки
    tp = t(k) + (t(k+1) - t(k)) * e(k) / (e(k) - e(k+1));   % линейная интерполяция момента входа в трубку (шаг переменный)
    sigma = max(0, (max(h) - hy) / abs(hy) * 100);   % перерегулирование в %, не меньше нуля
    R(i, :) = [eta, mu, tp_est, s_est, tp, sigma];   % записать результаты эксперимента в строку i
    S = stepinfo(W, 'SettlingTimeThreshold', delta);   % контрольный расчет показателей по передаточной функции
    fprintf('эксп. %d: tp = %.4f, sigma = %.2f%%  | stepinfo: tp = %.4f, sigma = %.2f%%\n', ...   
        i, tp, sigma, S.SettlingTime, S.Overshoot);   

    % отдельный рисунок для каждого эксперимента
    figure('Position', [100 100 800 500]); hold on; grid on;
    h1 = plot(t, h, 'LineWidth', 1.2, 'DisplayName', sprintf('h(t): \\sigma = %.2f%%', sigma));
    h2 = plot(t([1 end]), hy*(1 + delta)*[1 1], 'k--', 'DisplayName', sprintf('\\pm%g%% от h_{уст}', delta*100));
    plot(t([1 end]), hy*(1 - delta)*[1 1], 'k--', 'HandleVisibility', 'off');
    h3 = plot([tp tp], [0 hy*(1 + delta)], 'r:', 'LineWidth', 1.2, 'DisplayName', sprintf('t_п = %.3f c', tp));
    xlim([0 1.2]); xlabel('t, c'); ylabel('h(t)');
    ylim([0 1.1*max(max(h), 1)]);
    legend('Location', 'best');
    % print(sprintf('task2_exp%d.png', i), '-dpng', '-r150');  
end