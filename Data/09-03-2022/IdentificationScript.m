clear;

dir = "pos";
val = [0.06 0.065];

% dir = "neg";
% val = [0.04 0.045 0.0475 0.05 0.0525];

G = cell(2);

for m = 1:3
    figure;
    hold on;
    for n = 1:length(val)
        load("idf_" + dir +"_" + num2str(val(n)) + "_" + num2str(m) + ".mat");
        t = data(:, 1);
        w = data(:, 4);
        M = data(:, 5);
    
        idx = (t >= 9.5 & t < 10);
        MQ = mean(M(idx));
        wQ = mean(w(idx));
    
        idx = (t >= 9.5);
        t_idf = t(idx);
        t_idf = t_idf - t_idf(1);
        w_idf = w(idx) - wQ;
        M_idf = M(idx) - MQ;

%         N = 100;
%         t_idf = t_idf(1:end-N+1);
%         w_idf = moving_average(w_idf, N);
%         M_idf = M_idf(N/2:end-N/2);

        H = [-w_idf(1:end-1), M_idf(2:end) + M_idf(1:end-1)];
%         H = [w_idf(1:end-1), M_idf(2:end)];
        c{n} = (H'*H)\H'*w_idf(2:end);
%         Q = 1/2*(w_idf(2:end) - H*c{n})'*(w_idf(2:end) - H*c{n});
        T(m, n) = Tvz/2*(1 - c{n}(1))/(1 + c{n}(1));
        K(m, n) = (1 + 2*T(m, n)/Tvz)*c{n}(2);
%         T(m, n) = c{n}(1)*Tvz/(1 - c{n}(1));
%         K(m ,n) = c{n}(2)*(T(m, n)/Tvz + 1);

%         c = NLSM(t_idf-0.5, w_idf, M_idf);
%         K(m, n) = c(1);
%         T(m, n) = c(2);

        G = tf(K(m, n), [T(m, n) 1]);
%         Gz = tf(c{n}(2)*[1 1], [1 c{n}(1)], Tvz);
%         Gz = tf(c{n}(2)*[1], [1 -c{n}(1)], Tvz);

        data = iddata(w_idf, M_idf, Tvz);
        G = tfest(data, 1, 0);
        K(m, n) = G.Numerator/G.Denominator(2);
        T(m, n) = 1/G.Denominator(2);

        plot(t_idf, w_idf);
        plot(t_idf, lsim(G, M_idf, t_idf), 'r');
%         plot(t_idf, lsim(Gz, M_idf, t_idf), 'g');
%         plot(t_idf(2:end), H*c{n}, 'r');
    end
    hold off;
end

% mean(K)
% mean(K0)
% mean(T)
% mean(T0)

function [yf] = moving_average(y, dn)
    size_y = length(y);
    yf = zeros(size_y - dn, 1);
    n = 1;
    yf(1) = sum(y(1 : dn))/dn;

    while (n + dn <= size_y)
        yf(n + 1) = yf(n) - y(n)/dn + y(n + dn)/dn;
        n = n + 1;
    end
end

function c_best = NLSM(t, y, u)
    f = @(t, u, c) c(1)*u.*(1 - exp(-t/c(2)));
    Df = @(t, u, c) [u.*(1 - exp(-t/c(2))), -c(1)/c(2)^2*u.*t.*exp(-t/c(2))];

    c = [1/7.03e-5; 1.2e-4/7.03e-5];
    c = [1e5; 1.0];
    c_best = c;
    R = y - f(t, u, c);
    Q = sum(R.^2);
    Q_min = Q;
    stop = 0;

    for n = 1 : 1000
        H = Df(t, u, c);
        dc = (H'*H)\H'*R;
        c = c + dc;
        R = y - f(t, u, c);
        Q = sum(R.^2);
        
        if (Q < Q_min)
            Q_min = Q;
            c_best = c;
            stop = 0;
        end
        
        if (stop >= 100)
            break;
        else
            stop = stop + 1;
        end
    end
end