clc; clear; close all

g_x = [1 0 1 0 0 1 1 0 1 0 1 1 1 1 0 0 1];
k = 8;                  % Длина сообщения
r = length(g_x) - 1;    % Длина CRC части
num_messages = 2^k;     % Количество возможных сообщений
[H, G] = hammgen(4);    % Проверочная и порождающая матрицы

messages = de2bi((0:num_messages-1)');     % Все возможные сообщения
codewords = zeros(num_messages, k+r);    % Все кодовые слова с CRC
mes_3_11 = zeros(num_messages, 3, 11);     % Все кодовые слова, разделенные 
                                           % на три части по 8 бит + 3 нуля
ham_3_15 = zeros(num_messages, 3, 15);     % Закодированные по Хэммингу
channel_3_12 = zeros(num_messages, 3, 12); % Замодулированные BPSK (-1, 1)
                                           % без добавленных 3х нулей,
% Генерация                                  которые отправляются в канал
for j=1:num_messages
    [~, c] = gfdeconv([zeros(1, r), messages(j,:)], g_x);   % CRC часть
    codewords(j,:) = xor([zeros(1, r), messages(j,:)], ...  % Сообщение +
                         [c, zeros(1, k+r - length(c))]);   % + CRC
    mes_3_11(j, :, :) = [reshape(codewords(j, :), 8, 3)' ...% На 3 части +
                         [0 0 0; 0 0 0; 0 0 0]];            % + нули
    for i=1:3
        ham_3_15(j, i, :) = encode(reshape(mes_3_11(j, i, :), 1, 11), ...
                            15, 11, 'hamming/binary'); %mod(reshape(mes_3_11(j, i, :), 1, 11) * G, 2);
        channel_3_12(j, i, :) = ham_3_15(j, i, 1:end-3)*(-2) + 1; 
    end
end
% messages(5, :)
% codewords(5, :)
% disp(reshape(mes_3_11(5, 1, :), 1, 11))
% disp(reshape(mes_3_11(5, 2, :), 1, 11))
% disp(reshape(mes_3_11(5, 3, :), 1, 11))
% disp(reshape(ham_3_15(5, 1, :), 1, 15))
% disp(reshape(ham_3_15(5, 2, :), 1, 15))
% disp(reshape(ham_3_15(5, 3, :), 1, 15))
% disp(reshape(channel_3_12(5, 1, :), 1, 12))
% disp(reshape(channel_3_12(5, 2, :), 1, 12))
% disp(reshape(channel_3_12(5, 3, :), 1, 12))

% Теоретический расчет
SNRdB = -10:0;
SNR = 10.^(SNRdB./10);
Pe_bit_theor = qfunc(sqrt(2.*SNR));         
Ped_theor_up = ones(1, length(SNR)).*1/2^r; 
Ped_theor = zeros(1,length(SNRdB));
w = sum(codewords(2:end, :)');
d = min(w);
for i=1:length(SNR)
    Pe_theoretical = 0;
    for j=d:(k + r)
        Pe_theoretical = Pe_theoretical + sum(w == j) * ...
                (Pe_bit_theor(i)^j) * ((1 - Pe_bit_theor(i))^((k + r)-j));
    end
    Ped_theor(i) = Pe_theoretical;
end

% Моделирование
Pe_bit = zeros(1,length(SNRdB));  % Вероятность ошибки на бит
Pe_bit_without_Hamming = zeros(1,length(SNRdB));
Ped = zeros(1,length(SNRdB));     % Вероятность ошибки декодирования
Ped_soft = zeros(1,length(SNRdB));% Вероятность ошибки декодирования soft
T = zeros(1,length(SNRdB));       % Пропускная способность канала
for i=1:length(SNR)
    disp([i]);
    sigma = sqrt(1/(2*SNR(i)));
    nTests = 0; nSent = 0; 
    nErrDecode = 0; nErrBits = 0; nErrBits_withoutHamming = 0;
    for messages_to_send = 1:100*(SNRdB(i)+11)
        rand_ind = randi(num_messages, 1);
        c = reshape(channel_3_12(rand_ind, :, :), 3, 12);
        while 1
            nTests = nTests + 1;
            AWGN = c + sigma*randn(3, 12);
            unBPSK = AWGN < 0;
            nErrBits = sum(sum(xor(c, unBPSK)));
            plus_zeros = [unBPSK(:, :) [0 0 0; 0 0 0; 0 0 0]];
            
            corrected = zeros(3, 11);
            for part=1:3
                corrected(part, :) = decode(reshape(plus_zeros(part, :), 1, 15), 15, 11, ...
                                            'hamming/binary');%correctError(plus_zeros, H);
            end
            
            concatenated = [corrected(1, 1:end-3) ...
                            corrected(2, 1:end-3) ...
                            corrected(3, 1:end-3)];
            [~, S] = gfdeconv(concatenated, g_x);
            v = sum(xor(codewords(rand_ind), concatenated));
            nErrBits_withoutHamming = nErrBits_withoutHamming + v;
            if (bi2de(S) == 0)
                if (v > 0)
                    nErrDecode = nErrDecode + 1;
                end
                nSent = nSent + 1;
                break;
            end
        end
    end
    disp([i]);
    Ped(i) = nErrDecode/nTests;
    Pe_bit(i) = nErrBits/(nTests*(3*12));
    Pe_bit_without_Hamming(i) =  nErrBits_withoutHamming/(nTests*(k + r));
    T(i) = k * nSent / (36 * nTests);
end

figure;
semilogy(SNRdB, Ped, 'ko', ...
         SNRdB, Ped_theor_up, 'r.-', ...
         SNRdB, Ped_theor, 'b.-')
legend('Ped', 'Ped theor up', 'Ped theor');
xlabel('E/N_0, dB')

figure;
semilogy(SNRdB, Pe_bit, 'ko', ...
         SNRdB, Pe_bit_without_Hamming, 'r.-', ...
         SNRdB, Pe_bit_theor, 'b.-');
legend('Pe bit', 'Pe bit without Hamming', 'Pe bit theor');
xlabel('E/N_0, dB')

figure;
semilogy(SNRdB, T, 'b.-');
xlabel('E/N_0, dB')