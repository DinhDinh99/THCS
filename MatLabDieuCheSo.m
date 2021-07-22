clc
clear all
text = repmat('CASABLANCA',[1 50]);
t=double(text);
N=length(t)
frequency=zeros(1,128);
for k=0:127
    count=0;
    for i=1:N
        if(t(i)==k)
            count=count+1;
        end
    end
    frequency(k+1)=count;
end
sym=find(frequency)-1;
sym_count=frequency(sym+1);
sym_prob=sym_count/N;
[dict,avglen]=huffmandict(sym,sym_prob);
% Day tin sau ma hoa
enco=huffmanenco(t,dict)
P=length(enco)
% Hieu suat ma hoa
compressionratio=P/(N*8)
%---------------------------------------------------------------------
% Khoi dieu che 16-QAM, Ve chom sao
M = 16; % So diem chom sao
k = log2(M); % S? bit trên môt mau 
sps = 1; % S? l??ng m?u trên m?i ký hi?u (h? s? l?y m?u quá m?c)
dataInMatrix = reshape(enco,length(enco)/k,k)
dataSymbolsIn = bi2de(dataInMatrix);
dataMod = qammod (dataSymbolsIn, M); % Mã hóa nh? phân v?i ?? l?ch pha c?a 0 
scatterplot(dataMod);
hold on
%pho tin hieu sau dieu che
signal = dataMod; %tin hieu sau dieu che
n = length(dataMod);
fs = 65;
fnyquyst = fs/2;
figure;
plot(abs(fft(signal)));
grid on
axis tight
%------------------------------------------------------------
% Khoi kenh truyen
EbNo = 10;
SNR = EbNo+10*log10(k)-10*log10(sps);
receivedSignal = awgn(dataMod,SNR,'measured');
sPlotFig = scatterplot(receivedSignal,1,0,'g.');
hold on
scatterplot(dataMod,1,0,'k*',sPlotFig)
%-----------------------------------------------------------
%pho tin hieu sau kenh truyen co nhieu AWGN
signal1 = receivedSignal;
n = length(receivedSignal);
fs = 65;
fnyquyst = fs/2;
figure;
plot(abs(fft(signal1)));
%-----------------------------------------------------
% Khoi giai dieu che
dataSymbolsOut = qamdemod(receivedSignal,M);
dataOutMatrix = de2bi(dataSymbolsOut,k);
dataOut = dataOutMatrix(:); % Return data in column vector
%-----------------------------------------------------
% Khoi giai ma nguon
fprintf('Nguon tin thu duoc:')
output = huffmandeco(enco,dict) %decoding the data
%-----------------------------------------------------
% Tinh BER
dataIn = reshape(enco,length(enco),1);    % chuyen mtran thanh vector
[numErrors,ber] = biterr(dataIn,dataOut);
fprintf('Ty le loi bit ma hoa nhi phan: %5.2e, Dua tren %d loi.\n',ber,numErrors)

