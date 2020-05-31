function  buildFFTmaskSSB(N, Fs, loCut, fftsize)
%function buildFFTmaskSSB(N, Fs, loCut, fftsize)
%                                              
% Scrive su file la maschera per una fast convolution fatta con la FFT.
% La maschera � shiftata rispetto all'asse f=0 verso le frequenze positive
% in modo che il passband inizia a zero Hz e termina a +fpass Hz, in
% modo da poter essere usata per isolare la USB di un segnale. La LSB
% verr� preventivamente specchiata rispetto all'asse f=0
%
% N       = numero di punti del kernel, pari a fftsize - bufsize + 1 
% Fs      = sampling frequency
% loCut   = inizio (lato basso) della banda passante
% hiCut   = termine (lato alto) della banda passante
% fftsize = size della FFT, che determina la estensione con zeri del
%           kernel
%                        Alberto I2PHD  Giugno 2013
%                        rivisto Febbraio 2014, Marzo 2014 e Giugno 2015
%                      ========== usare questo ==========
Wpass = 1;     % Passband Weight
Wstop = 10;    % Stopband Weight
dens  = 20;    % Density Factor

NumFilt = 2;
BW=[2200 2800]; ...

realf = zeros(NumFilt, 1024);
imagf = zeros(NumFilt, 1024);

fnameR = sprintf('E:\\Ham\\STM32F429ZIT6\\Projects\\Keil_Contest\\FFTmaskSSB_R.txt');
fnameI = sprintf('E:\\Ham\\STM32F429ZIT6\\Projects\\Keil_Contest\\FFTmaskSSB_I.txt');
fidR=fopen(fnameR, 'w');
fidI=fopen(fnameI, 'w');

for iter=1:1:NumFilt
    hiCut=BW(iter);
msg = sprintf('Costruzione del filtro SSB da %d Hz', hiCut);
disp(msg);

transwidth = 200;            % larghezza della zona di transizione
flp = (hiCut - loCut)/2 - 40;     % larghezza del prototipo low pass
flpst = flp + transwidth;    % inizio stop band per il LP prototipo
fshift = flp + loCut + 40;   % ammontare dello shift verso destra

% Calcola i coefficienti usando la funzione FIRPM
b  = firpm(N-1, [0 flp flpst Fs/2]/(Fs/2), [1 1 0 0], [Wpass Wstop], ...
           {dens});
%Hd = dfilt.dffir(b);

z=zeros(1, fftsize - N);
bz=[b z];           % espandi il kernel fino al size della FFT
w=2*pi*fshift/Fs;   % calcola l'argomento dell'esponenziale complesso 
                    % per lo shift in frequenza
idx=0:1:fftsize-1;
bs=bz.*exp(1i*w*idx); % shifta verso l'alto la risposta in frequenza
fvtool(bs)
f=fft(bs);
%plot(10*log(abs(f)))   ;pause;
realf(iter,:)=real(f); imagf(iter,:)=imag(f);
end    

for iter=1:1:NumFilt
if(fidR > 0)
    if(iter == 1)
      fprintf(fidR, 'const float FFTmaskSSB_R[%d][1024] =\n{\n', NumFilt);
    end
    
    fprintf(fidR, '// SSB_R PB width %d Hz\n{\n', BW(iter));
    fprintf(fidR, '%13.10ff, %13.10ff, %13.10ff, %13.10ff, %13.10ff, %13.10ff, %13.10ff, %13.10ff, \n', realf(iter,:));
    if(iter == NumFilt) fprintf(fidR, '}\n};\n');
    else  fprintf(fidR, '},\n');
    end    
else
    disp('Errore aprendo il file E:\\Ham\\STM32F429ZIT6\\Projects\\Keil_Contest\\FFTmaskSSB_R.txt');
end
end
fclose(fidR);

for iter=1:1:NumFilt
if(fidI > 0)
    if(iter == 1)
      fprintf(fidI, 'const float FFTmaskSSB_I[%d][1024] =\n{\n', NumFilt);
    end  
    fprintf(fidI, '// SSB_I PB width %d Hz\n{\n', BW(iter));
    fprintf(fidI, '%13.10ff, %13.10ff, %13.10ff, %13.10ff, %13.10ff, %13.10ff, %13.10ff, %13.10ff, \n', imagf(iter,:));
    if(iter == NumFilt) fprintf(fidI, '}\n};\n');
    else  fprintf(fidI, '},\n');
    end    
else
    disp('Errore aprendo il file E:\\Ham\\STM32F429ZIT6\\Projects\\Keil_Contest\\FFTmaskSSB_I.txt');
end    
end
fclose(fidI);
end
% [EOF]
