clear all; clc;

C2Min = 0.8;
C2Max = 1;
C1 = 1;
C2 = C2Min;
CIndex = 1;

while C2 <= C2Max+0.02

    Adb = 140;
    Av = 10.^(Adb/20);
    N = 12;

    vrefp = 1;
    vrefm = -1;
    vref = vrefp - vrefm;
    vh = vref/4;
    vl = -vref/4;
    Fs = 10e6;

    Amp=1.999;
    freq = 1.777*10^6;
    time=100*10^-9;
    sample = int64(2^14);
    scount = 1;

    while scount <= sample
        Vin(scount) = Amp*sin(2*pi*freq*time*(scount-1));
        scount = scount + 1;
    end

    binaryArray = cell(1,sample);
    binary1 = zeros(1,N);
    binary = cell(1,N);
    count = 1;
    while count <= sample 

        i = 1;
        j = 2;

        while i <= N

            if i > 1
                Vin(count) = Vx(i-1);
            end

            if Vin(count) > vh
                if i < N
                    Vx(i) = ((1 + C2/C1) * Vin(count) - C2/C1 * vref)/(1+(1+C2/C1)/Av);
                end
                a(N-i+1) = 1;
                b(N-i+1) = 0;
            else
                if Vin(count) < vh && Vin(count) > vl
                    if i < N
                        Vx(i) = ((1 + C2/C1) * Vin(count))/(1+(1+C2/C1)/Av);
                    end
                    a(N-i+1) = 0;
                    b(N-i+1) = 1;
                else
                    if Vin(count) < vl
                        if i < N
                            Vx(i) = ((1 + C2/C1) * Vin(count) + C2/C1 * vref)/(1+(1+C2/C1)/Av);
                        end
                        a(N-i+1) = 0;
                        b(N-i+1) = 0;
                    end
                end
            end
            i = i +1;
        end

        while j <= N 
            binary1(N-j+2) = b(j) + a(j-1); %Add them up as decimal. Each bit is represented as one number in an array.

            j = j + 1; 
        end

        binary1(1) = a(N);

        j=N;

        while j >= 1
            if binary1(j) > 1;
                binary1(j) = binary1(j)-2;   %Make sure no number is greater than 1
                if j ~= 1
                    binary1(j-1) = binary1(j-1)+1;    %Carry left if > 1
                end
            end
            j = j - 1;
        end

        j = 1;

        while j <= N
            binary{j} = num2str(binary1(j));   %Convert integer into char cell by cell.
            j = j + 1;
        end

        binaryArray{count} = binary; 

        count = count + 1;
    end

    count1 = 1;
    word = cell(1,sample);
    decimalArray = zeros(1,sample);
    while count1 <= sample
        count2 = 1;
        while count2 <= N
            word{count1} = [word{count1} binaryArray{count1}{count2}]; %Combine the binaries (that are one cell per bit) into one single cell.
            count2 = count2 + 1;
        end
        decimalArray(count1) = bin2dec(word{count1}); %Convert binary into decimal.
        count1 = count1 + 1;
    end

    wind(1,:) = blackmanharris(2^14);
    vo = decimalArray.*wind;
    vof = fft(vo, 2^14);
    vofp = abs(vof).^2/(2^14)^2;
    vofps(1) = vofp(1);
    vofps(2:2^14/2) = 2*vofp(2:2^14/2);
    [Max, SigLoc] = max(vofps(6:length(vofps)));
    SigLoc = SigLoc + 5;
    SignalBins = [SigLoc-4:SigLoc+4];
    s = norm(vofps(SignalBins),1);
    noiseBins = 1:length(vofps);
    noiseBins(SignalBins)=[];
    noiseBins(1:5)=[];
    n = norm(vofps(noiseBins),1);
    %freq1 = [0:Fs/2^14:(Fs/2-Fs/2^14)];
    %figure(2);
    %plot(freq1, 10*log10(vofps));
    %xlabel('Frequency (Hz)');
    %ylabel('Output Power Spectrum (dB)');
    %title('FFT of the decimal output signal');

    snr(CIndex) = 10*log10(s/n);
    enob(CIndex) = (snr(CIndex)-1.76)/6.02;
    CMatching(CIndex) = (1-C2/C1);%*1e6;
    
    C2 = C2 + 0.02;
    CIndex = CIndex + 1;
end

figure(1);
plot(CMatching,enob);
title('Capacitor Matching vs ENOB');

% y2=decimalArray;
% xdft1 = fft(y2,N);
% xdft1 = xdft1(1:N/2+1);
% psdx1 = (1/(N*N)).*abs(xdft1).^2;
% psdx1(2:end-1) = 2*psdx1(2:end-1)
% freq = 0:Fs/length(y2):Fs/2;
% figure(1);
% plot(freq,10*log10(psdx1));
% 
% 
% 
% xlabel('frequency');
% ylabel('power(dB)')
% 
% 
% fn=find(psdx1==max(psdx1));
% signalbins=[fn-3:1:fn+3];
% signalpower=norm(psdx1(signalbins),1);
% noisebins=[1:1:length(xdft1)];
% noisebins(fn-3:1:fn+3)=[];
% noisebins(1:1)=[];
% noisepower=norm(psdx1(noisebins),1);
% snqr=10*log10(signalpower/noisepower)

       
        