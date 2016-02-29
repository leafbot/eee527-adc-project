clear all; clc;

Av = 10E7;
N = 6;%10;
C2 = 1;%0.98;
C1 = 1; 
vrefp = 1;
vrefm = -1;
vref = vrefp - vrefm;
vh = vref/8;%4;
vl = -vref/8;%4;
sample = 1;%int64(2^N*13.127901);
Vin = -0.6;%linspace(2*vrefm, 2*vrefp, sample);

%plot(Vin);

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
word = cell(1, sample);
while count1 <= sample
    count2 = 1;
    while count2 <= N
        word{count1} = [word{count1} binaryArray{count1}{count2}]; %Combine the binaries (that are one cell per bit) into one single cell.
        count2 = count2 + 1;
    end
    decimalArray(count1) = bin2dec(word{count1}); %Convert binary into decimal.
    count1 = count1 + 1;
end

histr = hist(decimalArray, 2^N);
decimalArray(decimalArray==0) = [];
decimalArray(decimalArray==2^N-1) = [];
ideal = length(decimalArray)/(2^N-2);
dnl = (histr-ideal)/ideal;

figure(1);
plot(dnl);
title('DNL');

i = 1;
while i < length(dnl)+1
    inl(i) = sum(dnl(1:i));
    i=i+1;
end

figure(2);
plot(inl);
title('INL');

% ramp=decimalArray;
% code=[0:1023];
% HPC(1:1024,1)=0;
% for j=1:1024,
%     HPC(j)=length(find(ramp==j-1));
% end
% %plot(code,HPC);
% HPC_ideal=(length(ramp)-HPC(1)-HPC(end))/(length(code)-2); %calculate ideal HPC
% 
% Input=linspace(-1,1,1022);
% 
% DNL(1:1022,1)=0;
% INL(1:1022,1)=0;
% code2=[0:1021];
% for j=1:1022,
%     DNL(j)=(HPC(j)-HPC_ideal)/HPC_ideal;
%     for m=1:j,
%         INL(j)=INL(j)+DNL(m);
%     end
% end
% plot(code2,DNL);


       
        