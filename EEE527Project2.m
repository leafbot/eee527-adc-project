clear all; clc;

Vin = 0.432;
Av = 10E7;
N = 6;
C2 = 1;
C1 = 1; 
vh = 0.25;
vl = -0.25;
vrefp = 1;
vrefm = -1;
vref = vrefp - vrefm;
i = 1;
j = 2;


while i <= N
    
    if i > 1
        Vin = Vx(i-1);
    end
    
    if Vin > vh
        if i < N
            Vx(i) = ((1 + C2/C1) * Vin - C2/C1 * vref)/(1+(1+C2/C1)/Av);
        end
        a(N-i+1) = 1;
        b(N-i+1) = 0;
    else
        if Vin < vh && Vin > vl
            if i < N
                Vx(i) = ((1 + C2/C1) * Vin)/(1+(1+C2/C1)/Av);
            end
            a(N-i+1) = 0;
            b(N-i+1) = 1;
        else
            if Vin < vl
                if i < N
                    Vx(i) = ((1 + C2/C1) * Vin + C2/C1 * vref)/(1+(1+C2/C1)/Av);
                end
                a(N-i+1) = 0;
                b(N-i+1) = 0;
            end
        end
    end
   i = i +1;
end

while j <= N 
    binary(N-j+2) = b(j) + a(j-1);
    
    j = j + 1;
    
    binary(1) = a(N); 
end

binary

xs = linspace(0,N,N-1)
plot (xs,Vx)
Vx 
a
b

       
        