% Porous Substrate Electroporation Circuit Model
% By Justin Brooks
% University of Nebraska-Lincoln
% Department of Mechanical and Materials Engineering
% Used in the publication "An Equivalent Circuit Model for Localized
% Electroporation on Track Etched Membranes" by Brooks et al.

clear
clc
close all

%% USER DEFINED INPUTS

%WAVEFORM PARAMETERS
square = true; %true or false, use square wave parameters? if false, exponential decay parameters will be used
HV = 25; %V, high voltage
LV = 0; %V, lower voltage in bilevel pulse, should be set to 0 for unilevel pulses
HVdur = 10; %ms, high voltage duration
LVdur = 0; %ms, low voltage duration
freq = 20; %Hz, frequency of pulses
pulses = 1; %pulses per pulse train
trains = 1; %number of pulse trains
trainInt = 0; %ms, interval between pulse trains
tConstant = 0.01; %exponential decay pulse time constant
Vpts = 100; %number of voltages exponential decay pulse is evaluated at

%FOURIER TRANSFORM ACCURACY
tPts = 500; %number of time points
fterms = 100; %fourier terms

%WELL DIMENSIONS
dw = 5*10^-3; %m, well diameter

%ELECTROLYTE
electrolyte = 'PBS'; %options are "PBS" and "DI" (DI water)

%FARADAIC IMPEDANCE
%model parameters can be found under "FARADAIC IMPEDANCE SELECTION"
ZfpFrac = 0.5; %0-1, fraction of faradaic impedance on the anode (positive) side
%(determined by electrode materials, geometry)

%SUBSTRATE
%model parameters can be found under "SUBSTRATE SELECTION"
substrate = 1200; %nm, options are "20" (DI only), "50" (PBS only), "200", "1200", and "2000"

%CELL MONOLAYER AT 100% CONFLUENCY
cells = true; %true or false, include cell monolayer in model?
Rc = 90;
Cc = 0;
TMPf = 0.8; %TMP as a fraction of the total cell monolayer voltage

%PARAMETRIC STUDIES
parametric = false; %true or false, perform parametric study? (used to create Fig. 3D-F)
parameter = 'voltage'; %options are "voltage", "width", or "frequency", parameter for study
paraRange = [1,2.5,5,10,25,50,100,250]; %parametric range
%recommended values for parametric ranges:
%voltage: 1,2.5,5,10,25,50,100,250
%width: 0.5,1,2.5,5,10,25,50,100
%frequency: 0.25,0.5,1,2.5,5,10,25,50

%% CONVERSION TO SECONDS AND WAVEFORM CALCULATIONS

HVdur = HVdur/1000; %s, convert from ms to s
LVdur = LVdur/1000; %s, convert from ms to s
trainInt = trainInt/1000; %s, convert from ms to s

tInt = (1/freq)/tPts; %s, time interval
duty = HVdur/(1/freq); %duty cycle of square unilevel wave

%% FARADAIC IMPEDANCE SELECTION

syms Vvar

if strcmp(electrolyte,'PBS')
    R1f = 89.844;
    R2f = 35.441*Vvar;
    R3f = 2532.7*Vvar^-1.1760;
    Cf =  0.0000020365*Vvar^0.55943;
elseif strcmp(electrolyte,'DI')
    R1f = 0;
    R2f = 73038*Vvar^-0.14135;
    R3f = 28728*exp(-0.14319*Vvar);
    Cf = 0.00000011491*log(Vvar)+0.00000041601;
end

%% SUBSTRATE SELECTION

subList = [20, 50, 200, 1200, 2000]; %list of substrates
subIndex = find(subList == substrate);

LsList = [25, 25, 25, 24, 23] * 10^-6; %m, channel length
dsList = [20, 50, 200, 1200, 2000] * 10^-9; %m, channel diameter
psList = [6*10^8, 6*10^8, 5*10^8, 2.2*10^7, 3*10^6] * 100^2; %channels/m^2, channel density

if strcmp(electrolyte,'PBS')
    R1sList = [0,1142.4, 16.055, 54.655, 83.890];
    R2sList = [0,11130, 61.580, 133.15, 1927.7];
    CsList = [0,0.000014852, 0.00014074, 0.000095155, 0.000088297];
elseif strcmp(electrolyte,'DI')
    R1sList = [131174.6784, 0, 38222.62, 115653.7, 39153.27];
    R2sList = [279825.511866644, 0, 0, 0, 0];
    CsList = [0.000000000924835823327622, 0, 0, 0, 0];
end

Ls = LsList(subIndex);
ds = dsList(subIndex);
ps = psList(subIndex);

R1s = R1sList(subIndex);
R2s = R2sList(subIndex);
Cs = CsList(subIndex);

%% PARAMETRIC LOOP

VeqAvg = zeros(1,length(paraRange)); %preallocate
VfpAvg = zeros(1,length(paraRange)); %preallocate
VsAvg = zeros(1,length(paraRange)); %preallocate
VcAvg = zeros(1,length(paraRange)); %preallocate
VfnAvg = zeros(1,length(paraRange)); %preallocate

if parametric == false
    paraRange = 0;
end

for paraInd = 1:length(paraRange)
    if parametric == true
        if strcmp(parameter,'voltage')
            HV = paraRange(paraInd);
            tInt = (1/freq)/tPts; %time interval, s
        elseif strcmp(parameter,'width')
            HVdur = paraRange(paraInd);
            freq = duty/HVdur; %constant duty cycle
            tInt = (1/freq)/tPts; %time interval, s
        elseif strcmp(parameter,'frequency')
            freq = paraRange(paraInd);
            tInt = 0.1/1000; %time interval, s
        end
    end
    
    %% GENERATE WAVEFORM
    
    if square == true
        HVtpulse = zeros(1,pulses*trains); %pulse time points
        for i = 1:trains %train number
            for j = 1:pulses %pulse number
                if i == 1 && j == 1 %1st train and 1st pulse
                    HVtpulse(pulses*(i-1)+j) = HVdur/2; %HV pulse
                elseif i > 1 && j == 1 %additional trains, 1st pulse per train
                    HVtpulse(pulses*(i-1)+j) = HVtpulse(pulses*(i-1)+j-1) - HVdur/2 + 1/freq + trainInt + HVdur/2;
                else %additional pulses
                    HVtpulse(pulses*(i-1)+j) = HVtpulse(pulses*(i-1)+j-1) - HVdur/2 + 1/freq + HVdur/2;
                end
            end
        end
        LVtpulse = HVtpulse + HVdur/2 + LVdur/2;
        
        tMax = max(HVtpulse) - HVdur/2 + 1/freq;
        t = 0:tInt:tMax;
        V = HV*pulstran(t,HVtpulse,@rectpuls,HVdur)+LV*pulstran(t,LVtpulse,@rectpuls,LVdur);
        
    else %if using exponential
        HVtpulse = zeros(1,pulses*trains); %pulse time points
        for i = 1:trains
            for j = 1:pulses
                if i == 1 && j == 1 %1st train and 1st pulse
                    HVtpulse(pulses*(i-1)+j) = 0; %HV pulse
                elseif i > 1 && j == 1 %additional trains, 1st pulse per train
                    HVtpulse(pulses*(i-1)+j) = HVtpulse(pulses*(i-1)+j-1) + 1/freq + trainInt;
                else %additional pulses
                    HVtpulse(pulses*(i-1)+j) = HVtpulse(pulses*(i-1)+j-1) + 1/freq;
                end
            end
        end
        
        tMax = max(HVtpulse) + 1/freq;
        t = 0:tInt:tMax;
        
        n = 1;
        V = zeros(1,length(t));
        for i = 1:length(t)
            if t(i) >= HVtpulse(n)
                tn = HVtpulse(n);
                if n < length(HVtpulse)
                    n = n+1;
                end
            end
            V(i) = HV*exp(-(t(i)-tn)/tConstant);
        end
    end
    
    %% CALCULATE IMPEDANCES
    
    %determine the voltage values to evaluate at
    if square == true
        if LV == 0 %unilevel square pulses
            VVal = [0, HV];
        else %bilevel square pulses
            VVal = [0, LV, HV];
        end
    else %exponential pulses
        offset = (HV/Vpts)/2;
        VVal = linspace(0+offset,HV-offset,Vpts);
    end
    
    %calculate circuit impedances that are voltage dependent at each voltage
    R2fVal = zeros(1,length(VVal)); %preallocate
    R3fVal = zeros(1,length(VVal)); %preallocate
    CfVal = zeros(1,length(VVal)); %preallocate
    for i = 1:length(VVal)
        if VVal(i) == 0 && strcmp(electrolyte,'PBS')
            R3fVal(i) = 62796; %calculated R3 at V = 40 mV, value goes to inf at 0
        elseif VVal(i) == 0 && strcmp(electrolyte,'DI')
            R2fVal(i) = 92353.4262; %calculated R2 at V = 50 mV, value goes to inf at 0
        else
            R2fVal(i) = double(subs(R2f,Vvar,VVal(i)));
            R3fVal(i) = double(subs(R3f,Vvar,VVal(i)));
        end
        CfVal(i) = double(subs(Cf,Vvar,VVal(i)));
    end
    
    syms s
    
    %calculate faradaic impedance equations in the frequency domain
    Zf = cell(length(VVal),1); %preallocate
    for i = 1:length(VVal)
        Zf{i} = (CfVal(i)*s*R3fVal(i)*(R1f+R2fVal(i))+R3fVal(i))/(CfVal(i)*s*(R1f+R2fVal(i)+R3fVal(i))+1);
    end
    
    if R2s == 0 || Cs == 0
        Zs = R1s; %substrate impedance
    else
        Zs = R2s/(Cs*s*R2s+1)+R1s; %substrate impedance
    end
    
    if cells == true
        if Cc ~= 0
            Zc = 1/(Cc*s) + Rc; %cell impedance
        else
            Zc = Rc;
        end
    else
        Zc = 0;
    end
    
    %% CALCULATE FOURIER TRANSFORM OF WAVEFORM
    fprintf('Parameter %1.0f/%1.0f\n',paraInd,length(paraRange))
    disp('Calculating Fourier transform of applied waveform.')
    
    [a,b,yfit] = Fseries(t,V,fterms);
    
    syms tvar
    a0 = a(1);
    a = a(2:length(a));
    VfL = cell(length(VVal),1); %preallocate
    VsL = cell(length(VVal),1); %preallocate
    Vf = cell(length(VVal),1); %preallocate
    Vs = cell(length(VVal),1); %preallocate
    for i = 0:length(b)
        if i == 0
            Veq = a0/2;
            VeqL = a0/(2*s);
            for j = 1:length(VVal)
                VfL{j} = (Zf{j}*VeqL)/(Zf{j}+Zs+Zc);
                VsL{j} = (Zs*VeqL)/(Zf{j}+Zs+Zc);
                Vf{j} = ilaplace(VfL{j}, s, tvar);
                Vs{j} = ilaplace(VsL{j}, s, tvar);
            end
        else
            Veq = Veq + a(i)*cos((i)*(tvar-tMax/2)*2*pi/tMax) + b(i)*sin((i)*(tvar-tMax/2)*2*pi/tMax);
            VeqL = laplace(a(i)*cos((i)*(tvar-tMax/2)*2*pi/tMax) + b(i)*sin((i)*(tvar-tMax/2)*2*pi/tMax));
            for j = 1:length(VVal)
                clc
                fprintf('Parameter %1.0f/%1.0f\n',paraInd,length(paraRange))
                disp('Calculating Laplace transform of Fourier series.')
                fprintf('Term %1.0f/%1.0f\nVoltage %1.0f/%1.0f\n',...
                    i,length(b),j,length(VVal))
                VfL{j} = (Zf{j}*VeqL)/(Zf{j}+Zs+Zc);
                VsL{j} = (Zs*VeqL)/(Zf{j}+Zs+Zc);
                Vf{j} = Vf{j} + ilaplace(VfL{j}, s, tvar);
                Vs{j} = Vs{j} + ilaplace(VsL{j}, s, tvar);
            end
        end
    end
    
    %% CALCULATE VOLTAGES AT EACH TIME POINT
    
    VeqMat = zeros(1,length(t)); %preallocate
    VfpMat = zeros(1,length(t)); %preallocate
    VsMat = zeros(1,length(t)); %preallocate
    VcMat = zeros(1,length(t)); %preallocate
    VfnMat = zeros(1,length(t)); %preallocate
    for i = 1:length(t)
        clc
        fprintf('Parameter %1.0f/%1.0f\n',paraInd,length(paraRange))
        disp('Calculating circuit voltages in time domain.')
        fprintf('t = %2.1f/%2.1f ms\n',t(i)*1000,tMax*1000)
        
        [minValue,VIndex] = min(abs(V(i)-VVal));
        Vft = double(subs(Vf{VIndex},[tvar,Vvar],[t(i),V(i)])); %voltage magnitude across Zf
        Vst = double(subs(Vs{VIndex},[tvar,Vvar],[t(i),V(i)])); %voltage magnitude across Zm
        
        Vfp = Vft*ZfpFrac; %voltage magnitude across Zfp
        Vfn = Vft*(1-ZfpFrac); %voltage magnitude across Zfn
        
        if cells == true
            Vc = V(i) - Vft - Vst; %voltage magnitude across Zc
        else
            Vc = 0;
        end
        
        V1 = V(i); %voltage at positive node of Zfp
        V2 = V1 - Vfp; %voltage at positive node of Zm
        V3 = V2 - Vst; %voltage at positive node of Zc
        V4 = V3 - Vc; %voltage at positive node of Zfn
        V5 = 0; %voltage at negative node of Zfn (ground)
        
        VeqMat(i) = subs(Veq,tvar,t(i)); %applied voltage at all times
        VfpMat(i) = V1 - V2; %anodic faradaic voltage at all times
        VsMat(i) = V2 - V3; %substrate voltage at all times
        VcMat(i) = V3 - V4; %cell monolayer voltage at all times
        VfnMat(i) = V4 - V5; %cathodic faradaic voltage at all times
    end
    
    VeqAvg(paraInd) = mean(VeqMat);
    VfpAvg(paraInd) = mean(VfpMat);
    VsAvg(paraInd) = mean(VsMat);
    VcAvg(paraInd) = mean(VcMat);
    VfnAvg(paraInd) = mean(VfnMat);
end

%% WAVEFORM VOLTAGE AND ELECTRIC FIELD PLOT

if parametric == false
    figure
    [ax, h] = plot2axes(t*1000,VeqMat,'yscale',1/(Ls*100),'xlim',...
        [min(t)*1000,max(t)*1000],'ylim',[0-0.25*max(V),1.25*max(V)]);
    hold on
    if ZfpFrac == 0.5
        plot(t*1000,VfpMat*2)
    else
        plot(t*1000,VfpMat*2)
    end
    hold on
    plot(t*1000,VsMat)
    if cells == true
        hold on
        plot2axes(t*1000,VcMat,'yscale',1/(Ls*100))
        plot2axes(t*1000,VcMat*TMPf,'yscale',1/(Ls*100))
    end
    if ZfpFrac ~= 0.5
        hold on
        plot(t*1000,VfnMat)
    end
    
    %LABELS
    title('Voltage and Electric Field vs. Time')
    xlabel('Time (ms)')
    yLabel1 = ylabel(ax(1),'Voltage (V)');
    yLabel2 = ylabel(ax(2),'Substrate Channel Electric Field (V/cm)','Rotation',270);
    %yLabel2.Position(1) = yLabel2.Position(1) + 10; %move yaxis label to the right
    
    %LEGEND
    leg = {'Applied Waveform','Anode Faradaic','Substrate',...
        'Cell Monolayer','TMP','Cathode Faradaic'};
    if ZfpFrac == 0.5
        idx = strcmp(leg,'Cathode Faradaic');
        leg(idx) = [];
        idx = strcmp(leg,'Anode Faradaic');
        leg(idx) = {'Faradaic and Bulk'};
    end
    if cells == false
        idx = strcmp(leg,'Cell Monolayer Voltage');
        leg(idx) = [];
        idx = strcmp(leg,'TMP');
        leg(idx) = [];
    end
    legend(leg)
end

%% PARAMETRIC PLOT

if parametric == true
    clear title
    clear xlabel
    clear ylabel
    figure
    
    if ZfpFrac == 0.5
        semilogx(paraRange,VfpAvg*2./VeqAvg*100) %plot faradaic and bulk
    else
        semilogx(paraRange,VfpAvg./VeqAvg*100) %plut anodic faradaic
    end
    
    hold on
    semilogx(paraRange,VsAvg./VeqAvg*100) %plot substrate
    
    if cells == true
        hold on
        semilogx(paraRange,VcAvg./VeqAvg*100) %plot cell monolayer
        semilogx(paraRange,VcAvg*TMPf./VeqAvg*100) %plot TMP
    end
    
    if ZfpFrac ~= 0.5
        hold on
        semilogx(paraRange,VfnAvg./VeqAvg*100) %plot cathodic faradaic
    end
    
    %LABELS
    title('% Applied Voltage vs. Parameter of Interest');
    xlabel('Parameter of Interest');
    ylabel('% Applied Voltage');
    axis([paraRange(1),paraRange(length(paraRange)),0,100]);
    
    %LEGEND
    leg = {'Anode Faradaic','Substrate',...
        'Cell Monolayer','TMP','Cathode Faradaic'};
    if ZfpFrac == 0.5
        idx = strcmp(leg,'Cathode Faradaic');
        leg(idx) = [];
        idx = strcmp(leg,'Anode Faradaic');
        leg(idx) = {'Faradaic and Bulk'};
    end
    if cells == false
        idx = strcmp(leg,'Cell Monolayer Voltage');
        leg(idx) = [];
        idx = strcmp(leg,'TMP');
        leg(idx) = [];
    end
    legend(leg)
end