clear

templength = 0.015;  %s
binlength = 0.005;  %s

ename{1} = 'ECR';
ename{2} = 'ECU';

figure
load(['unit.mat']);
unitSR = SampleRate;
unit = Data;

load(['Movement onset-E.mat'])
MO = Data;
stampSR = SampleRate;

[hdat,p] = psth_ST(unit,unitSR, MO,stampSR,binlength,-1,3);
subplot(3,4,1),bar(hdat,'k')
anum = [0 1 2 3 4] .* unitSR./(unitSR.*binlength);
xlim([0 4*unitSR/(unitSR.*binlength)+1])
set(gca,'XTick',anum+1)
set(gca,'XTickLabel',{'-1','MO','+1','+2','+3s'})

axcomb = 0;
for emg = 1:2
    
    load(['STA (unit, ',ename{emg},').mat']);
    
    templSR = SampleRate;
    spiketime = find((XData(1:end-1).*XData(2:end))<0)+1;
    if isempty(spiketime)
        spiketime = find(XData==0);
    end
    tempbase = pseL1.BaseLine;
    templ = YData(spiketime:spiketime+templength*templSR)-tempbase(spiketime:spiketime+templength*templSR);
    
    subplot(3,4,emg*4+1), plot(XData,YData,'k')
    if emg == 1
        title(['STA n=',num2str(nTrials)])
    end
    hold on
    plot([0 0],[min(YData)-1 max(YData)+1],':k')
    plot([templength templength],[min(YData)-1 max(YData)+1],':k')
    plot(XData,tempbase,':k')
    [comdata] = estEMG(templ,templSR,hdat,binlength);
    
    [B,A] = butter(3,10*(2/templSR),'low');
    comdata = filtfilt(B,A,double(comdata));
    
    subplot(3,4,emg*4+2), plot(comdata,'k')
    anum = [0 1 2 3 4] .* templSR;
    xlim([0 4*templSR+1])
    set(gca,'XTick',anum+1)
    set(gca,'XTickLabel',{'-1','MO','+1','+2','+3s'})
    if emg == 1
        title('Reconstructed EMG ECR')
    else
        title(ename{emg})
    end
    
    load([ename{emg},'.mat']);
    
    emgSR = SampleRate;
    emgdata = Data;
    [B,A] = butter(3,[10 200]*(2/emgSR));
    smEMG_E = filtfilt(B,A,double(emgdata));
    smEMG_E = smEMG_E-mean(smEMG_E);
    asmEMG_E = abs(smEMG_E);
    [B2,A2] = butter(3,10*(2/emgSR));
    low10 = filtfilt(B2,A2,double(asmEMG_E));
    
    
    [mEMG] = avEMG(low10,emgSR,MO,stampSR,-1,3);
    pEMG = max(mEMG);
    mEMG2 = mEMG./pEMG.*100;
    subplot(3,4,emg*4+3), plot(mEMG2,'k')
    anum = [0 1 2 3 4] .* emgSR;
    xlim([0 4*emgSR+1])
    set(gca,'XTick',anum+1)
    set(gca,'XTickLabel',{'-1','MO','+1','+2','+3s'})
    
    if emg == 1
        title(['Actual EMG'])
    end
    
    hold on
    
    %%% EMG onset/offset
    base = mean(mEMG(1:0.5*emgSR));
    sdv = std(mEMG(1:0.5*emgSR));
    if sdv<1
        sdv = 3;
    end
    flag = 0; emgON = 0;
    for t = [0.5*emgSR:1.5*emgSR]
        if mEMG(t)>base+2*sdv && flag<0.2*emgSR
            flag = flag+1;
        elseif flag>=0.2*emgSR
            emgON = t-flag;
            break
        else
            flag = 0;
        end
    end % for t
    flag = 0; emgOFF = 0;
    for t = [3*emgSR:-1:1.5*emgSR]
        if mEMG(t)>base+2*sdv && flag<0.1*emgSR
            flag = flag+1;
        elseif flag>=0.1*emgSR
            emgOFF = t+flag;
            break
        else
            flag = 0;
        end
    end % for t
    
    
    %%%comparison between actual and estimated EMGs
    comdata2 = comdata./pEMG*100;
    
    if emgON>0 && emgOFF>0
        plot([emgON emgON],[0 101],':k')
        plot([emgOFF emgOFF],[0 101],':k')
        
        R = corrcoef(mEMG2(emgON:emgOFF),comdata2(emgON:emgOFF));
        R = R(1,2);
        rate = sum(comdata2(emgON:emgOFF))/sum(mEMG2(emgON:emgOFF))*100;
        %%for figure
        Rsh = round(R*1000)/1000;
        rash = round(rate*1000)/1000;
        %%%
        subplot(3,4,emg*4+4), plot(mEMG2([emgON:100:emgOFF]),comdata2([emgON:100:emgOFF]),'k.')
        if emg == 1
            title(['Correlation R=',num2str(Rsh),' ',num2str(rash),'%'])
        else
            title(['R=',num2str(Rsh),' ',num2str(rash),'%'])
        end
        hold on
        lsline
        
    end
        
    
end % for emg





