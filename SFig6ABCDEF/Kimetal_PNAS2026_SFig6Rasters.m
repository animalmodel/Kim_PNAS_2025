binwidth = 10; % ms
range = [2000 3000]; % ms

for data = 1:2
    if data == 1
        load('dataOKT03901')
    else
        
        load('dataOKT04906')
    end
    
    [bb,aa] = butter(3,5/(TQsRate/2),'low');  
    fTQ = filtfilt(bb,aa,TQ);
%     dfTQ = diff(fTQ);
%     fdfTQ = filtfilt(bb,aa,dfTQ);
    
    EMG = [];
    %%%%%subtraction of mean value and rectified
    EMG(:,1) = (EMG_ECU)';  %% wrist extensor (Radial nerve)
    EMG(:,2) = (EMG_FDS)';  %% finger flexsor (Median nerve)
    
    eEMG=double(EMG);
    [B,A] = butter(3,[1 200]/(EMGsRate/2));  %%% 1 means half of sampleRate (5kHz) = 2.5kHz
    smEMG = filtfilt(B,A,eEMG);
    for sss = 1:size(smEMG,2)
        smEMG(:,sss) = smEMG(:,sss)-mean(smEMG(:,sss),1);
    end
    asmEMG = abs(smEMG);

    ref = MovOnEon;
    
    EMGonset = cell(2,1);
    restEMG = [];
    for rr = 1:size(ref,2)
        restEMG = [restEMG;asmEMG(ref(rr)-1.3*EMGsRate:ref(rr)-0.6*EMGsRate,:)];
    end
    baseline = mean(restEMG,1);
    sdEMG = std(restEMG,1);
    
    
    %calc EMG onset
    if ~isempty(ref)
        
        mnum = 1;
        dEMG = asmEMG(:,mnum); onset = [];
        SD = sdEMG(mnum); BL = baseline(mnum);
        
        for mo = 1:size(ref,2)
            tqs = round(ref(1,mo)/(unitsample/TQsRate));
            if length(find(fTQ(tqs:tqs+2*TQsRate)<-500))>300
                
                stime = round(ref(1,mo)/(unitsample/EMGsRate));
                flag = 0;
                for j = [stime:-1:stime-1*EMGsRate]
                    if dEMG(j)<BL+SD*2
                        flag = flag+1;
                    else
                        flag = 0;
                    end
                    if flag == 0.05*EMGsRate
                        onset = [onset, j+0.05*EMGsRate];
                        break
                    end
                    
                end % for j
            end % if
        end % for mo
        EMGonset{1,1} = onset;
        onset = onset./(EMGsRate/unitsample);
    end % for ee
    
    ind = [1:1:size(onset,2)];
    
    ind = randperm(size(onset,2),50);
    
    TQ_row=[];emg_row=[];
    position = 0;
    spk_row=[];
    for trial = ind
        position = position+1;
        sto=onset(1,trial)-range(1)*unitsample/1000;
        eto=onset(1,trial)+range(2)*unitsample/1000;  %%% total 12,500 points for 5kHz
        stp=[]; etp=[];
        
        if EMGsRate~=unitsample
            st = sto/(unitsample/EMGsRate);
            et = eto/(unitsample/EMGsRate);
        else
            st = sto;   et = eto;
        end
        emg_row(:,:,trial) = asmEMG(st:et,:);
        TQ_row(trial,:) = fTQ(st/(unitsample/TQsRate):et/(unitsample/TQsRate));
        
        subplot(6,2,[6+data 8+data]),
        hold on
        for pp=1:size(unitdata,2)-1
            if unitdata(pp)<sto && unitdata(pp+1)>sto
                stp=pp+1;
            end
            if unitdata(pp+1)>eto && unitdata(pp)<eto
                etp=pp;
            end
        end
        if isempty(stp) || isempty(etp)
            spike=0;
        else
            spike=unitdata(stp:etp)-st;
        end
        
        if ~isempty(spike) && size(spike,2)>2
            plot([spike' spike'],[position,position+0.8],'k')
        elseif ~isempty(spike)
            spike = [spike,0,0];
            plot([spike' spike'],[position,position+0.8],'k')
        end
        spk_row=[spk_row;spike'];
        
    end %%%for trial    
    
    set(gca,'XTick',[0:5].*unitsample+1)
    axis([0 5*unitsample+1 0 50])
    set(gca,'XTickLabel',{'-2000','-1000','EMGon','1000','2000','3000'})    
    
    avmus=mean(emg_row(:,1,:),3);
    avmus2=filter(ones(1,10),10,avmus);  %%moving average for illustaration

    subplot(6,2,data), plot(avmus2,'k')
    title('ECU')
    ht = (round(max(avmus2)/100)+1)*100;
    axis([0 5*EMGsRate 0 ht])
    
    avmus=mean(emg_row(:,2,:),3);
    avmus2=filter(ones(1,10),10,avmus);  %%moving average for illustaration
    
    subplot(6,2,2+data), plot(avmus2,'k')
    title('FDS')
    axis([0 5*EMGsRate 0 ht])
    
    avTQ = mean(TQ_row);
    subplot(6,2,4+data), plot(avTQ,'k')
    title('torque')
    axis([0 5*TQsRate -1000 1500])
    
    spk_row = setdiff(spk_row,0);
    HI = hist(spk_row,[binwidth/2*unitsample/1000:binwidth*unitsample/1000: sum(range)*unitsample/1000]);
    HI = HI./50.*1000/binwidth; % Hz
    
    binlength = binwidth/1000;
    movh = movmean(HI,5); % mean of (bin*5)ms
    movh = resample(movh,1,5);
    subplot(6,2,10+data),bar(([-2:0.05:2.95].*1000),movh), hold on
    
end % for i


