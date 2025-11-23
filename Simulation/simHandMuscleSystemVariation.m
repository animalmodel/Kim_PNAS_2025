function simHandMuscleSystemVariation()
clear all
close all hidden
global mdl mdlWks sampFreq triggerTime

aGainRate = 0.5;    % aGain is fiexed to be 0.5 times bGain
evalParamDiffSP =1; % Evaluate results of "set point" variation
evalParamDiff1b =1; % Evaluate results of "1b gain" variation
evalParamDiff1a =1; % Evaluate results of "1a gain" variation
outFigs = 1;
writeTableFile = 1;

if(evalParamDiffSP || evalParamDiff1b || evalParamDiff1a)
    sampFreq = 1000;
    [emgFile, dataLabel,~] = xlsread('avEMG_Ext_AgINs_notAgINs', 'AgINs_rec_data', 'A1:X3997');

    triggerTime = 1;
    refTimePre  = (0:1/sampFreq:(triggerTime+emgFile(1,1)))';
    emgData.time = [refTimePre(1:end-1); (emgFile(:,1)+triggerTime)];
    emgRefFile = emgFile(:,2:end);
    emgRefPre  = repmat(emgRefFile(1,:), length(refTimePre)-1, 1);
    emgRef = [emgRefPre;emgRefFile];

    for i=1:length(dataLabel)
        labelMat = strsplit(dataLabel{i});
        dataLabel{i} = labelMat{1};
    end

    normalizeEMG=1;
    if(normalizeEMG)
        emgRefMP = zeros(size(emgRef));
        for i=1:size(emgRef,2)
            %- EMG normalization -%
            %- Remove base-line -%
            rangeEnd  = 0.5*sampFreq+length(refTimePre); % 0.5-1 sec
            rangeInit = 1+length(refTimePre); % Init (-1 sec)
            baseVal = mean(emgRef(rangeInit:rangeEnd, i)); % average value of -1 to 0.5 sec
            emgRefMP(:,i) = emgRef(:,i)-baseVal;
            %--%
            %- Normarize -%
            emgRef(:,i) = emgRefMP(:,i)/max(emgRefMP(:,i));
            %--%
        end
    end
end

mdl = 'modelHandMuscleSystem';
open_system(mdl, 'loadonly');
resultParamCSVName = 'outData/searchResult_Opt.csv';
rParamOrg = readtable(resultParamCSVName);
bGain   = rParamOrg.bGain;
setPoint= rParamOrg.set_point;

if(evalParamDiffSP)
    iList = 1:size(emgRef,2);
    candList = [1.1, 1.05, 1, 0.95, 0.90];
    candName =num2str(candList(:)*100);
    %
    bGainListSP   = zeros(length(iList),length(candList));
    spMListSP     = zeros(length(iList),length(candList));
    durListSP     = zeros(length(iList),length(candList));
    rtListSP      = zeros(length(iList),length(candList));
    ftListSP      = zeros(length(iList),length(candList));
    maxMotListSP  = zeros(length(iList),length(candList));
    maxMotTListSP = zeros(length(iList),length(candList));
    dataLabelSearch = cell(1, length(iList));
    dataLabelSearchLabel = zeros(length(iList),1);
    num=1;
    for i=iList
        mdlWks = get_param(mdl,'ModelWorkspace');
        emgData.signals.values = emgRef(:,i);
        assignin(mdlWks,'emgData',emgData);
        spM0=setPoint(i);
        bGain0 = bGain(i);
        %- Set Parameters -%
        %- Variable Parameter: set point -%
        spMList   = spM0*candList;
        mNList = zeros(length(emgData.signals.values), length(candList));
        spList = zeros(length(emgData.signals.values), length(candList));
        durList = zeros(1,length(candList));
        rtList  = zeros(1,length(candList));
        ftList  = zeros(1,length(candList));
        maxMotList = zeros(1,length(candList));
        maxMotTList= zeros(1,length(candList));
        bGainList = repmat(bGain0, 1, length(candList));
        assignin(mdlWks,'bGain',bGain0)
        aGain0 = aGainRate*bGain0;       
        assignin(mdlWks,'aGain',aGain0); 
        for spI = 1:length(candList)
            disp([num2str(i) ': ' num2str(spMList(spI))]);
            spPattern = getSetPoint(spMList(spI));
            assignin(mdlWks,'spP',   spPattern);
            simOut = sim(mdl); % Variable: motoNeuron is written
            mN = simOut.motoNeuron.Data;
            sp = simOut.spPSim.Data;
            disp(['spI: ' num2str(spMList(spI))]);
            %
            [dur, rt, ft]    = getDuration(simOut.refEMG.time, mN);
            [maxMot,maxMotTI0] = max(mN);
            maxMotT = simOut.refEMG.time(maxMotTI0);
            mNList(1:length(mN),spI)  = mN;
            spList(1:length(sp),spI)  = sp;
            durList(spI) = dur;
            rtList(spI)  = rt;
            ftList(spI)  = ft;
            maxMotList(spI) =maxMot;
            maxMotTList(spI)=maxMotT;
        end
        bGainListSP(num,:)  = bGainList;
        spMListSP(num,:)    = spMList;
        durListSP(num,:)    = durList;
        rtListSP(num,:)     = rtList;
        ftListSP(num,:)     = ftList;
        maxMotListSP(num,:) = maxMotList;
        maxMotTListSP(num,:)= maxMotTList;
        %
        figure;
        legendStr = {'EMG','Optimal set ponit', '+10%','+5%','-5%','-10%'};
        resultFigName = ['outFigs/diffSp_' num2str(i,'%02d') '_' dataLabel{i+1}];
        drawGraphs(simOut.refEMG, mNList, spList, legendStr, outFigs, resultFigName);
        %
        dataLabelSearch{num} = dataLabel{i+1};
        dataLabelSearchLabel(num) = i;
        num=num+1;
    end
    %
    bGainListSPTable  = array2table(bGainListSP,  'VariableNames',cellstr([repmat('bGain_', length(candName),1), candName]));
    spMListSPTable    = array2table(spMListSP,    'VariableNames',cellstr([repmat('SP_', length(candName),1), candName]));
    durListSPTable    = array2table(durListSP,    'VariableNames',cellstr([repmat('Dur_', length(candName),1), candName]));
    rtListSPTable     = array2table(rtListSP,    'VariableNames',cellstr([repmat('rt_', length(candName),1), candName]));
    ftListSPTable     = array2table(ftListSP,    'VariableNames',cellstr([repmat('ft_', length(candName),1), candName]));
    maxMotListSPTable = array2table(maxMotListSP, 'VariableNames',cellstr([repmat('max_', length(candName),1), candName]));
    maxMotTListSPTable= array2table(maxMotTListSP,'VariableNames',cellstr([repmat('maxTime_', length(candName),1), candName]));
    labelTable = array2table(dataLabelSearchLabel, 'VariableNames', {'index'});
    resultSPTable = [labelTable, bGainListSPTable, spMListSPTable, durListSPTable, maxMotListSPTable, maxMotTListSPTable, rtListSPTable, ftListSPTable];
    resultSPTable.Properties.RowNames = dataLabelSearch;
    resultSPTable
    if(writeTableFile)
        resultTableFileName = 'outData/resultSPTable_small.xlsx';
        if(exist(resultTableFileName)>0)
            delete(resultTableFileName);
        end
        writetable(resultSPTable, resultTableFileName,'WriteRowNames',true);
    end
end

if(evalParamDiff1b)
    iList = 1:size(emgRef,2);
    candList = [1.1, 1.05 1, 0.95, 0.9];
    candName =num2str(candList(:)*100);
    %
    bGainList1B   = zeros(length(iList),length(candList));
    spMList1B     = zeros(length(iList),length(candList));
    durList1B     = zeros(length(iList),length(candList));
    rtList1B      = zeros(length(iList),length(candList));
    ftList1B      = zeros(length(iList),length(candList));
    maxMotList1B  = zeros(length(iList),length(candList));
    maxMotTList1B = zeros(length(iList),length(candList));
    dataLabelSearch = cell(1, length(iList));
    dataLabelSearchLabel = zeros(length(iList),1);
    num=1;
    for i=iList
        mdlWks = get_param(mdl,'ModelWorkspace');
        emgData.signals.values = emgRef(:,i);
        assignin(mdlWks,'emgData',emgData);

        spM0=setPoint(i);
        bGain0 = bGain(i);
        %- Variable Parameter: 1b gain -%
        bGainList = bGain0*candList;
        mNList = zeros(length(emgData.signals.values), length(candList));
        spList = zeros(length(emgData.signals.values), length(candList));
        durList = zeros(1,length(candList));
        rtList  = zeros(1,length(candList));
        ftList  = zeros(1,length(candList));
        maxMotList = zeros(1,length(candList));
        maxMotTList = zeros(1,length(candList));
        spMList = repmat(spM0, 1,length(candList));
        spPattern0 = getSetPoint(spM0);
        assignin(mdlWks,'spP',   spPattern0);
        n=1;
        for aGI = 1:length(candList)
            bGainI = bGainList(aGI);
            assignin(mdlWks,'bGain',bGainI);
            aGain = aGainRate*bGainI;
            assignin(mdlWks,'aGain',aGain)
            simOut = sim(mdl); % Variable: motoNeuron is written
            mN = simOut.motoNeuron.Data;
            sp = simOut.spPSim.Data;
            disp(['bGainI: ' num2str(bGainI)]);
            %
            [dur, rt, ft]    = getDuration(simOut.refEMG.time, mN);
            [maxMot,maxMotTI0] = max(mN);
            maxMotT = simOut.refEMG.time(maxMotTI0);
            mNList(1:length(mN),aGI)  = mN;
            spList(1:length(sp),aGI)  = sp;
            durList(aGI) = dur;
            rtList(aGI)  = rt;
            ftList(aGI)  = ft;
            maxMotList(aGI) =maxMot;
            maxMotTList(aGI)=maxMotT;
            n=n+1;
        end
        bGainList1B(num,:)  = bGainList;
        spMList1B(num,:)    = spMList;
        durList1B(num,:)    = durList;
        rtList1B(num,:)     = rtList;
        ftList1B(num,:)     = ftList;
        maxMotList1B(num,:) = maxMotList;
        maxMotTList1B(num,:)= maxMotTList;
        %
        figure;
        legendStr = {'EMG','Optimal 1bGain', '+10%','+5%','-5%','-10%'};
        resultFigName = ['outFigs/diff1BGain_' num2str(i,'%02d') '_' dataLabel{i+1}];
        drawGraphs(simOut.refEMG, mNList, spList, legendStr, outFigs, resultFigName);
        %
        dataLabelSearch{num} = dataLabel{i+1};
        dataLabelSearchLabel(num) = i;
        num=num+1;
    end
    %
    bGainList1BTable  = array2table(bGainList1B,  'VariableNames',cellstr([repmat('bGain_', size(candName,1),1), candName]));
    spMList1BTable    = array2table(spMList1B,    'VariableNames',cellstr([repmat('SP_', size(candName,1),1), candName]));
    durList1BTable    = array2table(durList1B,    'VariableNames',cellstr([repmat('Dur_', size(candName,1),1), candName]));
    rtList1BTable     = array2table(rtList1B,    'VariableNames',cellstr([repmat('rt_', size(candName,1),1), candName]));
    ftList1BTable     = array2table(ftList1B,    'VariableNames',cellstr([repmat('ft_', size(candName,1),1), candName]));
    maxMotList1BTable = array2table(maxMotList1B, 'VariableNames',cellstr([repmat('max_', size(candName,1),1), candName]));
    maxMotTList1BTable= array2table(maxMotTList1B,'VariableNames',cellstr([repmat('maxTime_', size(candName,1),1), candName]));
    labelTable = array2table(dataLabelSearchLabel, 'VariableNames', {'index'});
    %
    result1BTable = [labelTable, bGainList1BTable, spMList1BTable, durList1BTable, maxMotList1BTable, maxMotTList1BTable, rtList1BTable, ftList1BTable];
    result1BTable.Properties.RowNames = dataLabelSearch;
    result1BTable
    if(writeTableFile)
        resultTableFileName = 'outData/result1BTable_small.xlsx';
        if(exist(resultTableFileName)>0)
            delete(resultTableFileName);
        end
        writetable(result1BTable, resultTableFileName,'WriteRowNames',true);
    end
end

if(evalParamDiff1a)
    iList = 1:size(emgRef,2);
    candList = [1.1, 1.05 1, 0.95, 0.9];
    candName =num2str(candList(:)*100);
    %
    aGainList1A   = zeros(length(iList),length(candList));
    spMList1A     = zeros(length(iList),length(candList));
    durList1A     = zeros(length(iList),length(candList));
    rtList1A      = zeros(length(iList),length(candList));
    ftList1A      = zeros(length(iList),length(candList));
    maxMotList1A  = zeros(length(iList),length(candList));
    maxMotTList1A = zeros(length(iList),length(candList));
    dataLabelSearch = cell(1, length(iList));
    dataLabelSearchLabel = zeros(length(iList),1);
    num=1;
    for i=iList
        mdlWks = get_param(mdl,'ModelWorkspace');
        emgData.signals.values = emgRef(:,i);
        assignin(mdlWks,'emgData',emgData);

        spM0=setPoint(i);
        bGain0 = bGain(i);
        aGain0 = bGain(i)*aGainRate;
        %- Variable Parameter: 1a gain -%
        aGainList = aGain0*candList;
        mNList = zeros(length(emgData.signals.values), length(candList));
        spList = zeros(length(emgData.signals.values), length(candList));
        durList = zeros(1,length(candList));
        rtList  = zeros(1,length(candList));
        ftList  = zeros(1,length(candList));
        maxMotList = zeros(1,length(candList));
        maxMotTList = zeros(1,length(candList));
        spMList = repmat(spM0, 1,length(candList));
        spPattern0 = getSetPoint(spM0);
        assignin(mdlWks,'spP',   spPattern0);
        assignin(mdlWks,'bGain',bGain0)
        n=1;
        for aGI = 1:length(candList)
            aGainI = aGainList(aGI);
            assignin(mdlWks,'aGain',aGainI);
            simOut = sim(mdl); % Variable: motoNeuron is written
            mN = simOut.motoNeuron.Data;
            sp = simOut.spPSim.Data;
            disp(['aGainI: ' num2str(aGainI)]);
            %
            [dur, rt, ft]    = getDuration(simOut.refEMG.time, mN);
            [maxMot,maxMotTI0] = max(mN);
            maxMotT = simOut.refEMG.time(maxMotTI0);
            mNList(1:length(mN),aGI)  = mN;
            spList(1:length(sp),aGI)  = sp;
            durList(aGI) = dur;
            rtList(aGI)  = rt;
            ftList(aGI)  = ft;
            maxMotList(aGI) =maxMot;
            maxMotTList(aGI)=maxMotT;
            n=n+1;
        end
        aGainList1A(num,:)  = aGainList;
        spMList1A(num,:)    = spMList;
        durList1A(num,:)    = durList;
        rtList1A(num,:)     = rtList;
        ftList1A(num,:)     = ftList;
        maxMotList1A(num,:) = maxMotList;
        maxMotTList1A(num,:)= maxMotTList;
        %
        figure;
        legendStr = {'EMG','Optimal 1aGain', '+10%','+5%','-5%','-10%'};
        resultFigName = ['outFigs/diff1AGain_' num2str(i,'%02d') '_' dataLabel{i+1}];
        drawGraphs(simOut.refEMG, mNList, spList, legendStr, outFigs, resultFigName);
        %
        dataLabelSearch{num} = dataLabel{i+1};
        dataLabelSearchLabel(num) = i;
        num=num+1;
    end
    %
    aGainList1ATable  = array2table(aGainList1A,  'VariableNames',cellstr([repmat('aGain_', size(candName,1),1), candName]));
    spMList1ATable    = array2table(spMList1A,    'VariableNames',cellstr([repmat('SP_', size(candName,1),1), candName]));
    durList1ATable    = array2table(durList1A,    'VariableNames',cellstr([repmat('Dur_', size(candName,1),1), candName]));
    rtList1ATable     = array2table(rtList1A,    'VariableNames',cellstr([repmat('rt_', size(candName,1),1), candName]));
    ftList1ATable     = array2table(ftList1A,    'VariableNames',cellstr([repmat('ft_', size(candName,1),1), candName]));
    maxMotList1ATable = array2table(maxMotList1A, 'VariableNames',cellstr([repmat('max_', size(candName,1),1), candName]));
    maxMotTList1ATable= array2table(maxMotTList1A,'VariableNames',cellstr([repmat('maxTime_', size(candName,1),1), candName]));
    labelTable = array2table(dataLabelSearchLabel, 'VariableNames', {'index'});
    %
    result1ATable = [labelTable, aGainList1ATable, spMList1ATable, durList1ATable, maxMotList1ATable, maxMotTList1ATable, rtList1ATable, ftList1ATable];
    result1ATable.Properties.RowNames = dataLabelSearch;
    result1ATable
    if(writeTableFile)
        resultTableFileName = 'outData/result1ATable_small.xlsx';
        if(exist(resultTableFileName)>0)
            delete(resultTableFileName);
        end
        writetable(result1ATable, resultTableFileName,'WriteRowNames',true);
    end
end
end

%%
function spPattern = getSetPoint(pointM)
global sampFreq triggerTime

time = (0:1/sampFreq:4)';
spDur=0.1;

spData = zeros(size(time));
a=pointM/spDur;
x=0:(1/sampFreq):spDur;
y1=a*x;
y2=pointM-a*x;
%
pulseF1 = triggerTime*sampFreq+1-length(y1);
pulseF2 = triggerTime*sampFreq;
spData((1+pulseF1):(length(y1)+pulseF1))=y1;
spData((1+pulseF2):(length(y2)+pulseF2))=y2;

spPattern.time = time;
spPattern.signals.values =spData;
end


%% Evaluate Activation Duration
function [dur, tInitList, tEndList] = getDuration(tData, mNData)
triggerTime  =1;
sampFreq     =1000;

%- Normarize -%
mNDataN = mNData/max(mNData);
%--%
%- Evaluate Activation Duration -%
threshold = 0.5;
% Start Time of Activation
iZero = triggerTime*sampFreq;
mNDataNPre  = mNDataN(1:iZero);
mNDataNPost = mNDataN(iZero:end);
if(mNDataN(iZero)>threshold)
    actInit = find(mNDataNPre  <threshold, 1, 'last');
else
    actInit = find(mNDataNPost >threshold, 1, 'first')+iZero;
end
% End time of Activation
actEnd  = find(mNDataN>threshold, 1, 'last');
%--%
tInitList = tData(actInit);
tEndList  = tData(actEnd);
dur = tEndList-tInitList;

end

%%
function drawGraphs(refEMG, mNList, spList, legendStr, outFigs, resultFigName, maxMotT, ft)
subplot(2,1,1)
plot(refEMG.time, refEMG.Data, 'k-', 'LineWidth',2);
hold on
if(size(mNList,2)==5)
    plot(refEMG.time, mNList(:,3)/max(mNList(:,3)), 'b-', 'LineWidth',2);
    plot(refEMG.time, mNList(:,[1,2,4,5])/max(mNList(:,3)));
else
    for i=1:size(mNList,2)
        plot(refEMG.time, mNList(:,i));
    end
    plot([maxMotT, maxMotT], ylim, 'k-');
    plot([ft, ft], ylim, 'k--');
end
hold off
legend(legendStr)
subplot(2,1,2)
hold on
if(size(mNList,2)==5)
    plot(refEMG.time, spList(:,3), 'b-', 'LineWidth',2);
    plot(refEMG.time, spList(:,[1,2,4,5]));
else
    for i=1:size(mNList,2)
        plot(refEMG.time, spList(:,i));
    end
end
hold off
legend(legendStr(2:end))
if(outFigs)
    exportgraphics(gcf,[resultFigName '.png'],'Resolution',300);
    exportgraphics(gcf,[resultFigName '.eps']);
end
end

