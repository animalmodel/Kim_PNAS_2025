function simHandMuscleSystem()
clear all
close all hidden
global mdl mdlWks sampFreq triggerTime aGainRate

outDataFolder = 'outData';
outFigFolder  = 'outFigs';
if ~exist(outDataFolder, 'dir')
    mkdir(outDataFolder);
end
if ~exist(outFigFolder, 'dir')
    mkdir(outFigFolder);
end

sampFreq = 1000;
[emgFile, dataLabel,~] = xlsread('avEMG_Ext_AgINs_notAgINs', 'AgINs_rec_data', 'A1:X3997');

aGainRate = 0.5; % aGain is fiexed to be 0.5 times bGain
triggerTime = 1;
refTimePre  = (0:1/sampFreq:(triggerTime+emgFile(1,1)))';
emgData.time = [refTimePre(1:end-1); (emgFile(:,1)+triggerTime)];
figure
subplot(2,1,1)
plot(emgData.time(1:end-1))
subplot(2,1,2)
plot(diff(emgData.time))
ylim([-1/sampFreq, 2/sampFreq]);
emgRefFile = emgFile(:,2:end);
emgRefPre  = repmat(emgRefFile(1,:), length(refTimePre)-1, 1);
emgRef = [emgRefPre;emgRefFile];

for i=1:length(dataLabel)
    labelMat = strsplit(dataLabel{i});
    dataLabel{i} = labelMat{1};
end

dispEMG(emgData.time, emgRef, dataLabel);

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
        emgRef(:,i) = emgRefMP(:,i)/max(emgRefMP(:,i));
        %--%
    end
    dispEMG(emgData.time, emgRef, dataLabel);
end

doSim=1;
if(doSim)
    mdl = 'modelHandMuscleSystem';
    open_system(mdl);
    % for i=1
    for i=1:size(emgRef,2)
        mdlWks = get_param(mdl,'ModelWorkspace');
        emgData.signals.values = emgRef(:,i);
        assignin(mdlWks,'emgData',emgData);

        %- Setting Initial Values -%
        bGain0List = 1;
        spM0List   = 0.1;
        %
        resultParamList = zeros(length(bGain0List)*length(spM0List), 8);
        resultSeries  = cell(1,length(bGain0List)*length(spM0List));
        num = 1;
        for bGain0Ind = 1:length(bGain0List)
            for spM0Ind = 1:length(spM0List)
                bGain0 = bGain0List(bGain0Ind);
                aGain0 = aGainRate*bGain0;
                spM0   = spM0List(spM0Ind);
                assignin(mdlWks,'bGain',bGain0)
                assignin(mdlWks,'aGain',aGain0);
                spPattern0 = getSetPoint(spM0);
                assignin(mdlWks,'spP',   spPattern0);
                %--%
                simOut = sim(mdl); % Variable: motoNeuron is written
                %--%
                showFigInit =0;
                if(showFigInit)
                    figure;
                    subplot(2,1,1)
                    plot(simOut.refEMG.time, [simOut.refEMG.Data, simOut.motoNeuron.Data]);
                    plot(simOut.refEMG.time, [simOut.refEMG.Data, simOut.motoNeuron.Data/max(simOut.motoNeuron.Data)]);
                    subplot(2,1,2)
                    plot(simOut.spPSim.Time, simOut.spPSim.Data);
                end
                rmseMNeuron0 = rms(simOut.motoNeuron.Data/max(simOut.motoNeuron.Data)-simOut.refEMG.Data);
                disp([num2str(i) ': Initial bGain: ' num2str(bGain0) ', RMSE0(Normalized): ' num2str(rmseMNeuron0)]);
                %- Search for Optimal Parameters -%
                doSearch = 1;
                if(doSearch)
                    searchVar0 = [bGain0, spM0];
                    %-- Search Algorithm: fmincon --%
                    maxIte  = 500;
                    maxEval =1000000;
                    options = optimoptions('fmincon','StepTolerance', 1e-14, 'Display','none', 'MaxIterations', maxIte, 'MaxFunctionEvaluations', maxEval, 'OptimalityTolerance', 1e-8, 'ConstraintTolerance', 1e-12);
                    lb = [0,0]; % Limited for only Positive Values
                    ub = [];
                    [searchVarStar, feval, exitflag, output] = fmincon(@mNeuronSimRMSE,searchVar0, [], [], [], [], lb, ub, [], options);
                    exitflag
                    output.stepsize
                    output.constrviolation
                    %--%
                    if(exitflag==1||exitflag==2), disp('Parameter search finished!'); else, disp('Parameter search un-finished'); end
                    disp([num2str(i) ': Optimal bGain:' num2str(searchVarStar(1)) ', set point :' num2str(searchVarStar(2)) ', RMSE(Normalized): ' num2str(feval) ', Iteration number: ' num2str(output.iterations)]);
                    %--%

                    assignin(mdlWks,'bGain',searchVarStar(1))
                    assignin(mdlWks,'aGain',aGainRate*searchVarStar(1))
                    assignin(mdlWks,'spP',   getSetPoint(searchVarStar(2)));
                    simOut = sim(mdl); % Variable: motoNeuron is written
                    %--%

                    % Remove the first added data and subtract one from the time.
                    initI = length(refTimePre);
                    spTime = simOut.spPSim.Time(initI:end,1)-1;
                    motoNeuronData = simOut.motoNeuron.Data(initI:end,:);
                    refEMGTime = simOut.refEMG.time(initI:end,:);
                    refEMGData = simOut.refEMG.Data(initI:end,:);
                    spPSimData = simOut.spPSim.Data(initI:end,:);
                    resultAngData = simOut.resultAng.Data(initI:end,:);
                    %--%

                    saveResultAll =1;
                    if(saveResultAll)
                        resultFilename = [outDataFolder '/searchResult_All.csv'];
                        if(~exist(resultFilename, 'file'))
                            text ={'File Number', 'File Name', 'Evaluation Value', 'EMG Maximum', 'bGain', 'set-point', 'Iteration Number', 'exitflag', 'Initial 1bGain', 'Initial set-point'};
                            fid = fopen(resultFilename, 'w');
                            for tI = 1:(length(text)-1)
                                fprintf(fid, '%s, ', text{tI});
                            end
                            fprintf(fid, '%s\n', text{end});
                            fclose(fid);
                        end
                        fid = fopen(resultFilename, 'A');
                        fprintf(fid, '%s, ', num2str(i));
                        fprintf(fid, '%s, ', dataLabel{i+1});
                        fprintf(fid, '%s, ', num2str(feval,'%.12f'));   % Evaluation Value
                        fprintf(fid, '%s, ', num2str(max(refEMGData))); % EMG Maximum
                        fprintf(fid, '%s, ', num2str(searchVarStar(1),'%.12f')); % bGain
                        fprintf(fid, '%s, ', num2str(searchVarStar(2),'%.12f')); % spP
                        fprintf(fid, '%s, ', num2str(output.iterations)); % Iteration Number
                        fprintf(fid, '%s, ', num2str(exitflag)); % Was the search succeeded?
                        fprintf(fid, '%s, ', num2str(bGain0));   % Initial Value: 1b Gain
                        fprintf(fid, '%s\n', num2str(spM0));     % Initial Value: Set Point
                        fclose(fid);
                    end
                    resultParamList(num,:) = [feval, max(refEMGData), searchVarStar(1), searchVarStar(2), output.iterations, exitflag, bGain0, spM0];
                    resultSeries{num} = [spTime, motoNeuronData, refEMGData, spPSimData, refEMGTime, resultAngData];
                    num=num+1;
                end
            end
        end
        [~, fevalOptInd] = min(resultParamList(:,1));
        fevalOpt = resultParamList(fevalOptInd,:);
        saveResultOpt =1;
        if(saveResultOpt)
            resultFilename = [outDataFolder '/searchResult_Opt.csv'];
            if(~exist(resultFilename, 'file'))
                text ={'File Number', 'File Name', 'Evaluation Value', 'EMG Maximum', 'bGain', 'set-point', 'Iteration Number', 'exitflag', 'Initial 1bGain', 'Initial set-point'};
                fid = fopen(resultFilename, 'w');
                for tI = 1:(length(text)-1)
                    fprintf(fid, '%s, ', text{tI});
                end
                fprintf(fid, '%s\n', text{end});
                fclose(fid);
            end
            fid = fopen(resultFilename, 'A');
            fprintf(fid, '%s, ', num2str(i));
            fprintf(fid, '%s, ', dataLabel{i+1});
            fprintf(fid, '%s, ', num2str(fevalOpt(1),'%.12f')); % Evaluation Value
            fprintf(fid, '%s, ', num2str(fevalOpt(2)));         % EMG Maximum
            fprintf(fid, '%s, ', num2str(fevalOpt(3),'%.12f')); % bGain
            fprintf(fid, '%s, ', num2str(fevalOpt(4),'%.12f')); % spP
            fprintf(fid, '%s, ', num2str(fevalOpt(5))); % Iteration Number
            fprintf(fid, '%s, ', num2str(fevalOpt(6))); % Was the search succeeded?
            fprintf(fid, '%s, ', num2str(fevalOpt(7))); % Initial Value: 1b Gain
            fprintf(fid, '%s\n', num2str(fevalOpt(8))); % Initial Value: Set Point
            fclose(fid);
        end
        spTime         = resultSeries{fevalOptInd}(:,1);
        motoNeuronData = resultSeries{fevalOptInd}(:,2);
        refEMGData     = resultSeries{fevalOptInd}(:,3);
        spPSimData     = resultSeries{fevalOptInd}(:,4);
        refEMGTime     = resultSeries{fevalOptInd}(:,5);
        resultAngData  = resultSeries{fevalOptInd}(:,6);

        figure;
        subplot(3,1,1)
        plot(refEMGTime, [refEMGData, motoNeuronData]);
        legend('EMG','MotoNeuron');
        subplot(3,1,2)
        plot(refEMGTime, resultAngData);
        legend('hand angle');
        subplot(3,1,3)
        plot(spTime, spPSimData);
        legend('Set Point');
        saveFig=1;
        if(saveFig)
            resultFigName = [outFigFolder '/resultFig' num2str(i,'%02d') '_' dataLabel{i+1}];
            exportgraphics(gcf,[resultFigName '.png'],'Resolution',300);
            exportgraphics(gcf,[resultFigName '.eps']);
        end
        figure;
        subplot(3,1,1)
        plot(refEMGTime, [refEMGData/max(refEMGData), motoNeuronData/max(motoNeuronData)]);
        legend('EMG','MotoNeuron (Norm)');
        subplot(3,1,2)
        plot(refEMGTime, resultAngData);
        legend('hand angle');
        subplot(3,1,3)
        plot(spTime, spPSimData);
        legend('Set Point');
        if(saveFig)
            resultFigName = [outFigFolder '/resultFig' num2str(i,'%02d') '_' dataLabel{i+1}];
            exportgraphics(gcf,[resultFigName 'Norm.png'],'Resolution',300);
            exportgraphics(gcf,[resultFigName 'Norm.eps']);
        end
        saveCSV=1;
        if(saveCSV)
            exportMat   = [spTime, motoNeuronData, refEMGData, spPSimData];
            exportTable = array2table(exportMat, 'VariableNames', {'time', 'motoneuron', 'EMGref', 'set point'});
            resultCSVName = [outDataFolder '/' num2str(i,'%02d') '_' dataLabel{i+1} '.csv'];
            writetable(exportTable, resultCSVName,'Delimiter',',','QuoteStrings',true);
        end
    end
end
end

%% Execute simulation from bGain0, and returns the (normalized) RMSE from EMGref
function rmseMNeuron=mNeuronSimRMSE(searchVar)
global mdl mdlWks sampFreq aGainRate

assignin(mdlWks,'bGain',searchVar(1))
aGain = aGainRate*searchVar(1);
assignin(mdlWks,'aGain',aGain);
assignin(mdlWks,'spP', getSetPoint(searchVar(2)));
%--%
simOut=sim(mdl); % Variable: motoNeuron is written

if(simOut.refEMG.time(end)==4 && ~isempty(simOut.motoNeuron.Data)) % Simulation finished
    % results after 1 sec are used for evaluation
    evalMotoNeuron = simOut.motoNeuron.Data(1*sampFreq:end);
    evalRefEMG     = simOut.refEMG.Data(1*sampFreq:end);
    rmseMNeuron = rms(evalMotoNeuron/max(evalMotoNeuron)-evalRefEMG);

    % Evaluation of motion angle %
    resultAngMin = min(simOut.resultAng.Data(1*sampFreq:end)); % minimum angle
    if(resultAngMin>29.5), rmseMNeuron = rmseMNeuron+1; end    % penalty: if joint moved less than 0.5 degree

    % Evaluation of oscillation %
    evalMotoNeuronDot = simOut.motoNeuronDot.Data(1*sampFreq:end);
    weightVel   = 0.0001;
    rmseMNeuron = rmseMNeuron + weightVel*max(evalMotoNeuronDot);
else % Simulation finished with error
    rmseMNeuron = 10;
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

%%
function dispEMG(time, emgRef, dataLabel)
figure;
for i=1:12
    subplot(6,2,i)
    plot(time, emgRef(:,i));
    ylabel(dataLabel{i+1});
end
figure;
for i=13:size(emgRef,2)
    subplot(6,2,i-12)
    plot(time, emgRef(:,i));
    ylabel(dataLabel{i+1});
end
end
