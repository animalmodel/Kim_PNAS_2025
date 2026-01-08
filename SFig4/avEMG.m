%%%make averaged EMG
% SR, stampSR = kHz
% sttime,endtime = ms
% 

function [edata] = avEMG(data,SR,stamp,stampSR,sttime,endtime)
if SR == stampSR
else
    stamp = round(stamp.*SR./stampSR);
end 
emg_row = [];
for trial = 1:max(size(stamp))
    st=stamp(trial)+sttime*SR; 
    et=stamp(trial)+endtime*SR; 
    if st>1 && et<length(data)
        emg_row = [emg_row;data(st:et)];
    end
end
edata = mean(emg_row,1);
