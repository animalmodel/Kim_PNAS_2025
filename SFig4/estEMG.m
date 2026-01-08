%%%make reconstructed EMG
function [edata] = estEMG(temp,SR,hdat,binlength)
comh = [hdat;zeros(binlength*SR-1,length(hdat))]; %% hdat=1*N?
comh = reshape(comh,1,binlength*SR*length(hdat));  %for comvolution

edata = conv(comh,temp);
