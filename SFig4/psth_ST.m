%%%make PSTH
function [h,p,spk_row,ind_row] = psth_ST(unitdata,unitSR, stampdata,stampSR,binlength,sttime,endtime)
if unitSR == stampSR
else
    stampdata = round(stampdata.*unitSR./stampSR);
end 
spk_row=[]; ind_row = [];
for trial = 1:max(size(stampdata))
     st=stampdata(trial)+sttime*unitSR; 
     et=stampdata(trial)+endtime*unitSR; 
     stp=[]; etp=[];
     for pp=1:length(unitdata)-1
         
         if unitdata(pp)<st && unitdata(pp+1)>st
               stp=pp+1;
         end
         if unitdata(pp+1)>et && unitdata(pp)<et
               etp=pp;
         end
     end % for pp
     if isempty(stp) || isempty(etp)
            spike=0;
     else
            spike=unitdata(stp:etp)-stampdata(trial); 
     end
     if spike ~= 0
     spk_row=[spk_row;spike'];
     ind_row = [ind_row;ones(size(spike')).*trial];
     end
end % for trial
% ind_row = ind_row(find(spk_row ~= 0));
% spk_row = setdiff(spk_row,0);
p = [binlength/2*unitSR:binlength*unitSR:(endtime-sttime)*unitSR]+sttime*unitSR;
h = hist(spk_row,p);
h = h./trial;





