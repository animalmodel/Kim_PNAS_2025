load('Spinalgain_EMG_relationship')

figure
   for perf = 1:4
      dat = emg_para{1,perf}(:);
      [N,edges] = histcounts(dat,10);
      subplot(2,2,perf), bar(edges(2:end),N./length(dat))
      hold on, plot([0.05 1.05],[0.1 0.1],':k')
      title([emg_para{2,perf},' ',emg_para{3,perf}])
      p_row = ones(1,10);
      for b = 10:10
      p=myBinomTest(N(b),length(dat),0.1);
      p_row(b) = p;
      end
      if sum(p_row<0.05) > 0
          ind = find(p_row<0.05);
          for i = ind
              if p_row(i)<0.01
                    text(ind*0.1,0.15,'**')
              else
                    text(ind*0.1,0.15,'*')
              end
          end
      end

   end
