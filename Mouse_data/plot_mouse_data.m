Mou = [1,2,3,4,5];

addpath('Results\')

for mm = 1:5
    mouse = num2str(Mou(mm));    
    data{mm} = importdata(['../Mouse_data/pH',mouse,'_1024.dat']);
end

figure(1)
plot(cell2mat(data),'LineWidth',2)
xlim([0,1024])
title('Hypoxic mouse pressure data')
    name = ['MouseDataPressure.png'];
    set(gcf,'PaperPositionMode','auto');
    print(gcf,name,'-dpng','-r0')       

