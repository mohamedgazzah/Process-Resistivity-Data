function [] = PlotMRandHallFunc(filename,Ch1,Ch2,Ch3,Hysteresis,Diff)
    DataAll = ProcessResData(filename,Ch1,Ch2,Ch3,Hysteresis,Diff);
    l1 = 0.826; l2 = 0.821; w = 0.411; t = 0.195;
    figure(1)
    for i=1:length(DataAll.loopdata)
        title("Hall")
        plot(DataAll.loopdata{i}.MagneticField,DataAll.loopdata{i}.Resistance1*10^8,"DisplayName",num2str(round(mean(DataAll.loopdata{i}.Temperature))),'LineWidth',1.5)
        legend show
        xlabel("{\it B} (T)")
        ylabel("{\it \rho} (\mu\Omega cm)" )
        xlim([-9 9])
        grid on
        hold on 
    end
    figure(2)
    for i=1:length(DataAll.loopdata)
        title("Hall")
        plot(DataAll.loopdata{i}.MagneticField,DataAll.loopdata{i}.Resistance2*10^8,"DisplayName",num2str(round(mean(DataAll.loopdata{i}.Temperature))),'LineWidth',1.5)
        legend show
        xlabel("{\it B} (T)")
        ylabel("{\it \rho} (\mu\Omega cm)" )
        xlim([-9 9])
        grid on
        hold on 
    end
    figure(3)
    for i=1:length(DataAll.loopdata)
        title("MR")
        plot(DataAll.loopdata{i}.MagneticField,DataAll.loopdata{i}.Resistance3*10^8,"DisplayName",num2str(round(mean(DataAll.loopdata{i}.Temperature))),'LineWidth',1.5)
        legend show
        xlabel("{\it B} (T)")
        ylabel("{\it \rho} (\mu\Omega cm)" )
        xlim([-9 9])
        grid on
        hold on 
    end
    figure(4)
    for i=1:length(DataAll.loopdata)
        title("MR")
        plot(DataAll.loopdata{i}.MagneticField,DataAll.loopdata{i}.DeltaMR3,"DisplayName",num2str(round(mean(DataAll.loopdata{i}.Temperature))),'LineWidth',1.5)
        legend show
        xlabel("{\it B} (T)")
        ylabel("MR (%)" )
        xlim([-9 9])
        grid on
        hold on 
    end