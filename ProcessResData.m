%% Process Data Assym and Symm
%[Data] = ProcessResData(filename,Ch1,Ch2,Ch3,Hysteresis,Diff)
%HOW TO USE
%filename is name of the file
%Ch1 is a string either Hall or MR for your first channel
%Ch2 is a string either Hall or MR for your second channel
%Ch3 is a string either Hall or MR for your third channel
%Hysteresis is the temperature at which hysteresis stops

%OUTPUT of this function is a data structure with all the data for each
%temperature labeled T*Temp*K containing all resistivities, the magnetic
%field and raw temperature file for debugging
%loopdata is a structure within Data that contains all the temperatures
%data to be able to loop through more easily 

%IMPORTANT!!!!:IF MEASURMENT WAS RUN BELOW 2.5K BUT NOT AT 1.8K, TEMP
%DISPLAYED WILL BE 1.8K

function [Data] = ProcessResData(filename,Ch1,Ch2,Ch3,Hysteresis,Diff)
%% Load in Data
[Temp, Mag, Res1, Res2, Res3] = LoadResData(filename);

%% If Sequence hasn't finished running


%% If magnetic field starts at 0T

%start at 0 going up to 9
if (Mag(2))<20000 && Mag(2)>0
    IndDelMax = min(find(Mag>89500));
    Mag = Mag([IndDelMax+1:end])/10000;% units T
    %% RES1
    Res1 = Res1([IndDelMax+1:end])*100*1000*1000;% units muOhm cm
    %% RES2
    Res2 = Res2([IndDelMax+1:end])*100*1000*1000;% units muOhm cm
    %% RES3
    Res3 = Res3([IndDelMax+1:end])*100*1000*1000;% units muOhm cm
    %% Temp
    Temp = Temp([IndDelMax+1:end]);
end
%start at 0 going down to -9
if (Mag(2))>-20000 && Mag(2)<0
    IndDelMax=min(find(Mag<-89500));
    Mag=Mag([IndDelMax+1:end])/10000;% units T
    %% RES1
    Res1 = Res1([IndDelMax+1:end])*100*1000*1000;% units muOhm cm
    %% RES2
    Res2 = Res2([IndDelMax+1:end])*100*1000*1000;% units muOhm cm
    %% RES3
    Res3 = Res3([IndDelMax+1:end])*100*1000*1000;% units muOhm cm
    %% Temp
    Temp = Temp([IndDelMax+1:end]);
end
%% Set Data Structures
TemporaryStruct = struct();
Data = struct();
%% Seperate Temperatures
if max(abs(Mag))>100
    Mag=Mag/10000;
end
MagAll=Mag;Res1all=Res1;Res2all=Res2;Res3all=Res3;%units of micro ohm cm 
difftempindices = find(diff(Temp)>Diff);
j = 1;
for i = 1:length(difftempindices) + 1
    if isempty(difftempindices) == 1
        k = length(Temp);
    else
        if i==length(difftempindices) + 1
            k=length(Temp);
        else
            k = difftempindices(i);
        end
    end
    T = round(mean(Temp(j:k)));
    if T == 2
        T = 1.8;
        Tstring = "T1p8K";
    else
        Tstring = "T"+T+"K";
    end
    Res1=Res1all(j:k);Res2=Res2all(j:k);Res3=Res3all(j:k);Mag=MagAll(j:k);
    %% IF in HYSTERESIS temp range Find 4 Quadrants
    if T<=Hysteresis
        
        zeroindex = find(Mag>=-0.001 & Mag<=0.001);
        negnineindex = find(Mag<-8.95);
        posnineindex = find(Mag>8.95);

        %start at 9T
        if Mag(1)>0
            %% RES1
            firstquad1 = [Res1(1:zeroindex(1))];
            secondquad1 = [Res1(zeroindex(1)+1:negnineindex(1))];
            thirdquad1 = [Res1(negnineindex(2):zeroindex(2))];
            fourthquad1 = [Res1(zeroindex(2)+1:end)];

            %% RES2
            firstquad2 = [Res2(1:zeroindex(1))];
            secondquad2 = [Res2(zeroindex(1)+1:negnineindex(1))];
            thirdquad2 = [Res2(negnineindex(2):zeroindex(2))];
            fourthquad2 = [Res2(zeroindex(2)+1:end)];

            %% RES3
            firstquad3 = [Res3(1:zeroindex(1))];
            secondquad3 = [Res3(zeroindex(1)+1:negnineindex(1))];
            thirdquad3 = [Res3(negnineindex(2):zeroindex(2))];
            fourthquad3 = [Res3(zeroindex(2)+1:end)];
        end
        %start at -9T
        if Mag(1)<0
            %% RES1
            firstquad1 = [Res1(1:zeroindex(1))];
            secondquad1 = [Res1(zeroindex(1)+1:posnineindex(1))];
            thirdquad1 = [Res1(posnineindex(2):zeroindex(2))];
            fourthquad1 = [Res1(zeroindex(2)+1:end)];

            %% RES2
            firstquad2 = [Res2(1:zeroindex(1))];
            secondquad2 = [Res2(zeroindex(1)+1:posnineindex(1))];
            thirdquad2 = [Res2(posnineindex(2):zeroindex(2))];
            fourthquad2 = [Res2(zeroindex(2)+1:end)];

            %% RES3
            firstquad3 = [Res3(1:zeroindex(1))];
            secondquad3 = [Res3(zeroindex(1)+1:posnineindex(1))];
            thirdquad3 = [Res3(posnineindex(2):zeroindex(2))];
            fourthquad3 = [Res3(zeroindex(2)+1:end)];
        end

        %% CH1
        if Ch1 == "Hall"
            %% Assym
            firstthirdquad1 = cat(1,firstquad1,flipud(thirdquad1));
            secondfourthquad1 = cat(1,flipud(secondquad1),fourthquad1);

            assym1and3_1 = (firstthirdquad1 - flipud(firstthirdquad1))/2;
            assym2and4_1 = (secondfourthquad1 - flipud(secondfourthquad1))/2;

            finalfirstquad1 = assym1and3_1(1:floor(length(assym1and3_1)/2));
            finalsecondquad1 = assym2and4_1(1:floor(length(assym2and4_1)/2));
            finalthirdquad1 = assym1and3_1(floor(length(assym1and3_1)/2 + 1):end);
            finalfourthquad1 = assym2and4_1(floor(length(assym2and4_1)/2 + 1):end);
            FinalRes1 = cat(1,finalfirstquad1,flipud(finalsecondquad1),flipud(finalthirdquad1),finalfourthquad1);
        elseif Ch1 == "MR"
            %% Sym
            firstthirdquad1 = cat(1,firstquad1,flipud(thirdquad1));
            secondfourthquad1 = cat(1,flipud(secondquad1),fourthquad1);

            assym1and3_1 = (firstthirdquad1 + flipud(firstthirdquad1))/2;
            assym2and4_1 = (secondfourthquad1 + flipud(secondfourthquad1))/2;

            finalfirstquad1 = assym1and3_1(1:floor(length(assym1and3_1)/2));
            finalsecondquad1 = assym2and4_1(1:floor(length(assym2and4_1)/2));
            finalthirdquad1 = assym1and3_1(floor(length(assym1and3_1)/2 + 1):end);
            finalfourthquad1 = assym2and4_1(floor(length(assym2and4_1)/2 + 1):end);
            FinalRes1 = cat(1,finalfirstquad1,flipud(finalsecondquad1),flipud(finalthirdquad1),finalfourthquad1);
            
            R10 = Interp1NonUnique(Mag,FinalRes1,0);
            DeltaMR1 = (FinalRes1-R10)/R10*100;
        end
        %% CH2
        if Ch2 == "Hall"
            %% Assym
            firstthirdquad2 = cat(1,firstquad2,flipud(thirdquad2));
            secondfourthquad2 = cat(1,flipud(secondquad2),fourthquad2);

            assym1and3_2 = (firstthirdquad2 - flipud(firstthirdquad2))/2;
            assym2and4_2 = (secondfourthquad2 - flipud(secondfourthquad2))/2;

            finalfirstquad2 = assym1and3_2(1:floor(length(assym1and3_2)/2));
            finalsecondquad2 = assym2and4_2(1:floor(length(assym2and4_2)/2));
            finalthirdquad2 = assym1and3_2(floor(length(assym1and3_2)/2 + 1):end);
            finalfourthquad2 = assym2and4_2(floor(length(assym2and4_2)/2 + 1):end);
            FinalRes2 = cat(1,finalfirstquad2,flipud(finalsecondquad2),flipud(finalthirdquad2),finalfourthquad2);
        elseif Ch2 == "MR"
            %% Sym
            firstthirdquad2 = cat(1,firstquad2,flipud(thirdquad2));
            secondfourthquad2 = cat(1,flipud(secondquad2),fourthquad2);

            assym1and3_2 = (firstthirdquad2 + flipud(firstthirdquad2))/2;
            assym2and4_2 = (secondfourthquad2 + flipud(secondfourthquad2))/2;

            finalfirstquad2 = assym1and3_2(1:floor(length(assym1and3_2)/2));
            finalsecondquad2 = assym2and4_2(1:floor(length(assym2and4_2)/2));
            finalthirdquad2 = assym1and3_2(floor(length(assym1and3_2)/2 + 1):end);
            finalfourthquad2 = assym2and4_2(floor(length(assym2and4_2)/2 + 1):end);
            FinalRes2 = cat(1,finalfirstquad2,flipud(finalsecondquad2),flipud(finalthirdquad2),finalfourthquad2);
            
            R20 = Interp1NonUnique(Mag,FinalRes2,0);
            DeltaMR2 = (FinalRes2-R20)/R20*100;
        end
        %% CH3
        if Ch3 == "Hall"
            %% Assym
            firstthirdquad3 = cat(1,firstquad3,flipud(thirdquad3));
            secondfourthquad3 = cat(1,flipud(secondquad3),fourthquad3);

            assym1and3_3 = (firstthirdquad3 - flipud(firstthirdquad3))/2;
            assym2and4_3 = (secondfourthquad3 - flipud(secondfourthquad3))/2;

            finalfirstquad3 = assym1and3_3(1:floor(length(assym1and3_3)/2));
            finalsecondquad3 = assym2and4_3(1:floor(length(assym2and4_3)/2));
            finalthirdquad3 = assym1and3_3(floor(length(assym1and3_3)/2 + 1):end);
            finalfourthquad3 = assym2and4_3(floor(length(assym2and4_3)/2 + 1):end);
            FinalRes3 = cat(1,finalfirstquad3,flipud(finalsecondquad3),flipud(finalthirdquad3),finalfourthquad3);
            
%             R30 = Interp1NonUnique(Mag,FinalRes3,0);
%             DeltaMR3 = (FinalRes3-R30)/R30*100;
        elseif Ch3 == "MR"
            %% Sym
            firstthirdquad3 = cat(1,firstquad3,flipud(thirdquad3));
            secondfourthquad3 = cat(1,flipud(secondquad3),fourthquad3);

            assym1and3_3 = (firstthirdquad3 + flipud(firstthirdquad3))/2;
            assym2and4_3 = (secondfourthquad3 + flipud(secondfourthquad3))/2;

            finalfirstquad3 = assym1and3_3(1:floor(length(assym1and3_3)/2));
            finalsecondquad3 = assym2and4_3(1:floor(length(assym2and4_3)/2));
            finalthirdquad3 = assym1and3_3(floor(length(assym1and3_3)/2 + 1):end);
            finalfourthquad3 = assym2and4_3(floor(length(assym2and4_3)/2 + 1):end);
            FinalRes3 = cat(1,finalfirstquad3,flipud(finalsecondquad3),flipud(finalthirdquad3),finalfourthquad3);
            R30 = Interp1NonUnique(Mag,FinalRes3,0);
            DeltaMR3 = (FinalRes3-R30)/R30*100;
        end
    end
    %% IF NOT in HYSTERESIS temp range One direction Assym/Sym
    if T>Hysteresis
        %% CH1
        if Ch1 == "MR"
            FinalRes1 = (Res1 + flipud(Res1))/2;

            R10 = Interp1NonUnique(Mag,FinalRes1,0);
            DeltaMR1 = (FinalRes1-R10)/R10*100;
        elseif Ch1 == "Hall"
            FinalRes1 = (Res1 - flipud(Res1))/2;
        elseif Ch1 == ""
            FinalRes1 = zeros(length(Res2),1);
        end
        %% CH2
        if Ch2 == "MR"
            FinalRes2 = (Res2 + flipud(Res2))/2;
            R20 = Interp1NonUnique(Mag,FinalRes2,0);
            DeltaMR2 = (FinalRes2-R20)/R20*100;
        elseif Ch2 == "Hall"
            FinalRes2 = (Res2 - flipud(Res2))/2;
        elseif Ch2 == ""
            FinalRes2 = zeros(length(Res3),1);
        end
        %% CH3
        if Ch3 == "MR"
            FinalRes3 = (Res3 + flipud(Res3))/2;
            R30 = Interp1NonUnique(Mag,FinalRes3,0);
            DeltaMR3 = (FinalRes3-R30)/R30*100;
        elseif Ch3 == "Hall"
            FinalRes3 = (Res3 - flipud(Res3))/2;
        elseif Ch3 == ""
            FinalRes3 = zeros(length(Res2),1);
        end
    end
    %% Finding Chai Squared
    STD = [std(FinalRes1),std(FinalRes2),std(FinalRes3)];
    STD = ones(length(Mag),3)*[STD(1) STD(2) STD(3); 0 0 0; 0 0 0];
    %% Putting data in Data Structure
    TemporaryStruct.("Temperature") = Temp(j:k);
    TemporaryStruct.("MagneticField") = Mag;%Mag(j:k);
    TemporaryStruct.("Resistance1") = FinalRes1;%(j:k);
    TemporaryStruct.("Resistance2") = FinalRes2;%(j:k);
    TemporaryStruct.("Resistance3") = FinalRes3;%(j:k);
    TemporaryStruct.("STDV") = STD;
    if Ch1 == "MR"
    TemporaryStruct.("DeltaMR1") = DeltaMR1;
    elseif Ch2 == "MR"
    TemporaryStruct.("DeltaMR2") = DeltaMR2;
    elseif Ch3 == "MR"
    TemporaryStruct.("DeltaMR3") = DeltaMR3;
    end
    
    Data.(Tstring) = TemporaryStruct;
    Data.loopdata{i}=Data.(Tstring);
    j = k + 1;
    
    
    
end




