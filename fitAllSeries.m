function none = fitAllSeries(segmentType)
    none = 0;    
    if segmentType == 1
        balancedFit();
        
    elseif segmentType == 2
        dynamicFit();
    end

end

function none = dynamicFit()
    none = 0; %ÿ������Ϊ         1*length
    %��Ϻ������Ϊ     1*��length/segmentSize��*3
    global fitSeries;           %��Ͻ��       cell����
    global disturbedSeries;     %���Һ������   cell����
    global seriesLength;
    global seriesCount;
    
    none = 1;
    fitSeries = [];
    seriesSize =size(disturbedSeries{1});
    seriesLength=seriesSize(2); %���г���
    for i = 1:seriesCount
        currentFitSeries = getMeanDynamicSegment(disturbedSeries{i}(:,1:seriesLength));
        fitSeries{i} = currentFitSeries;
    end
end


function none = balancedFit()
    %ÿ������Ϊ         1*length
    %��Ϻ������Ϊ     1*��length/segmentSize��*3
    global fitSeries;           %��Ͻ��       cell����
    global disturbedSeries;     %���Һ������   cell����
    global seriesLength;
    global seggedLength;
    global seriesCount;
    global segmentSize;
    global segmentIndex;
    global reconstructedSeries;
    global fitseriesMP;
    global reconstructedVariance;
    global reconstructedMean;
    none = 1;
    fitSeries = [];
    seriesSize =size(disturbedSeries{1});
    seriesLength=seriesSize(2); %���г���
    seriesLength = seriesLength-rem(seriesLength,segmentSize);
    seggedLength = seriesLength/segmentSize;
    segmentIndex = 1:1:segmentSize;%[1 2 3 ... segmentSize]
    for i = 1:seriesCount
        currentFitSeries = getBalancedSegment(disturbedSeries{i}(:,1:seriesLength));
        fitSeries{i} = currentFitSeries;
        reconstructedSeries{i}=reconstruct(currentFitSeries);
    end
    [reconstructedMean,reconstructedVariance] = GetMpDestribution(reconstructedSeries);
end

function reconstructed=reconstruct(string)
    global segmentSize;
    global seriesLength;
    reconstructed=[];
    [dimension,length]= size(string);
    tempSegmentIndex = [1:segmentSize];
    for  i = 1:length
        reconstructed = [reconstructed polyval(string(:,i),tempSegmentIndex)];
    end
end


function segments = getMeanDynamicSegment(currentSeries)
    global seriesLength;
    global dimension;
    global errorThreshold;
    errorThreshold = 0.2;
    segments = [];
    tempSegment = [];
    tempSegmentIndex = [];
    lastIndex = 0;
    for j =1:seriesLength
%         temp = mean(currentSeries(:,j));
        temp = max(currentSeries(:,j))+min(currentSeries(:,j));
        tempSegment = [tempSegment,temp];
        tempSegmentIndex = [tempSegmentIndex,j-lastIndex];
        tempFit = polyfit(tempSegmentIndex,tempSegment,dimension);
        errorValue = calcError(tempSegment,tempSegmentIndex,tempFit);
        if errorValue > errorThreshold
            segments = [segments,tempFit'];
            tempSegment=[];
            tempSegmentIndex=[];
            lastIndex=j;
        end
    end
end
%     for time = 1:5
%             
%         for j = 1:seriesLength
%             randonSelection = round(unifrnd(1,5));
%             temp = currentSeries(randonSelection,j);
%             tempSegment = [tempSegment, temp];
%             tempFit = polyfit(tempSegmentIndex,tempSegment,dimension);
%             errorValue = getError(tempSegment,tempSegmentIndex,tempFit);
%             if errorValue > errorthreshold
%                 segments = [segments,tempFit];
%                 funcParaMx = [funcParaMx;tempFit];
%                 tempSegment=[];
%                 tempSegmentIndex=[];
%             end
%         end
%     end
function segments = getBalancedSegment(currentSeries)
    global seriesLength;
    global segmentSize;
        
    global segmentIndex;
    global dimension;
    global method;
    segments=[];
    for j = 1:segmentSize:seriesLength
        subSeries = currentSeries(:,j:j+segmentSize-1);
        meanSeries = mean(subSeries);%����ÿ��ƽ��ֵ
%         meanSeries = max(subSeries)+min(subSeries);
        fit = polyfit(segmentIndex,meanSeries,dimension);%��ֵ���

        if method ~= 1
            %�������ѡ����󣬼������������Ͻ����ȡƽ��
            randomSelection = round(unifrnd(1,5,6,segmentSize));
            for selection = 1:6
                selectedSeries = [];
                tempSegment=[];
                for selectIndex = 1:segmentSize
                    selectedSeries = [selectedSeries subSeries(selectIndex,randomSelection(selection,selectIndex))];
                end
                tempSegment = tempSegment + polyfit(segmentIndex,selectedSeries,dimension);
            end
            fit = fit * 0.5 + tempSegment/12.0;
        end
%         ��ԭ��Ͻ��
        for step = 1:segmentSize
            
        end
        segments = [segments,fit'];%��ϵ��������ת��Ϊ������

    end

end

function errorValue = calcError(segment,index,func)
    SSE = 0;
    avg = mean(segment);
    length = size(index(1,:),2);
    for i=1:length
        fitValue = polyval(func,index(1,i));
        SSE = SSE + (fitValue - segment(i))^2;
    end
    errorValue = sqrt(SSE)/abs(avg);
end

function errorValue = getError(value,index,func)
    global ts;
    global currentTs;
    SSE = 0;
    avg = mean(value,2);
    length = size(index(1,:),2);
    for i=1:length
        fitValue = polyval(func,index(1,i));
        SSE = SSE + (fitValue - ts(currentTs,index(1,i)))^2;
    end
    errorValue = sqrt(SSE)/abs(avg);
end