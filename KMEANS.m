function [ri,jaccard,fm,pa,na,aa,meanvalue,runningtume,sc] = KMEANS(cheat,distanceType,eucliMode)
    %参数：：：：：：：：：距离类型：：序列集：：标准差：作弊模式
    fprintf('KMEANS:\n');
    global k;
    global seriesLength;    %序列长度
    global seriesCount;     %序列个数
    global disttype;
    global Labels;
    global seggedLength;
    global disturbedSeries;
    global fitSeries;
    global RawTimeSeries;
    global reconstructedSeries;
    global regionalSeries;
    disttype = distanceType;
    Times = 10;                     %迭代次数上限
    rows = seriesCount;              %序列个数
    length = seggedLength;
    start = 1;                      %类别开始的下标

    if eucliMode == "fit"
        regionalSeries = fitSeries;
    else
        for i = 1:seriesCount
            if eucliMode=="mean"
                regionalSeries{i} = mean(disturbedSeries{i}); 
            elseif eucliMode=="raw"
                regionalSeries{i} =  RawTimeSeries(i,:);
            elseif eucliMode=="fit"
                return;
                regionalSeries{i}= fitSeries{i};
            elseif eucliMode=="reconstructed"
                regionalSeries{i}= reconstructedSeries{i};
            else 
                regionalSeries{i} = disturbedSeries{i}; 
            end
        end
        length = size(regionalSeries{1},2);
    end
    

    for i = 1:seriesCount
        distanceMatrix(i,i)=0;
    end
    %-----开始计时-----%
    StartTime = clock;
    %-----初始簇心-----%
        %计算均值簇心
    if cheat
        for i = 1:k
            meanseries = zeros(1,length);
            count = 0;
            for j = 1:rows
                if Labels(j) == i + start - 1
                    meanseries = meanseries + regionalSeries{j};
                    count = count + 1;
                end
            end
            meanseries = meanseries/count;
            means{i} = {i+start-1 meanseries };
            means{i}{3} = 1;
            means{i}{4} = 0;
        end
    else
        %初始化簇心为随机选择
        for i = 1:k
            selection = round(rows/k*i);%round(rand(0,1)*rows);
            means{i}={i+start-1 regionalSeries{selection} 1 selection};
        end
    end
    
    %-----进入迭代-----%
    %   保存簇心分配
    ClusterAssignment = zeros(1,rows);
    
    for time = 1:Times
      
        fprintf('迭代%f\n',time);
        %保留上次迭代的簇心
        for i = 1:k
            count = means{i}{3};
            LastIteration{i} = {means{i}{2}/count count means{i}{4}};
        end
        
        
        ClusterChange = 0;
        %分配簇心
        for i = 1:rows
            mindistance = 100000000.0;
            mincluster = -1;
            for j=1:k
                distance = DISTANCE(LastIteration{j}{1},regionalSeries{i});
                if distance < mindistance
                   mindistance = distance;
                   mincluster = means{j}{1};
                end
            end
            %如果簇心发生变化
            if ClusterAssignment(i) ~= mincluster
                ClusterChange = ClusterChange +1;
            end
            ClusterAssignment(i) = mincluster;
                           
            means{mincluster}{2} = means{mincluster}{2} + regionalSeries{i};
            means{mincluster}{3} = means{mincluster}{3} + 1;            
            

        end

        %重新计算均值簇心
        for i=1:k
            means{i}{2} = means{i}{2}/means{i}{3};
            means{i}{3} = 1;
            % if distancetype == 3
            %     means{i}{2} = round(means{i}{2});
            % end
             %寻找medoids
             tempmindistance = 99999999;
             tempminindex = -1;
             fit = 0;
             for temp = 1:rows
                if ClusterAssignment(temp) == i
                    fit = fit+1;
                   tempdistance = DISTANCE(means{i}{2},regionalSeries{temp});
                   if tempdistance<tempmindistance
                       tempmindistance = tempdistance;
                       tempminindex = temp;%记录当前距离簇心最近的时间序列
                   end 
                end
             end
             if fit == 0 
               time = Times;%如果遇到空簇 直接退出循环
               break;
             end
 
             if fit ~= 0 
                 means{i}{2} = regionalSeries{tempminindex};
                 means{i}{4} = tempminindex;
             end
        end
%         ClusterAssignment
        %进入迭代后，如果利用分配好的均值簇心，则在第一次循环中退出
        %如果簇心分配没有改变 退出循环
        if cheat || ClusterChange==0
            break;
        end
    end
    
    %结束计时
    EndTime = clock;
    a = 0;
    b = 0;
    c = 0;
    d = 0;

    % 计算距离矩阵,然后计算轮廓系数
    distanceMatrix = zeros(seriesCount,seriesCount);
    distanceMatrix(:,:) = -1;
    for i=1:seriesCount
        for j = 1:seriesCount
            distanceMatrix(i,j) = DISTANCE(regionalSeries{i},regionalSeries{j});
            distanceMatrix(j,i) = distanceMatrix(i,j);
        end
    end
    sc = mean(getSC(distanceMatrix,ClusterAssignment));


    %-----计算结果-----%
    for i = 1:rows
        for j = 1:rows
           if i~=j
               c = true;
               p = true;
                %  P
               if ClusterAssignment(i) == ClusterAssignment(j)
                  p = true;
              else
                  p = false;
              end
                %  C
              if  Labels(i) ==  Labels(j)
                  c = true;
              else
                  c = false;
              end
               
              if c
                 if p
                     a=a+1;
                 else
                     b=b+1;
                 end
              else
                  if p
                     c=c+1;
                 else
                     d=d+1;
                 end
              end
           end
        end
    end
    ri = (a+d)/(a+b+c+d);
    jaccard = a/(a+b+c);
    fm = sqrt(a/(a+b)*a/(a+c));
    pa = a/(a+c);
    na = d/(b+d);
    aa = (pa + na)/2;
    meanvalue = (ri+jaccard+fm)/3;
    runningtume = etime(EndTime,StartTime);

   

end



function distance = DISTANCE(stringa, stringb)
    global disttype;
    if disttype== 1
        distance = EUCLIDEAN(stringa,stringb);
    elseif disttype == 2
        distance = PUREEUCLIDEAN(stringa,stringb);
    elseif disttype == 3
        distance = myDTW(stringa,stringb);
    elseif disttype == 4
        distance = dtw(stringa,stringb);
    elseif disttype == 5
        distance = Pearson(stringa,stringb);
    elseif disttype == 6
        distance = segmentPearson(stringa,stringb);
	elseif disttype == 7
        global reconstructedMean;
        global reconstructedVariance;
        distance = MpDistance(stringa,stringb,reconstructedMean,reconstructedVariance);
    elseif disttype == 8
        distance = FastDTW(stringa,stringb);
    end
end

%% 对每个系数分别计算DTW
function dist = myDTW(stringa,stringb)
    global dimension;
    dist = 0;
    for tempDimension = 1:dimension+1
        dist = dist + dtw(stringa(tempDimension,:),stringb(tempDimension,:));
    end
end

function dist = FastDTW(stringa,stringb)
    global fuzzyThres;
    UB1 = max(stringa);
    UB2 = max(stringb);
    LB1 = min(stringa);
    LB2 = min(stringb);
    fuzzyDist = (UB1+LB1-UB2-LB2)/2;
    if fuzzyDist < 1.5
        dist = fuzzyThres;
    else
        dist = dtw(stringa,stringb);
    end
end

function distance = MpDistance(stringa,stringb,meanMatrix,varianceMatrix)
    global seriesLength;
    distance = 0;
    for i = 1:seriesLength
        distance = distance + abs(normcdf(stringa(i),meanMatrix(1,i),varianceMatrix(1,i))-normcdf(stringb(i),meanMatrix(1,i),varianceMatrix(1,i)));
    end
end


function euclidean =PUREEUCLIDEAN(stringa,stringb)
    
    euclidean = sqrt(sum((stringa-stringb).^2));
end

function euclidean = EUCLIDEAN(stringa,stringb)
    
    euclidean = sqrt(sum(sum((stringa-stringb).^2)));       %??
end


function pearson=Pearson(stringa,stringb)
     [dimension,seriesLength] = size(stringa(:,:));
    pearson = 0;
    for index=1:seriesLength
        eucli = EUCLIDEAN(stringa(:,index),stringb(:,index));
%         temp = corr((stringa(:,index)),(stringb(:,index)),'type','Pearson');
        temp = myPearson(stringa(:,index),stringb(:,index));
        if temp >= -1 &&temp<=1
        else
            temp = 0;
        end
        pearson = pearson + (2-temp)*eucli;
    end
end


function pearson=segmentPearson(stringa,stringb)
    global segmentSize;
     [dimension,seriesLength] = size(stringa(:,:));
    pearson = 0;
    for index = 1:segmentSize:seriesLength+1-segmentSize
        eucli = EUCLIDEAN(stringa(:,index:index+segmentSize-1),stringb(:,index:index+segmentSize-1));
        temp = myPearson((stringa(:,index:index+segmentSize-1))',(stringb(:,index:index+segmentSize-1))');
%         temp = corr((stringa(:,index:index+segmentSize-1))',(stringb(:,index:index+segmentSize-1))','type','Pearson');
        if temp >= -1 &&temp<=1
        else
            temp = 0;
        end
        pearson = pearson + (2-temp)*eucli;
    end
end


% 计算轮廓系数
function sc = getSC(distanceMap,ClusterAssignment)
    global seriesCount;     %序列个数
    nodeInfo = zeros(1,seriesCount);
    for i = 1:seriesCount
        interiorTotal = 0.0;
        interiorCount = 0;
        exteriorTotal = 0.0;
        exteriorCount = 0;
        for j = 1:seriesCount
            if ClusterAssignment(i) == ClusterAssignment(j)
                interiorCount = interiorCount+1;
                interiorTotal = interiorTotal + distanceMap(i,j);
            else
                exteriorCount = exteriorCount+1;
                exteriorTotal = exteriorTotal + distanceMap(i,j);
            end
        end
        b = exteriorTotal/exteriorCount;
        a = interiorTotal/interiorCount;
        nodeInfo(1,i) = (b-a)/max(a,b);
    end
    sc = mean(nodeInfo);
end

function mp = MP(stringa,stringb)
    global regionalSeries;
    global seriesLength;
    global deviances;
    for i = 1: seriresLength
        a = stringa(i);
        b = stringb(i);


    end



end