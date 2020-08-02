function none = Entry()
    clear all;
    %   ????
    none = 1;
    global FileName;
    global cheat;
    global RawTimeSeries;   %原始
    global Variance;        %方差
    global Variances;       %方差
    global DisturbedFile;   %?????????????始数据
    global disturbedSeries; %扰乱的数据
    global seriesLength;    %序列长度
    global seriesCount;     %序列个数
    global segmentSize;     %分段大小
    global Labels;          %标签
    global k;               %?????
    global dimension;       %多项式维度
    global fitSeries;
    global reconstructedSeries;
    global method;          %1:取所有的均值进行拟合
                            %2:随机选择+均值
    global errorThreshold;
    global fuzzyThres;
    %-----数据准备-----%
     
    %   聚类方法
    kmeans = true;
    dp = false;
    cheat = false;
    
    dimensionFloor= 1;
    dimensionCeil = 3;
    segmentSizeFloor = 3;
    segmentSizeCeil = 24;
    segmentSizeStep = 7;
    repeatClusterTime = 2;
    %   数据集文件
    FileName = 'data/FaceFour_TEST';
    k = 4;
    segmentSize = 6;
    dimension = 2;
    fuzzyThres = 150;
    %   载入数据集
    RawTimeSeries = load(FileName);
    OriginalTimeSeries = load(FileName);
    [seriesCount,seriesLength] = size(RawTimeSeries(:,:));
    random = randperm(ceil(seriesCount*0.5));
    for i = 1:seriesCount*0.5
        RawTimeSeries(random(i),:)=[];
    end
    [seriesCount,seriesLength] = size(RawTimeSeries(:,:));
    random = randperm(ceil(seriesCount*0.5));
    for i = 1:seriesCount*0.5
        RawTimeSeries(random(i),:)=[];
    end

    %   ??????????
    Labels = RawTimeSeries(:,1);
    RawTimeSeries(:,1)=[];
    [seriesCount,seriesLength] = size(RawTimeSeries(:,:));
    RawTimeSeriesCell = {};
    for i = 1:seriesCount 
        RawTimeSeriesCell{i}(:)=RawTimeSeries(i,:);
    end
    %-----删除结果文件-----%
    delete('KmeansResult');
    delete('DPResult');
    
    e = 0.0;
    d = 0.0;
    for i = 1:10
        for j = 1:seriesLength
            e = e + RawTimeSeries(i,j);
        end
    end
    e = e/(seriesLength*10);
   
    for i = 1:10
        for j = 1:seriesLength
            d = d+(RawTimeSeries(i,j) - e)^2;
        end
    end
    d = d/(seriesLength*10);
    sqrtd = sqrt(d);    
    
    bestMethod = -1;
    bestVariance = -1;
    bestSegment = -1;
    bestSC = -1;

    %未扰乱聚类
    undisturbedResult = zeros(2,13);  
                [ri,jaccard,fm,pa,na,aa,mean,runningtume,sc] = KMEANS(cheat,1,"raw");
                undisturbedResult(2,:) = undisturbedResult(2,:)+ [-1,2,dimension,ri,jaccard,fm,pa,na,aa,mean,runningtume,segmentSize,sc];
%                 [ri,jaccard,fm,pa,na,aa,mean,runningtume] = clustering(RawTimeSeriesCell,6);
                undisturbedResult(1,:) = undisturbedResult(1,:)+ [-1,1,dimension,ri,jaccard,fm,pa,na,aa,mean,runningtume,segmentSize,sc];
                
    %-----Variance Iteration-----%
%     for Variance=0:0.05:0.4
    for Variance=0.1
        for i = 2:2
             TOFILE(undisturbedResult(i,:),2);
        end
        euclideanResult = zeros(4,13);     
        result = zeros(2,2,segmentSizeCeil-segmentSizeFloor+1,dimensionCeil,13);    
            for clustertime = 1:repeatClusterTime
                %   采用随机化的标准差：每一条时间序列的标准差在0——variance中随机生成
                % Variances = unifrnd(Variance*2/3,Variance*4/3,rows,1);
                %   采用统一的标准差：每一条时间序列的标准差相同
                for i = 1:seriesCount
                    Variances(i) = sqrtd*Variance;
                end

                for i = 1:seriesCount
                    for j = 1:seriesLength
                        %   根据方差随机生成五个数据
                        disturbedSeries{i}(:,j)= normrnd(RawTimeSeries(i,j),Variances(i),1,5)';
                    end
                end
                %拟合所有序列
                % fitAllSeries();
                
                %计算距离
                dimension = 0;
                [ri,jaccard,fm,pa,na,aa,mean,runningtume,sc] = KMEANS(cheat,1,"mean");
                euclideanResult(2,:) = euclideanResult(2,:)+ [0,2,dimension,ri,jaccard,fm,pa,na,aa,mean,runningtume,segmentSize,sc];
%                 [ri,jaccard,fm,pa,na,aa,mean,runningtume] = clustering(disturbedSeries,5);
                %[ 分段类型，聚类方法 ]
                euclideanResult(1,:) = euclideanResult(1,:)+ [0,1,dimension,ri,jaccard,fm,pa,na,aa,mean,runningtume,segmentSize,sc];
                
                
                %segmenttype 1 平均分段 2 动态分段
                %distancetype 1 参数欧氏距离 2 原始序列欧式距离 3 dtw
                
                dimension = 1;
                for segmentSize = segmentSizeFloor:segmentSizeStep:segmentSizeCeil
                    dimensionTop = dimensionCeil;
                    if dimensionCeil >= segmentSize
                       dimensionTop = segmentSize-1 ;
                    end
                    for dimension = dimensionFloor:dimensionTop
                        fitAllSeries(1);%平均分段
                            temp = zeros(13,1);
                            %[ri,jaccard,fm,pa,na,aa,mean,runningtume] = clustering(fitSeries,3);
                            temp(:) = result(1,1,segmentSize-segmentSizeFloor+1,dimension-dimensionFloor+1,:);
                            result(1,1,segmentSize-segmentSizeFloor+1,dimension-dimensionFloor+1,:) =temp(:) + [1;1;dimension;ri;jaccard;fm;pa;na;aa;mean;runningtume;segmentSize;sc];
                            [ri,jaccard,fm,pa,na,aa,mean,runningtume,sc] = KMEANS(cheat,1,"fit");
                            if sc>bestSC
                                bestSC = sc;
                                bestMethod = 1;
                                bestSegment = segmentSize;
                                bestVariance = Variance;
                                bestDimension = dimension;
                            end
                            temp(:) =  result(1,2,segmentSize-segmentSizeFloor+1,dimension-dimensionFloor+1,:);
                            result(1,2,segmentSize-segmentSizeFloor+1,dimension-dimensionFloor+1,:) = temp(:) + [1;2;dimension;ri;jaccard;fm;pa;na;aa;mean;runningtume;segmentSize;sc];
%                         fitAllSeries(2);%动态分段        重建
%                             [ri,jaccard,fm,pa,na,aa,mean,runningtume] = clustering(reconstructedSeries,3);
                            temp(:) =  result(2,1,segmentSize-segmentSizeFloor+1,dimension-dimensionFloor+1,:);
                            result(2,1,segmentSize-segmentSizeFloor+1,dimension-dimensionFloor+1,:) = temp(:) + [2;1;dimension;ri;jaccard;fm;pa;na;aa;mean;runningtume;segmentSize;sc];
                            [ri,jaccard,fm,pa,na,aa,mean,runningtume,sc] = KMEANS(cheat,8,"reconstructed");
                            if sc>bestSC
                                bestSC = sc;
                                bestMethod = 2;
                                bestSegment = segmentSize;
                                bestVariance = Variance;
                                bestDimension = dimension;
                            end
                            temp(:) =  result(2,2,segmentSize-segmentSizeFloor+1,dimension-dimensionFloor+1,:);
                            result(2,2,segmentSize-segmentSizeFloor+1,dimension-dimensionFloor+1,:) = temp(:) + [2;2;dimension;ri;jaccard;fm;pa;na;aa;mean;runningtume;segmentSize;sc];
                    end
                end
            end
            
        
            
        for clusterType = 2:2
           TOFILE(euclideanResult(clusterType,:)/repeatClusterTime,2);
           for segmentType = 1:2
                for temp = segmentSizeFloor:segmentSizeStep:segmentSizeCeil     
                    for temp2 = dimensionFloor:dimensionCeil
                        TOFILE(result(segmentType,clusterType,temp-segmentSizeFloor+1,temp2-dimensionFloor+1,:)/repeatClusterTime,3);
                    end
                end
           end
        end
            
       
         % 结果导入excel
        % 绘制三维图
%         picResult(:,:,:) = result(:,:,8);
%         bar(picResult);
%         set(gca,'YDir','normal');  %从下到上递增
%         xlabel('分段长度');
%         ylabel('jaccard');
%         legend('1次','2次','3次','4次','5次','6次');
%         title('face four 多项式次数与分段长度');

    end
       

    RawTimeSeries = OriginalTimeSeries;
    Labels = RawTimeSeries(:,1);
    RawTimeSeries(:,1)=[];
    [seriesCount,seriesLength] = size(RawTimeSeries(:,:));
    RawTimeSeriesCell = {};
    for i = 1:seriesCount 
        RawTimeSeriesCell{i}(:)=RawTimeSeries(i,:);
    end
    for i = 1:seriesCount
        Variances(i) = sqrtd*bestVariance;
    end

    for i = 1:seriesCount
        for j = 1:seriesLength
            %   根据方差随机生成五个数据
            disturbedSeries{i}(:,j)= normrnd(RawTimeSeries(i,j),Variances(i),1,5)';
        end
    end
    segmentSize = bestSegment;
    dimension = bestDimension;
    fitAllSeries(1);%平均分段
    if bestMethod == 1
        [ri,jaccard,fm,pa,na,aa,mean,runningtume,sc] = KMEANS(cheat,1,"fit");
    elseif bestMethod == 2
        [ri,jaccard,fm,pa,na,aa,mean,runningtume,sc] = KMEANS(cheat,1,"reconstructed");
    end
    TOFILE([bestMethod,-1,dimension,ri,jaccard,fm,pa,na,aa,mean,runningtume,segmentSize,sc],1);
    ClusterResult = load('KmeansResult');
     xlswrite('KmeansResult.xlsx',ClusterResult,1,'B3');

end

function none = TOFILE(result,distanceType)
        global Variance;
        global dimension;
        none = 0;
        segmentMethod = result(1);
        clusterMethod = result(2);
        dimension = result(3);
        ri = result(4);
        jaccard = result(5);
        fm = result(6);
        pa = result(7);
        na = result(8);
        aa = result(9);
        mean = result(10);
        runningtume = result(11);
        segmentSize = result(12);
        sc = result(13);
        resultfile = fopen("KmeansResult","a");
        fprintf(resultfile,'%i,%i,%i,%i,%f,%i,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',segmentMethod,clusterMethod,dimension,segmentSize,Variance,distanceType,ri,jaccard,fm,pa,na,aa,mean,runningtume,sc);
        fclose(resultfile);   
end
    
function [ri,jaccard,fm,pa,na,aa,mean,runningtume] = getClusterIndex(ClusterAssignment)
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

