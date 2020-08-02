function [ri,jaccard,fm,pa,na,aa,meanvalue,runningtime] = clustering(allSeries,distanceType)
    fprintf('Cluster:\n');
    global k;
    global seriesLength;    %序列长度
    global seriesCount;     %序列个数
    global disttype;
    global Labels;
    disttype = distanceType;
    StartTime = clock;
    distMap=zeros(seriesCount,seriesCount);
     for x = 1:seriesCount
        for y = 1:seriesCount
            distance = DISTANCE(allSeries{x},allSeries{y});
            distMap(x,y) = distance;
            distMap(y,x) = distance;
        end
     end
     cl=linkage(distMap,'average');
    cluSeq=cluster(cl,'maxclust',k);
    
    EndTime = clock;
    a=0;
    b=0;
    c=0;
    d=0;
    for i=1:1:seriesCount
        for j=(i+1):1:seriesCount
            if cluSeq(i)==cluSeq(j)
                if Labels(i)==Labels(j)
                    a=a+1;
                else
                    c=c+1;
                end
            else 
                if Labels(i)==Labels(j)
                    b=b+1;
                else
                    d=d+1;
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
    runningtime = etime(EndTime,StartTime);
end

                
function distance = DISTANCE(stringa, stringb)
    global disttype;
    if disttype==1
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
    end
end


function euclidean =PUREEUCLIDEAN(stringa,stringb)
    euclidean = sqrt(sum(sum((stringa-stringb).^2)));
end

function euclidean = EUCLIDEAN(stringa,stringb)
    euclidean = sqrt(sum(sum((stringa-stringb).^2)));       %??
end

%% 对每个系数分别计算DTW
function dist = myDTW(stringa,stringb)
    global dimension;
    dist = 0;
    for tempDimension = 1:dimension
        dist = dist + dtw(stringa(tempDimension,:),stringb(tempDimension,:));
    end

end


function dist = DTW(ts1,ts2)
    global para1;
    global para2;
    p11 = para1(ts1,:);
    p12 = para1(ts2,:);
    p21 = para2(ts1,:);
    p22 = para2(ts2,:);
    %去除0元素
    p11(find(p11==0))=[];
    p12(find(p12==0))=[];
    p21(find(p21==0))=[];
    p22(find(p22==0))=[];
    %正则化
    p11=mapminmax(p11);
    p12=mapminmax(p12);
    p21=mapminmax(p21);
    p22=mapminmax(p22);
    %分别计算距离
    dist1 = dtw(p11,p12);
    dist2 = dtw(p21,p22);
    dist = dist1 + dist2;
end


function pearson=Pearson(stringa,stringb)
    global seriesLength;
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
    global seriesLength;
    global segmentSize;
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

