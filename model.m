clear all;
global ts;
global label;
global currentTs;
global errorthreshold;
global segmentMx;
global funcParaMx;
global tempSegment;
global tempSegmentIndex;
global fitValueDisplay;
global X;

ts = load('Meat_TRAIN.tsv');
label = ts(:,1);
ts(:,1) = [];
errorthreshold = 0.5;

%% 整体过程
% for currentTs=1:size(ts,1)
%     segmentMx = getSegment(currentTs);
% end

currentTs=1;
segmentMx = getSegment(currentTs);
for in=1:size(ts(currentTs,:),2) 
    X=[X,in];
end
plot(X,ts(currentTs,:));
grid on
hold on
plot(X,fitValueDisplay);


%% 时间序列分段函数，输出为一个矩阵，每行为一个子序列
function segment = getSegment(i)
    global ts;
    global errorthreshold;
    global funcParaMx;
    global tempSegment;
    global tempSegmentIndex;
    global fitValueDisplay;
    segment=[];
    tempSegment=[];
    tempSegmentIndex=[];
    for j=1:size(ts(i,:),2)
        temp = ts(i,j);                                 
        tempSegment = [tempSegment,temp];               %当前段增加一个时间点的数据
        tempSegmentIndex = [tempSegmentIndex,j];        %保存当前段下标
        tempFit = fitting(tempSegment,tempSegmentIndex);%拟合
        errorValue = getError(tempSegment,tempSegmentIndex,tempFit);
        fitValueDisplay = [fitValueDisplay,polyval(65,j)]%j下标下的拟合的结果
        if errorValue > errorthreshold
            segment = [segment;j];
            funcParaMx = [funcParaMx;tempFit];
            tempSegment=[];
            tempSegmentIndex=[];
        end
    end
end

%% 动态拟合过程
function fit = fitting(value,index)
    fit = polyfit(index,value,2);
end

%% 误差函数定义
function errorValue = getError(value,index,func)
    global ts;
    global currentTs;
    global fitValueDisplay;
    SSE = 0;
    avg = mean(value,2);
    length = size(index(1,:),2);
    for i=1:length
        fitValue = polyval(func,index(1,i));
        SSE = SSE + (fitValue - ts(currentTs,index(1,i)))^2;
    end
    errorValue = sqrt(SSE)/abs(avg);
end



