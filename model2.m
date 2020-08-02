clear all;
global ts;
global label;
global currentTs;
global errorthreshold;
global segmentMx;
global funcParaMx;
global tempSegment;
global tempSegmentIndex;
% global fitValueDisplay;
% global X;
global para1;
global para2;
global para3;
global distMx;

ts = load('Plane_TRAIN.tsv');
label = ts(:,1);
ts(:,1) = [];
errorthreshold = 0.5;
para1=zeros(size(ts,1),20);
para2=zeros(size(ts,1),20);
para3=zeros(size(ts,1),20);

%% 整体过程
for currentTs=1:1:size(ts,1)
    segmentMx = getSegment(currentTs);
    para1(currentTs,1:size(funcParaMx,1)) = funcParaMx(:,1);
    para2(currentTs,1:size(funcParaMx,1)) = funcParaMx(:,2);
    para3(currentTs,1:size(funcParaMx,1)) = funcParaMx(:,3);
end

for i=1:1:size(ts,1)
    for j=i:1:size(ts,1)
        distMx(i,j) = calcDist(i,j);
    end
end

% currentTs=40;
% segmentMx = getSegment(currentTs);
% for in=1:size(ts(currentTs,:),2) 
%     X=[X,in];
% end
% plot(X,ts(currentTs,:),'r');
% grid on
% hold on
% plot(X,fitValueDisplay,'g');


%% 时间序列分段函数，输出为一个矩阵，每行为一个子序列
function segment = getSegment(i)
    global ts;
    global errorthreshold;
    global funcParaMx;
    global tempSegment;
    global tempSegmentIndex;
%     global fitValueDisplay;
    segment=[];
    funcParaMx=[];
    tempSegment=[];
    tempSegmentIndex=[];
%     fitValueDisplay=[];
    for j=1:size(ts(i,:),2)
        temp = ts(i,j);
        tempSegment = [tempSegment,temp];
        tempSegmentIndex = [tempSegmentIndex,j];
        tempFit = fitting(tempSegment,tempSegmentIndex);
        errorValue = getError(tempSegment,tempSegmentIndex,tempFit);
%         fitValueDisplay = [fitValueDisplay,polyval(tempFit,j)];
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
    SSE = 0;
    avg = mean(value,2);
    length = size(index(1,:),2);
    for i=1:length
        fitValue = polyval(func,index(1,i));
        SSE = SSE + (fitValue - ts(currentTs,index(1,i)))^2;
    end
    errorValue = sqrt(SSE)/abs(avg);
end

%% 对每个系数分别计算DTW
function dist = calcDist(ts1,ts2)
    global para1;
    global para2;
    global para3;
    p11 = para1(ts1,:);
    p12 = para1(ts2,:);
    p21 = para2(ts1,:);
    p22 = para2(ts2,:);
    p31 = para3(ts1,:);
    p32 = para3(ts2,:);
    %去除0元素
    p11(find(p11==0))=[];
    p12(find(p12==0))=[];
    p21(find(p21==0))=[];
    p22(find(p22==0))=[];
    p31(find(p31==0))=[];
    p32(find(p32==0))=[];
    %正则化
    p11=mapminmax(p11);
    p12=mapminmax(p12);
    p21=mapminmax(p21);
    p22=mapminmax(p22);
    p31=mapminmax(p31);
    p32=mapminmax(p32);
    %分别计算距离
    dist1 = dtw(p11,p12);
    dist2 = dtw(p21,p22);
    dist3 = dtw(p31,p32);
    dist = dist1 + dist2 +dist3;
end