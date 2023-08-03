
year = 2008;
pathString = ['./', num2str(year) '/'];

rawData = struct;

% delta_t stores the time difference between observations in minutes.
rawData.delta_t = importdata([pathString, 'delta_t.mat']);

% x and y are in meters from the upper-left corner of the images. Note that
% positive y-values are the distance down from the top of the image. The
% rows are the floes and the columns are the observations.
rawData.x = importdata([pathString, 'x3.mat']);
rawData.y = importdata([pathString, 'y3.mat']);

% These will be rewritten later once the dataset is filtered.
nFloes = size(rawData.x, 1);
nObs = size(rawData.x, 2);

% These will be rewritten later once the dataset is filtered.
tMinutesArray = nan(1, nObs);
xMetersArray = nan(nFloes, nObs);
yMetersArray = nan(nFloes, nObs);

% Note that delta_t has nObs-1 elements.
tMinutesArray(1) = 0;
for iObs = 1:nObs-1
    tMinutesArray(iObs+1) = tMinutesArray(iObs) + rawData.delta_t(iObs);
end

for iFloe = 1:nFloes
    for iObs = 1:nObs
        metersPerPixel = 250;
        xMetersArray(iFloe, iObs) = metersPerPixel*rawData.x(iFloe, iObs);
        yMetersArray(iFloe, iObs) = metersPerPixel*rawData.y(iFloe, iObs);
    end
end

properties = importdata([pathString, 'prop.mat']);  
prop_order = importdata([pathString, 'PROP_order.mat']);  
FLOE_LIBRARY = importdata([pathString, 'FLOE_LIBRARY_order.mat']);  
for i = 1:size(rawData.x,1)
    if any(rawData.x(i,:))
    day1=find(rawData.x(i,:)~=0);
    spaces=size(day1,2);
    
        for f= 1:spaces
        
        prop=properties{1,day1(f)};    
        [row1,col1]=find(rawData.x(i,day1(f))==prop(:,6) & rawData.y(i,day1(f))==prop(:,7));
        
        
        row11=find(prop(row1,2)==prop_order{1,day1(f)}(:,2));  
        mask11=FLOE_LIBRARY{row11,day1(f)};
        MASK{i,day1(f)}=mask11;
           
        end   
    end
end

radiiArray = nan(nFloes);
for iFloe = 1:nFloes
    dayIdx = find(rawData.x(iFloe, :)~=0);
    if isempty(dayIdx)
        continue
    end
    firstDayIdx = dayIdx(1);
    floeImage = MASK{iFloe, firstDayIdx};

    floePixelSize=250;
    A = floeImage;
    xA=1:size(A,2); yA=1:size(A,1);
    [xA,yA]=meshgrid(floePixelSize*xA,floePixelSize*yA); %% floePixelSize m is the pixel size from Earthdata NASA webpage
    xc=mean(xA(A)); yc=mean(yA(A));
    r=sqrt((xA-xc).^2+(yA-yc).^2);
    r_max=max(r(A));
    
    radiiArray(iFloe) = r_max;
end


% Only keep one observation per day.
keepObsIdx = 1:2:nObs;
xMetersArray = xMetersArray(:, keepObsIdx);
yMetersArray = yMetersArray(:, keepObsIdx);
tMinutesArray = tMinutesArray(keepObsIdx);
nObs = size(xMetersArray, 2);


% Remove empty trajectories and length 1 trajectories.
keepFloeIdx = [];
minObs = 2;
for iFloe = 1:nFloes
    if nnz(xMetersArray(iFloe, :)~=0) >= 2
        keepFloeIdx = [keepFloeIdx, iFloe];
    end
end

xMetersArray = xMetersArray(keepFloeIdx, :);
yMetersArray = yMetersArray(keepFloeIdx, :);
radiiArray = radiiArray(keepFloeIdx);
nFloes = size(xMetersArray, 1);


save('ProcessedData',  'xMetersArray', 'yMetersArray', 'tMinutesArray', "nFloes", 'nObs', "radiiArray")


