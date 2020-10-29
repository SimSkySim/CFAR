% OS CFAR 2D
% October 20th 2020
% Simyon Pinsky

clear all
close all

%Parameters
N = 141;%Raws
M = 151;%Colomn
detector.ProbabilityFalseAlarm = 10e-5;%Probablity of false alarm;
detector.GuardBandSize = [0 0];
detector.TrainingBandSize = [6 4];
detector.Method = 'OSCA';
plot_en = 1;

% Create a N-by-M image containing random complex data. Then, square the
% data to simulate a square-law detector.
rs = RandStream.create('mt19937ar','Seed',5);
x = 2/sqrt(2)*(randn(rs,N,M) + 1i*randn(rs,N,M));
x2 = abs(x).^2;

% Process all the cells in each image. To do this, find the row and column
% of each CUT cell whose training region falls entirely within each image.
Mgc = detector.GuardBandSize(2);
Ngr = detector.GuardBandSize(1);
Mtc = detector.TrainingBandSize(2);
Ntr = detector.TrainingBandSize(1);
cutidx = [];
colstart = Mtc + Mgc + 1;
colend = M - ( Mtc + Mgc);
rowstart = Ntr + Ngr + 1;
rowend = N - ( Ntr + Ngr);
for m = colstart:colend-1
    for n = rowstart:rowend-1
        cutidx = [cutidx,[n;m]];%Generate all the cut indexes.
    end
end

if plot_en
    to_be_analyzed = zeros(size(x2(:,:,1)));
    to_be_analyzed(cutidx(1,:),cutidx(2,:)) = 1;
    figure
    subplot(2,1,1)
    surf(x2(:,:,1))
    view(2)
    colorbar
    xlabel('M')
    ylabel('N')
    title('Image')
    xlim([1 M])
    ylim([1 N])
    axis equal
    subplot(2,1,2)
    surf(to_be_analyzed)
    view(2)
    colorbar
    xlabel('M')
    ylabel('N')
    xlim([1 M])
    ylim([1 N])
    title('Area to be processed by the CFAR')
    axis equal
end

% Perform the detection on all CUT cells. Return the detection classification and the
% threshold used to classify the cell.
[dets,th,noise,loc] = apply_detector(x2,cutidx,detector,N,M);
 
% Compute the Empirical Pfa
ncutcells = size(cutidx,2);
pfa_emp = sum(dets(:))/ncutcells
% The empirical and specified pfa must agree.

% Display the average empirical threshold value over all images.
mean(th(:))

% Compute the theoretical threshold factor for the required pfa.
threshfactor = npwgnthresh(detector.ProbabilityFalseAlarm,1,'noncoherent');
threshfactor = 10^(threshfactor/10);
disp(threshfactor)

% The theoretical threshold factor multiplied by the noise variance should
% agree with the measured threshold.
noisevar = mean(x2(:));
disp(threshfactor*noisevar);
% The theoretical threshold and empirical threshold agree to within
% an acceptable difference.

%Rearange and plot the results.
detimg = zeros(N,M);
th_val = zeros(N,M);
noise_val = zeros(N,M);
for k = 1:ncutcells
    detimg(cutidx(1,k),cutidx(2,k)) = dets(k);
    th_val(cutidx(1,k),cutidx(2,k)) = th(k);
    noise_val(cutidx(1,k),cutidx(2,k)) = noise(k);
end

figure
subplot(2,2,1)
imagesc(x2)
set(gca,'YDir','normal')
tstr = sprintf('CFAR Input N = %3.0f M = %3.0f',N,M');
title(tstr)
colorbar
subplot(2,2,2)
imagesc(loc)
set(gca,'YDir','normal')
tstr = sprintf('Detections Mgc = %3.3f Ngr = %3.3f Mtc = %3.3f Ntr = %3.3f',Mgc,Ngr,Mtc,Ntr');
title(tstr)
colorbar
subplot(2,2,3)
imagesc(th_val)
set(gca,'YDir','normal')
tstr = sprintf('Thresholds for Pfa = %1.3f',detector.ProbabilityFalseAlarm);
title(tstr)
colorbar
subplot(2,2,4)
imagesc(detimg)
set(gca,'YDir','normal')
tstr = sprintf('CUT pfa-emp = %3.5f Pfa = %3.5f Av-Emp-Thr = %3.3f TheTthXNoiseVar = %3.3f',pfa_emp,detector.ProbabilityFalseAlarm,mean(th(:)),threshfactor*noisevar');
title(tstr)
colorbar

function [dets,th,noise,loc] = apply_detector(x2,cutidx,detector,N,M)
numCUT = size(cutidx,2);
th = zeros(numCUT,1);
noise = zeros(numCUT,1);


for k = 1:numCUT
    loc = zeros(size(x2));
    if k == 1
        [trnindsl,trnindsr] = get2DTrainingInds(detector,x2,cutidx(:,k));%Same function as for GO and SO CFARs.
        Nc = length(trnindsl)+length(trnindsr);
        Pfa = detector.ProbabilityFalseAlarm;
        Rank = floor(Nc*3/4);
        ThFac = req_alpha_calc(Pfa,Nc,Rank);%Threshold factor Calculation
    else
        dIdx = cutidx(:,k)-cutidx(:,k-1);
        trnindsl = trnindsl  + dIdx(1) + size(x2,1)*dIdx(2);
        trnindsr = trnindsr  + dIdx(1) + size(x2,1)*dIdx(2);
    end
    % Form noise power estimate
    loc(trnindsl)=1;
    loc(trnindsr)=2;
    loc(cutidx(1,k),cutidx(2,k))=5;
    if (k==1)
        figure
        subplot(2,1,1)
        surf(loc);
        view(2)
        xlabel('M')
        ylabel('N')
        title('First CUT')
        xlim([1 M])
        ylim([1 N])
        axis equal
    end
    if (k==numCUT)
        subplot(2,1,2)
        surf(loc);
        view(2)
        xlabel('M')
        ylabel('N')
        title('Last CUT')
        xlim([1 M])
        ylim([1 N])
        axis equal
    end
    temp = sort([x2(trnindsl); x2(trnindsr)],1,'ascend');
    noisePowEst = temp(Rank);
    % Form the threshold for this CUT
    th(k) = noisePowEst * ThFac;
    noise(k) = noisePowEst;
end

linIdx = sub2ind(size(x2),cutidx(1,:)',cutidx(2,:)');
y = x2(linIdx) > th;
dets=y;

end

function varargout = get2DTrainingInds(obj,X,Idx)
%   get2DTrainingIdx Calculate the 2-D training cells indices
%   [FrontIndices, RearIndices]=get2DTrainingIdx(htr,Ncells,Idx)
%   returns the training cell indices in front of the cell under
%   test (FrontIndices) and after the cell under test
%   (RearIndices).

[ NumRows, NumColumns]  = size(X);

GuardRegionSize = getGuardRegionSize(obj);
TrainingRegionSize = getTrainingRegionSize(obj);

% Compute Training and Guard cell indices
indx1_r = linspace(Idx(1)-(TrainingRegionSize(1)-1)/2, ...
    Idx(1)+(TrainingRegionSize(1)-1)/2,TrainingRegionSize(1));
indx1_c = linspace(Idx(2)-(TrainingRegionSize(2)-1)/2, ...
    Idx(2)+(TrainingRegionSize(2)-1)/2,TrainingRegionSize(2));
indx2_r = linspace(Idx(1)-(GuardRegionSize(1)-1)/2, ...
    Idx(1)+(GuardRegionSize(1)-1)/2,GuardRegionSize(1));
indx2_c = linspace(Idx(2)-(GuardRegionSize(2)-1)/2, ...
    Idx(2)+(GuardRegionSize(2)-1)/2,GuardRegionSize(2));

% Generate the subscripts of the Guard and Training Regions
[indx1_C, indx1_R] = meshgrid(indx1_c,indx1_r);
[indx2_C, indx2_R] = meshgrid(indx2_c,indx2_r);

% Convert subscripts to linear indices
linIndx1 = sub2ind(size(X),indx1_R(:),indx1_C(:));
linIndx2 = sub2ind(size(X),indx2_R(:),indx2_C(:));
linIndxCUT = sub2ind(size(X),Idx(1),Idx(2));

% Compute indices of training cells, excluding the guard region and
% CUT.
TrainingCellIndx = setdiff(linIndx1,linIndx2);
chIdx = 1;
%Front training region are cells with indices less than CUT
FrontTrainingCellIdx = TrainingCellIndx(TrainingCellIndx > linIndxCUT);
RearTrainingCellIdx  = TrainingCellIndx(TrainingCellIndx < linIndxCUT);
FrontTrainingCellIdx = repmat(FrontTrainingCellIdx,1,size(chIdx,2)) + ...
    repmat((chIdx-1)*(NumRows*NumColumns),size(FrontTrainingCellIdx,1),1);
RearTrainingCellIdx = repmat(RearTrainingCellIdx,1,size(chIdx,2)) + ...
    repmat((chIdx-1)*(NumRows*NumColumns),size(RearTrainingCellIdx,1),1);
varargout{1} = FrontTrainingCellIdx;
varargout{2} = RearTrainingCellIdx;

end

function GuardRegionSize = getGuardRegionSize(obj)
if isscalar(obj.GuardBandSize)
    GuardRegionSize = [2*obj.GuardBandSize + 1 ...
        2*obj.GuardBandSize + 1];
else
    GuardRegionSize = 2*obj.GuardBandSize + 1;
end
end

function TrainingRegionSize = getTrainingRegionSize(obj)
GuardRegionSize = getGuardRegionSize(obj);
TrainingRegionSize = 2*obj.TrainingBandSize + GuardRegionSize;
end

function req_alpha = req_alpha_calc(req_pfa,Ntc,Rank)
%Threshold Factor Calculation
alpha = 4:20;
a=factorial(Ntc);
b=factorial(alpha+Ntc-Rank);
c=factorial(Ntc-Rank);
d=factorial(alpha+Ntc);
Pfa = (a*b)./(c*d);
req_pfa_log = log10(req_pfa);
[val ind]=min(abs((log10(Pfa)-req_pfa_log)));
req_alpha = alpha(ind);
if 1
    figure
    plot(alpha,log10(Pfa),'r')
    hold on
    plot(alpha(ind),log10(Pfa(ind)),'ok')
    ylim([-6 -2])
    title('Thr Fact')
    xlabel('Alpha')
    ylabel('log10(Pf)')
end
end

