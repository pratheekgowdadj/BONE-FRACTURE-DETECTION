Irgb = imread('dataset\new.jpg');
Igray = rgb2gray(Irgb);
I = im2double(Igray);
Gx = [-1 1];
Gy = Gx';
Ix = conv2(I,Gx,'same');
Iy = conv2(I,Gy,'same');
edgeFIS = mamfis('Name','edgeDetection');
edgeFIS = addInput(edgeFIS,[-1 1],'Name','Ix');
edgeFIS = addInput(edgeFIS,[-1 1],'Name','Iy');
sx = 0.1;
sy = 0.1;
edgeFIS = addMF(edgeFIS,'Ix','gaussmf',[sx 0],'Name','zero');
edgeFIS = addMF(edgeFIS,'Iy','gaussmf',[sy 0],'Name','zero');
edgeFIS = addOutput(edgeFIS,[0 1],'Name','Iout');
wa = 0.1;
wb = 1;
wc = 1;
ba = 0;
bb = 0;
bc = 0.7;
edgeFIS = addMF(edgeFIS,'Iout','trimf',[wa wb wc],'Name','white');
edgeFIS = addMF(edgeFIS,'Iout','trimf',[ba bb bc],'Name','black');
r1 = "If Ix is zero and Iy is zero then Iout is white";
r2 = "If Ix is not zero or Iy is not zero then Iout is black";
edgeFIS = addRule(edgeFIS,[r1 r2]);
edgeFIS.Rules
Ieval = zeros(size(I));
for ii = 1:size(I,1)
    Ieval(ii,:) = evalfis(edgeFIS,[(Ix(ii,:));(Iy(ii,:))]');
end
figure
image(Ieval,'CDataMapping','scaled')
colormap('gray')
title('Edge Detection Using Fuzzy Logic')

%%hough transform

[H,T,R] = hough(edge,'RhoResolution',0.5,'Theta',-90:2:89.5);
maxHough = max(H, [], 1);
HoughThresh = (max(maxHough) - min(maxHough))/2 + min(maxHough);
[~, HoughPeaks] = findpeaks(maxHough,'MINPEAKHEIGHT',HoughThresh, 'MinPeakDistance', MinHoughPeakDistance);

% Plot Hough detection results
figure(1)
plot(T, maxHough);
hold on
plot([min(T) max(T)], [HoughThresh, HoughThresh], 'r');
plot(T(HoughPeaks), maxHough(HoughPeaks), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
hold off
xlabel('Theta Value'); ylabel('Max Hough Transform');
legend({'Max Hough Transform', 'Hough Peak Threshold', 'Detected Peak'});

% Locate site of break
if numel(HoughPeaks) > 1
    BreakStack = zeros(size(img, 1), size(img, 2), numel(HoughPeaks));
    % Convolute edge image with line of detected angle from hough transform
    for m = 1:numel(HoughPeaks)

        boneKernel = strel('line', HoughConvolutionLength, T(HoughPeaks(m)));
        kern = double(bwmorph(boneKernel.getnhood(), 'dilate', HoughConvolutionDilate));
        BreakStack(:,:,m) = imfilter(edgeImg, kern).*edgeImg;
    end

    % Take difference between convolution images.  Where this crosses zero
    % (within tolerance) should be where the break is.  Have to filter out
    % regions elsewhere where the bone simply ends.
    brImg = abs(diff(BreakStack, 1, 3)) < BreakLineTolerance*max(BreakStack(:)) & edgeImg > 0;
    [BpY, BpX] = find(abs(diff(BreakStack, 1, 3)) < BreakLineTolerance*max(BreakStack(:)) & edgeImg > 0);
    brImg = bwmorph(brImg, 'dilate', breakPointDilate);
    brReg = regionprops(brImg, 'Area', 'MajorAxisLength', 'MinorAxisLength', ...
        'Orientation', 'Centroid');
    brReg(vertcat(brReg.Area) ~= max(vertcat(brReg.Area))) = [];

    % Calculate bounding ellipse
    brReg.EllipseCoords = zeros(100, 2);
    t = linspace(0, 2*pi, 100);
    brReg.EllipseCoords(:,1) = brReg.Centroid(1) + brReg.MajorAxisLength/2*cos(t - brReg.Orientation);
    brReg.EllipseCoords(:,2) = brReg.Centroid(2) + brReg.MinorAxisLength/2*sin(t - brReg.Orientation);

else
    brReg = [];

end

% Draw ellipse around break location
figure(2)
imshow(img)
hold on
colormap('gray')
if ~isempty(brReg)
    plot(brReg.EllipseCoords(:,1), brReg.EllipseCoords(:,2), 'r');
end
hold off


figure(3)
imshow(leval)