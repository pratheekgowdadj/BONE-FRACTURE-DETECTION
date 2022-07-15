

im=imread('dataset\new.jpg');
im=im2double(im);
I1 = rgb2gray(im);
%smoothening the image with a filter
gfilter= [0 0 1 0 0;
       0 1 2 1 0;
       1 2 -16 2 1;
       0 1 2 1 0;
       0 0 1 0 0];
   
smim=conv2(I1,gfilter);


% finding the zero crossings

[rr,cc]=size(smim);
zc=zeros([rr,cc]);

for i=2:rr-1
    for j=2:cc-1
        if (smim(i,j)>0)
             if (smim(i,j+1)>=0 && smim(i,j-1)<0) || (smim(i,j+1)<0 && smim(i,j-1)>=0)
                             
                zc(i,j)= smim(i,j+1);
                        
            elseif (smim(i+1,j)>=0 && smim(i-1,j)<0) || (smim(i+1,j)<0 && smim(i-1,j)>=0)
                    zc(i,j)= smim(i,j+1);
            elseif (smim(i+1,j+1)>=0 && smim(i-1,j-1)<0) || (smim(i+1,j+1)<0 && smim(i-1,j-1)>=0)
                  zc(i,j)= smim(i,j+1);
            elseif (smim(i-1,j+1)>=0 && smim(i+1,j-1)<0) || (smim(i-1,j+1)<0 && smim(i+1,j-1)>=0)
                  zc(i,j)=smim(i,j+1);
            end
                        
        end
            
    end
end


otpt=im2uint8(zc);
% thresholding
otptth= otpt>105;

%%figure;
 % subplot(2,2,1);imshow(im);title('Origional image');
  %subplot(2,2,2);imshow(smim);title('Smoothened image');
 % subplot(2,2,3);imshow(otpt);title('Output image');
 % subplot(2,2,4);imshow(otptth);title('Output image with threshold');%%

  % final result
  % figure, imshow(otptth);


% Do hough transform on edge image to find angles at which bone pieces are
% found
% Use max value of Hough transform vs angle to find angles at which lines
% are oriented.  If there is more than one major angle contribution there
% will be two peaks detected but only one peak if there is only one major
% angle contribution (ie peaks here = number of located bones = Number of
% breaks + 1)
[H,T,R] = hough(otpt,'RhoResolution',1,'Theta',-90:2:89.5);
maxHough = max(H, [], 1);
HoughThresh = (max(maxHough) - min(maxHough))/2 + min(maxHough);
[~, HoughPeaks] = findpeaks(maxHough,'MINPEAKHEIGHT',HoughThresh);

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



   


