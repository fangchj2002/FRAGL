%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo of "Fuzzy region-based active contours driven by hybrid fitted energy 
% with local and global information for image segmentation"
% Jiangxiong Fang
% East China University of Technology&&Nanchang University, Nanchang, China
% 6th, Jan., 2019
% Email: fangchj2002@163.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;
addpath 'images';
ImgID = 9706;
Img = imread([num2str(ImgID),'.jpg']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add different types of noise:Gaussian/speckle/salt & pepper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Img = imnoise(Img,'Gaussian',0,10^2/255^2);%10^2/25^2
Img = imnoise(Img,'speckle',0.002);%10^2/25^2
Img = imnoise(Img,'salt & pepper',0.01);%10^2/25^2
 
tic;
%setting the initial level set function 'u':
[M,N,L] = size(Img);
u = zeros(M,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the pseudo level set function (LSF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u(:,:) = 0.3;
u(40:60,40:80) = 0.7;

switch ImgID
    case 1
        Img_gray = Img;
        rad = 3;
        iterNum = 40;
        lambda1 = 1;
        lambda2 = 1;
        alpha = .5;
        belta = .5;
        m1 = 1;
        m2 = 1;
    case 2
        Img_gray = rgb2gray(Img);
        iterNum = 40;
        rad = 1;
        lambda1 = 2;
        lambda2 = 1;
        alpha = .5;
        belta = .5;
        m1 = 1;
        m2 = 1;
    case 3
        Img_gray = rgb2gray(Img);
        iterNum = 40;
        rad = 1;
        lambda1 = 1;
        lambda2 = 1;
        alpha = .5;
        belta = .5;
        m1 = 1;
        m2 = 1;
    case 4
        Img_gray = Img;
        iterNum = 40;
        rad = 3;
        lambda1 = 1.8;
        lambda2 = 1;
        alpha = .5;
        belta = .5;
        m1 = 1;
        m2 = 1;
    case 5
        Img_gray = Img;
        iterNum = 40;
        rad = 1;
        lambda1 = 1.5;
        lambda2 = 1;
        alpha = 0.5;
        belta = 0.5;
        m1 = 1;
        m2 = 1;
    otherwise
        Img_gray = rgb2gray(Img);
        rad = 1;
        iterNum = 40;
        lambda1 = 1;
        lambda2 = 1;
        alpha = 0.5;
        belta = .5; 
        m1 = 1;
        m2 = 1;
end

[Ix,Iy] = gradient(double(Img_gray));
f = Ix.^2+Iy.^2;
g = 1./(1+f);  % edge indicator function
diswght = disweight(rad);
if diswght==0
    saliency = Img;
else
    saliency = imfilter(Img_gray,diswght,'replicate');
end

figure;subplot(2,2,1);imshow(Img);hold on;%axis off,axis equal

title('Initial contour');
[c,h] = contour(u-0.5,[0 0],'r','LineWidth',2);

subplot(2,2,2);
imshow(saliency,[]);hold on;
subplot(2,2,3);

energy1 = [];
delF = [];

for n=1:iterNum
    [u,e,deltaF] = FRAGL_v1(double(saliency),u,diswght,lambda1,lambda2,alpha,belta,m1,m2,g);

    if mod(n,5)==0
        pause(0.1);
        imshow(Img, []);hold on;axis off,axis equal
        [c,h] = contour(u-0.5,[0 0],'r','LineWidth',2);
        iterNum=[num2str(n), 'iterations'];
        title(iterNum);
        hold off;
    end
end

seg = ((u-0.5)>0);
subplot(2,2,4),imshow(seg);

totalIterNum=[num2str(n), ' iterations'];
title(['Final contour, ', totalIterNum]);
time = toc;
figure;

mesh(u-0.1);
title('Final level set function');


