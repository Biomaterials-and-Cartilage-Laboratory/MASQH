% The 'homogeneity_one_image' function computes the homogeneity index of a
% single image through the MASQH method. It allows users to chose an image
% from a folder and returns a single index between 0 and 1.
%The functions takes one input : 'exclude_zone'. This parameter should be 1
% if zones need to be exluded and 0 if the image should be used without any
% zone exluded. Default value is 0.

% The algorithm is described and analysed in : F. Milano, A. Chevrier, 
% G. De Crescenzo and M. Lavertu, "Robust segmentation-free algorithm for 
% homogeneity quantification in images," in IEEE Transactions on Image 
% Processing, doi: 10.1109/TIP.2021.3086053.

function  [Index] = homogeneity_one_image(exclude_zone) %,tolerance, Chitosan, Erythrocites )
close all % close open matlab images

% Path of image to analyze
[path, ~] = imgetfile;
Initial_Picture_Color = imread(path);

%Initial_Picture_Color = imread('/Users/fionamilano/Documents/lab/2. Homogeneity_quantification/image_database/JAB8833 10x_zoom40.tif');

[r, c]=size(Initial_Picture_Color); % length and width of image to analyze

if nargin == 0 || isempty(exclude_zone)
    exclude_zone = 0;
    binaryImage_bis = zeros(r,c);
end

%% Initialisation of the main algorithm variables
% nbr_div = the number of square sizes used in the analysis (= k in the
% associated article)
nbr_div = 100;
% number of pixels taken per square ( = x in the associated article)
nbr_echantillon = 5000;
% number of squares used for the analysis ( = n in the associated article)
nbr_carre = 400;

%% Indices and vectors initialisation (to contain square sizes, p-values, proportions, ...)
% Creation  of the vector containing the various squares' lengths
vect = linspace(round(r/4),50,nbr_div);
vect = 2*round(vect);

% vector containig pourcentage of p-val > 0.05
pourcent = zeros(1,nbr_div);

% initialisation other indices
iter = 1;


%% Image channel extraction
% Take only the green channel of the image (can be changed)
Initial_Picture=im2double(Initial_Picture_Color(:,:,2));
imshow(Initial_Picture)

if exclude_zone == 1
    %% Selection the zones to be excluded if needed
    imshow(Initial_Picture)
    i = 1;
    finished = 'NO';
    binaryImage_bis = zeros(r,c);
    
    Start = questdlg('Do you want to remove some zone(s)?', ...
        'confirmation', ...
        'YES', 'NO', 'NO');
    if strcmpi(Start,'NO')
        finished = 'YES';
    else
        while strcmpi(finished,'NO')
            h(i) = imfreehand();
            finished = questdlg('Finished?', ...
                'confirmation', ...
                'YES', 'NO', 'UNDO', 'NO');
            if strcmpi(finished, 'UNDO')
                delete(h(i))
                finished = 'NO';
            else
                binaryImage = h.createMask();
                binaryImage_bis = binaryImage_bis+binaryImage;
                imshow(imcomplement(binaryImage_bis).*Initial_Picture);
            end
        end
    end
    
    % image avec zones mises en évidence
    imshow(Initial_Picture_Color);
    hold on
    blue = cat(3, zeros(r,c), zeros(r,c), ones(r,c));
    dark = imshow(blue);
    hold off
    set(dark, 'AlphaData', binaryImage_bis.*0.3);
end

% the rest of the code works in uint8
initial_picture = im2uint8(Initial_Picture);


%% Analysis at 'nbr_div' different box size
close all;

for i = vect % for each square size
    
    % Radom square placement and withdrawal of associated histograms
    [subhist,not_possible] = subimage(i, nbr_carre, initial_picture, binaryImage_bis, nbr_echantillon);
    
    if not_possible ==  0
        p_val = zeros(nbr_carre/2,1);
        for k = 1:2:nbr_carre
            % the histograms of two squares are compared. if they come from
            % different distributions, pval should be < 0.5
            sub_hist_1 = subhist(k,:);
            sub_hist_2 = subhist(k+1,:);
            
            % statistical test per se
            [~,p] = kstest2(sub_hist_1, sub_hist_2);
            p_val((k+1)/2) = p;
        end
        pourcent(iter) = sum(p_val>0.05)/length(p_val); % computation of pourcentage of p-val > 0.05
    end
    iter = iter+1;
end

Index = sum(pourcent)/100;

figure
plot(pourcent)
%legend({'Proportion on  p-value > 0.05'}, 'FontSize',14)
xlabel('iteration', 'FontSize', 14)
ylabel('proportion', 'FontSize', 14)
title('\fontsize{16} Evolution of the proportion of p-values  > 0.05')

figure
plot(movmean(pourcent,10))
%legend({'Proportion on  p-value > 0.05'}, 'FontSize',14)xlabel('iteration', 'FontSize', 14)
ylabel('proportion', 'FontSize', 14)
title('\fontsize{16} Evolution of the proportion of p-values > 0.05 (moving average)')

fprintf('|Index of analyzed image =    | %.2f| \n', Index)

end

% The subimage fonction randomly chose a number of points (nbr_echantillon)
% in n randomly positioned squares in the image 'image'. It output one
% histogram of the chosen points for every positioned squares.
% At the same time, it verifies that the chosen square are not placed in
% the part of the picture that should not be analysed (through the mask
% 'binaryImage'). If a square is misplaced, the function tries to place it
% elsewhere. After trying 10 times, if no square are found, the actual
% square size is skipped.
function [sub_hist, not_possible] = subimage(c, n, image, binaryImage, nbr_echantillon)
[row, col] = size(image);
sub_hist = zeros(n,255);
max_iter = 10;

couples = zeros(nbr_echantillon, 2);
couples(:,1) = round(1+rand([1 nbr_echantillon]).*(c-1));
couples(:,2) = round(1+rand([1 nbr_echantillon]).*(c-1));
not_possible = 0;

i = 1;
while i <= n
    nbr_ejected = 0;
    centres = rand(2,1); % n° of row and col of the square center
    centres(1) = round(1+centres(1).*(row-c-1));
    centres(2) = round(1+centres(2).*(col-c-1));
    for j = 1:nbr_echantillon
        val = image(centres(1)+couples(j,1), centres(2)+couples(j,2));
        nbr_ejected = nbr_ejected+binaryImage(centres(1)+couples(j,1), centres(2)+couples(j,2));
        sub_hist(i,val+1) = sub_hist(i,val+1)+1;
    end
    %max_eject = round(nbr_echantillon/20);
    if nbr_ejected==0
        i = i+1;
        max_iter = 0;
    else
        sub_hist(i,:) = 0;
        max_iter = max_iter+1;
        if max_iter > 50
            not_possible = 1;
            return
        end
    end
end

end




