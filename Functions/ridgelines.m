function [bgm2,fgm] = ridgelines(U,filename,outidx)

%%% Getting the Watershed Lines (Ridge Lines) to demarcate the different
%%% areas

U1 = imtophat(U,strel('disk',20));%figure , imshow(uint8(I));
%imshow(uint8(I))
U2 = gaussFiltDU(U1);
BW = mat2gray(U2);
%figure, imshow(uint8(BW));

%figure , imshow((mat2gray(U2,[70 80])).*mask1(:,1:2595));
%figure , imshow(imerode((mat2gray(U2,[70 80])).*mask1(:,1:2595)),se);

N = 10;
BW = [BW(:,end-(2*N-1):end),BW,BW(:,1:2*N)];

% BW = BW.*mask;
% figure, imshow(uint8(BW));

gmag = imgradient(BW);


se  = strel('disk',5);
img = mat2gray(BW);
Ie = imerode(img,se);
Iobr = imreconstruct(Ie,img);
%figure, imshow((Iobr)), title('Opening-by-Reconstruction');

Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
%figure, imshow(Iobrcbr), title('Opening-Closing by Reconstruction');


fgm = imregionalmax(Iobrcbr);
%figure, imshow(fgm), title('Regional Maxima of Opening-Closing by Reconstruction');


bw = imbinarize(Iobrcbr);
%figure, imshow(bw), title('Thresholded Opening-Closing by Reconstruction');


% D = bwdist(bw);
% DL = watershed(D);
% bgm2 = DL == 0;
% figure, imshow(bgm2), title('Watershed Ridge Lines)');

%D = bwdist(fgm);
D = bwdist(fgm,'quasi-euclidean');
DL = watershed(D);
bgm2 = DL == 0;
%figure, imshow(bgm2), title('Watershed Ridge Lines)');


% Saving the images of watershed ridge lines and regional maxima detected
% if exist('filename','var') && ~exist('outfilename','var')
%     kslash = union(strfind(filename,'\'),strfind(filename,'/'));
%     filename(kslash) = '_';
%     if ~isempty(kslash)
%         filenameshort = filename(kslash(end-1)+1:end);
%     end
% end
% 
% 
% figure('visib','off')
% set(gcf,'paperpositionmode','auto')
% imshow(bgm2+fgm)
% set(gca,'Position',[0 0 1 1])
% %saveas(gcf, ['Watershed Ridge Lines and regional maxima',filesep,'Regionalmaxn and watershed_',filenameshort(1:end-4),'_',num2strDU(outidx,4),'.jpg']);
% saveas(gcf, ['NucImages\Watershed Ridge Lines and regional maxima',filesep,'Regionalmaxn and watershed_',filenameshort(1:end-4),'_',num2strDU(outidx,4),'.jpg']);
% close(gcf)
end
