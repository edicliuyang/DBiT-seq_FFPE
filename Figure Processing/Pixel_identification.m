cd 'C:\Users\EDIC-THINKPAD\Desktop\FFPEs\lymph'
I = imread('lymph.jpg');
I = rgb2gray(I);
BW = imbinarize(I);

figure
imshow(BW)

pixel = 50;
pixel_count = 2*pixel-1;
[numRows,numCols] = size(BW);
pixel_w = numCols/pixel_count;
pixel_h = numRows/pixel_count;
str ="";

for i = 1:50
    y = round(2*(i-1)*pixel_h + 1);
    for j = 1:50
        x = round(2*(j-1)*pixel_w + 1);
        pixel = BW(y:round(y+pixel_h-1),x:round(x+pixel_w-1));
        C = sum(pixel,'all');
        if C > 0
            str = str+','+j+'x'+i;   
        end
    end    
end
fid = fopen('position.txt','wt');
fprintf(fid, str);
fclose(fid);