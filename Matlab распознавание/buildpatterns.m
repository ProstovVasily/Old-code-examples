function buildpatterns()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
ptrnt = uint16(zeros(25, 15, 10));
for j = 1:10
    for k = 1:100
        ptrnt(:,:,j) = ptrnt(:,:,j) + uint16(imresize(imread(strcat('digits/', int2str(j- 1), '/', int2str(k), '.bmp')), [25, 15]));
    end;
    ptrnt(:,:,j) = ptrnt(:,:,j)/100;
    imwrite (uint8(ptrnt(:,:,j)), strcat('patterns/', int2str(j - 1), '.bmp'), 'bmp');
end;

end

