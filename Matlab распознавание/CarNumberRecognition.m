function three_numbers = CarNumberRecognition(filename, show, recognition_method)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
filt = fspecial ('unsharp');
convar = zeros (3, 1);
three_numbers = [0, 0, 0];
im = imread(filename);
init = im;
left_ = 0;
right_ = 0.99;
[count] = imhist (im);
total = sum(count);
while (sum(count(1:uint8((right_-0.05)*256)))/10 > (total - sum(count(1:uint8((right_-0.05)*256)))))
    right_ = right_ - 0.05;
end;
while (sum(count(1:uint8((left_+0.05)*256)))/10 > (total - sum(count(1:uint8((left_+0.05)*256)))))
    left_ = left_ + 0.05;
end;
im = imadjust (im, [left_, right_], []);
im = imadjust (im);
im = imfilter (im, filt);
im = medfilt2(im);
aver = mean (im(:));
[count, n] = imhist (im);
mx1 = max(count(2:250));
pos1 = 0;
 for k = 2:250
    if (count(k) == mx1)
        pos1 = k;
        break;
    end;
end;
ps = aver;
issue = double(ps);
if (abs(pos1 - issue - (180 - issue)/3) > 25)
    issue = issue + (180 - issue)/3;
else
    issue = issue - 30;
end;
issue = double(issue);
if (aver - issue > 80)
    issue = aver + 10;
end;
fdig = 0;
shft = + 10;
stop = 0;
maxdig = 0;
ind = issue;
while ((fdig < 3) && (stop < 50))
fdig = 0;
stop = stop + 1;
sz = size(im);
bw = im2bw(im, issue/256);
bw = medfilt2(bw);
bw = ~bw;
[lab, n] = bwlabel(bw, 8);
stat = regionprops(lab, 'all');
%onlydigitsearch
for k = 1:n
    if ((stat(k).BoundingBox(4) < 31) && (stat(k).BoundingBox(3) < sz(2)/6)...
            && (stat(k).BoundingBox(3) > 5) && (stat(k).BoundingBox(4) > 19))
        fdig = fdig + 1;
    end;
end;
if (fdig >= maxdig)
    maxdig = fdig;
    ind = issue;
end;
if (issue > 240) || (issue < 10)
    shft = - shft;
end;
issue = issue + shft;       
end;
if (fdig < maxdig)
    issue = ind;
    bw = im2bw(im, issue/256);
    bw = medfilt2(bw);
    bw = ~bw;
    [lab, n] = bwlabel(bw, 8);
    stat = regionprops(lab, 'all');
end;
res_issue = issue;       
vect = zeros(n, 1);
vectsz = 0;
for k = 1:n
    if ((stat(k).BoundingBox(4) < 2*sz(1)/3) && (stat(k).BoundingBox(3) < sz(2)/6)...
            && (stat(k).BoundingBox(3) > 5) && (stat(k).BoundingBox(4) > 15))
        vectsz = vectsz + 1;
        vect(vectsz) = k;
    end;
end;
vect = vect(1:vectsz);
%sort
for k = 1:vectsz - 1
    mini = k;
    for p = k+1:vectsz
        if (stat(vect(p)).BoundingBox(1) < stat(vect(mini)).BoundingBox(1))
            mini = p;
        end;
    end;
    temp = vect(k);
    vect(k) = vect(mini);
    vect(mini) = temp;
end;
st = stat(1:vectsz);
st(:) = stat(vect(:));
correctness = zeros(vectsz, 1);
corr = zeros(vectsz, 2);
for k = 1:vectsz
    current = st(k).BoundingBox;
    if (current(4) < 18) 
        correctness(k) = correctness(k) - 10; 
    end;        
    if (current(3) > current(4))
        correctness(k) = correctness (k) - 20;
        continue;%definetly wrong
    end;
    for j = k+1:vectsz
        sec = st(j).BoundingBox;
        if ((abs(current(4) - sec(4)) > 4) || (abs(current(2) - sec(2)) > 5))
            continue;            
        end;
        if (sec(1) - current(1) > 40)%3*(current(3) + sec(3))
            correctness(k) = correctness(k) - 40;%not found next digit
            break;    
        end;
        if (sec(1) - current(1) > 5)
            correctness(k) = correctness (k) - ((sec(1) - current(1)) - 5);
        end;
        %found second digit
        corr(k, 1) = j;
        correctness(k) = correctness (k) - abs(current(4) - sec(4))/2 + 10 + sec(4);
        correctness(k) = correctness (k) - abs(current(3) - sec(3));
        break;
    end;
    if (corr(k,1) == 0)
        correctness(k) = correctness (k) - 20;
        continue;
    end;
    current = st(corr(k, 1)).BoundingBox;
    for i = corr(k,1) + 1:vectsz
        sec = st(i).BoundingBox;
        if ((abs(current(4) - sec(4)) > 4) || (abs(current(2) - sec(2)) > 5))
            continue;            
        end;
        if (sec(1) - current(1) > 2*(current(3) + sec(3)))
            %not found next digit
            break;    
        end;
        %found third digit
        corr(k, 2) = i;
        correctness(k) = correctness (k) - abs(current(4) - sec(4))/2 + sec(4);
        correctness(k) = correctness (k) - abs(current(3) - sec(3))/2;
        break;
    end;
    if (corr(k,2) == 0)
        correctness(k) = correctness (k) - 20;
        corr(k, 2) = corr(k, 1);
    end; 
    correctness(k) = correctness(k) + 30;
    %check for letters;
    %leftchecking
    current = st(k).BoundingBox;
    for j = k-1:-1:1
        left = st(j).BoundingBox;
        if ((current(2) > left(2) + left(4)) || (left(2) > current(2) + current(4)))
            continue;
        end;
        if  (current(1) - left(1) > 2*(current(3) + left(3)))
            break;
        end;
        if ((current(4) < 3*left(4)) && (current(4) > left(4)))
            correctness(k) = correctness(k) + 5;
            break;
        end;
    end;
    %rightchecking
    current = st(corr(k, 2)).BoundingBox;
    find = 0;
    for j = corr(k, 2)+1:vectsz
        left = st(j).BoundingBox;
        if ((current(2) > left(2) + left(4)) || (left(2) > current(2) + current(4)))
            continue;
        end;
        if  (current(1) - left(1) > 2*(current(3) + left(3)))
            break;
        end;
        if ((current(4) < 3*left(4)) && (current(4) > left(4)))
            correctness(k) = correctness(k) + 5;
            if (find == 0) find = 1;
            else break;
            end;           
        end;
    end;
end;
maxc = max(correctness);
for k = 1:vectsz
    if (maxc == correctness(k))
        d1 = k;
        break;
    end;
end;
d2 = corr(d1, 1);
d3 = corr(d1, 2);
if (d2 + d3 == 0)
    szd = 1;
    d = [d1];
elseif (d3 == d2)
    szd = 2;
    d = [d1, d2];
else
    d = [d1, d2, d3];
    szd = 3;
end;
%find third digit
if (szd == 2)
    stf = uint16(st(d(1)).BoundingBox);
    sts = uint16(st(d(2)).BoundingBox);
    if (sts(1) - stf(1) < stf (3) + min(stf(3), sts(3)))
        %try left or right
        %left
        where = 0;
        if (d1 == 1) || (stf(1) - stf(3) > st(d1-1).BoundingBox(1) + st(d1-1).BoundingBox(3))
            where = 1;
            search_zone = [stf(2) - 1, stf(1) - stf(3) - 3, ...
                  stf(2) + stf(4) + 1, stf(1) - 1];
        end;
        %right
        if (where == 0)
            if (d2 == vectsz) || (sts(1) + sts(3) < st(d2+1).BoundingBox(1))
                where = 3;
                search_zone = [sts(2) - 1, sts(1) + sts(3) + 1, ...
                      sts(2) + sts(4) + 1 , sts(1) + 2*sts(3) + 3];
            end;
        end;
    else
        %try between
        search_zone = [min(stf(2), sts(2)) - 2, stf(1) + stf(3) + 1, ...
                  max(stf(2) + stf(4), sts(2) + sts(4)) + 2, sts(1) - 1];
        %
        where = 2;
    end;
    if (where == 0)
        std = cat(1, st(d(1)), st(d(2)));
    else
    ims = im(search_zone(1):search_zone(3), search_zone(2):search_zone(4));
    onedig = 0;
    while (onedig ~= 1)
        bwn = im2bw(ims, issue/256);
        bwn = ~bwn;
        bwn = medfilt2(bwn);
        [lab, onedig] = bwlabel(bwn);
        issue = issue + 5;
        if (issue > 250)
            break;
        end;
    end;
    statt = regionprops(lab, 'all');
    st2 = statt(1);
    st2(1).BoundingBox = uint16(st2(1).BoundingBox);
    st2(1).BoundingBox(1) = st2(1).BoundingBox(1) + search_zone(2) - 1;
    st2(1).BoundingBox(2) = st2(1).BoundingBox(2) + search_zone(1) - 1;
    if (where == 2)
        std = cat(1, st(d(1)), st2, st(d(2)));
    elseif (where == 1)
        std = cat(1, st2, st(d(1)), st(d(2)));
    elseif (where == 3)
        std = cat(1, st(d(1)), st(d(2)), st2);
    end;
    szd = 3;
    end;
elseif (szd == 3)
    std = cat(1, st(d(1)), st(d(2)), st(d(3)));
else
    std = st(d(1));
end;
%add broken parts
mn = zeros(szd, 1);
for j = 1:szd
    mn(j) = std(j).BoundingBox(3);
end;
if (min(mn) < 8)
for j = 1:n
    cur = stat(j).BoundingBox;
    if ((cur(4) > 30) || (cur(3) > 8) || (cur(1) > std(szd).BoundingBox(1) + std(szd).BoundingBox(3) + 3) ...
            || (cur(1) + cur(3) < std(1).BoundingBox(1) - 1) ...
            || (cur(2) < std(1).BoundingBox(2) - 2) || (cur(2) + cur(4) >  std(1).BoundingBox(2) + std(1).BoundingBox(4) + 2))
        continue;
    elseif (szd > 1)
        if (cur(2) == std(2).BoundingBox(2))
            continue;
        end;
    end;
        for k = 1:szd
            f = k;
            if (cur(1) < std(k).BoundingBox(1));
                break;
            end;
        end;
        if (f == 1)
            side = 0;            
        elseif (f == szd) && (cur(1) > std(szd).BoundingBox(1))
            side = 1;
        elseif (std(f).BoundingBox(3) < std(f - 1).BoundingBox(3))
            side = 0;
        else
            side = 1;
            f = f - 1;
        end;
        sizes = std(f).BoundingBox;
        if (side == 1)
            std(f).BoundingBox(3) = cur(1) + cur(3) - sizes(1);
        else
            std(f).BoundingBox(3) = std(f).BoundingBox(3) + (sizes(1) - cur(1));
            std(f).BoundingBox(1) = cur(1);
        end;
        sizes = std(f).BoundingBox;
        sizes(2) = min(std(f).BoundingBox(2), cur(2));
        sizes(4) = max(std(f).BoundingBox(2) + std(f).BoundingBox(4), cur(2)+cur(4)) - sizes(2);
        std(f).BoundingBox = sizes;
end;
end;
for k = 1:szd-1
    if (std(k+1).BoundingBox(1) - (std(k).BoundingBox(1) + std(k).BoundingBox(3)) > 5)
        std(k+1).BoundingBox(1) = std(k+1).BoundingBox(1) - 4;
        std(k+1).BoundingBox(3) = std(k+1).BoundingBox(3) + 4;
    end;
end;

if (std(1).BoundingBox(3) < 10) && ((d(1) == 1) || (std(1).BoundingBox(1) > st(d(1) - 1).BoundingBox(1) + st(d(1) - 1).BoundingBox(1) + 5))
    std(1).BoundingBox(1) = std(1).BoundingBox(1) - 4;
    std(1).BoundingBox(3) = std(1).BoundingBox(3) + 4;
end;

if (show)
    imrgb = cat(3, init, init, init);
    for t = 1:szd
        cur = uint16(std(t).BoundingBox);
        imrgb(cur(2):cur(2)+cur(4), cur(1)-1, 1) = 255;
        imrgb(cur(2):cur(2)+cur(4), cur(1)-1, 2) = 0;
        imrgb(cur(2):cur(2)+cur(4), cur(1)-1, 3) = 0;
        imrgb(cur(2):cur(2)+cur(4), cur(1)+cur(3)+1, 1) = 255;
        imrgb(cur(2):cur(2)+cur(4), cur(1)+cur(3)+1, 2) = 0;
        imrgb(cur(2):cur(2)+cur(4), cur(1)+cur(3)+1, 3) = 0;
        imrgb(cur(2)-1, cur(1):cur(1)+cur(3), 1) = 255;
        imrgb(cur(2)-1, cur(1):cur(1)+cur(3), 2) = 0;
        imrgb(cur(2)-1, cur(1):cur(1)+cur(3), 3) = 0;
        imrgb(cur(2)+cur(4)+1, cur(1):cur(1)+cur(3), 1) = 255;
        imrgb(cur(2)+cur(4)+1, cur(1):cur(1)+cur(3), 2) = 0;
        imrgb(cur(2)+cur(4)+1, cur(1):cur(1)+cur(3), 3) = 0;
    end;
    figure, imshow(imrgb);
end;
%0 recogition method (invariants)
if (recognition_method == 0)
    for i = 1:szd
        onedig = 0;
        stc = std(i);
        issue = res_issue - 5;
        sizes = uint16(stc.BoundingBox);
        ims = im(sizes(2):sizes(2)+sizes(4), sizes(1):sizes(1)+sizes(3));
        shft = 5;
        while (onedig ~= 1)
            bwn = im2bw(ims, issue/256);
            bwn = ~bwn;
            bwn = medfilt2(bwn);
            [lab, onedig] = bwlabel(bwn);
            issue = issue + shft;
            if (issue > 250)
                shft = -5;
            end;
            if (issue < 10)
                break;
            end;
        end;
        stt = regionprops(lab, 'all');
        stc = stt(1);
        sizes = stc.BoundingBox;
        sq = sizes(3)*sizes(4);
        convar(i) = double((sq - stc.ConvexArea))/sq;
        %try connect smth
        EL = stc.EulerNumber;
        corrected = false;
        dif = 0;
        if ((EL == 1) || ((EL == 0) && (convar(i) < 0.2)))
            bn = stc.Image;
            se = strel('diamond', 1);
            for t = 1:2
                bn = imdilate(bn, se);
                [lb, ca] = bwlabel(bn, 4);
                s = regionprops(lb, 'all');
                dif = s(1).FilledArea - s(1).Area;
                if (EL > s(1).EulerNumber)
                    if (dif > 5 + 2*EL)
                        if (EL == 1) && (dif > 30)
                            EL = 0;
                            corrected = true;
                            break;
                        else
                        EL = s(1).EulerNumber;
                        corrected = true;
                        break;
                        end;
                    else
                        continue;
                    end;
                elseif (s(1).EulerNumber > EL)
                    break;
                end; 
                if (EL == -1)
                    break;
                end;
            end;              
        end;
        %endtry
        cm = stc.Centroid(2);
        sizes = stc.BoundingBox;
        c = sizes (2) + sizes(4)/2;
        eps = 1;
        if (EL == 1)
            bnr = stc.Image(:,uint16(sizes(3)/2)+1:end);
            [~, ca] = bwlabel(bnr);
            smthru = sum(sum(bnr(uint16(sizes(4)/4):uint16(sizes(4)/3), :)));
            bnl = stc.Image(:,1:uint16(sizes(3)/2));
            smthlb = sum(sum(bnl(uint16(2*sizes(4)/3):uint16(3*sizes(4)/4), :)));
            smthlu = sum(sum(bnl(uint16(sizes(4)/4):uint16(sizes(4)/3), :)));
            if (smthru + smthlb < 5)
                three_numbers(i) = 5;
            elseif (convar(i) < 0.2) || (~((stc.Extrema(7,2) > sizes(2) + sizes(4)/2) && (stc.Extrema(7,2) < sizes(2) + 3*sizes(4)/4 + 2)) && (cm >= c))%2, 3, 5
                if (ca == 1) && (smthlu + smthlb  < 10)
                    three_numbers(i) = 3;
                elseif(smthru > 5)
                    three_numbers(i) = 2;
                elseif (cm > c + 1) 
                    three_numbers(i) = 6;
                else
                    three_numbers(i) = 5;
                end;
            else %1, 4, 7
                if ((stc.Extrema(7,2) > sizes(2) + sizes(4)/2) && (stc.Extrema(7,2) < sizes(2) + 3*sizes(4)/4 + 2)) %left-bottom
                    three_numbers(i) = 4;
                elseif (stc.Extrema(5, 1) < sizes(1) + 3*sizes(3)/4)%bottom-right
                    three_numbers(i) = 7;
                elseif (cm > c + 1) && (stc.Extrema(2, 1) < sizes(1) + sizes(3)/2)
                    three_numbers(i) = 6;
                else
                    three_numbers(i) = 1;
                end;
            end;
        elseif (EL == 0) %6, 9, 0
            if (convar(i) < 0.2) && (4*stc.FilledArea < 5*stc.Area) && (~corrected)
                three_numbers(i) = 8;
            elseif(dif > 30) && corrected && (convar(i) < 0.15)
                three_numbers(i) = 0;
            elseif (cm > c + eps)
                three_numbers(i) = 6;
            elseif (cm < c - eps)
                three_numbers(i) = 9;
            else
                three_numbers(i) = 0;
            end;
        elseif (EL == -1)
            three_numbers(i) = 8;
        else
            %'warning EulerNumber' 
            three_numbers(i) = 0;
            %just random
        end;
    end;
%1 recognition method pattern matching
elseif(recognition_method == 1)
    %create patterns
    ptrnt = uint16(zeros(25, 15, 10));
    %buildpatterns();
    ptrnt = int16(ptrnt);
    for j = 1:10
        ptrnt(:, :, j) = imread(strcat('patterns/', int2str(j-1), '.bmp'));
        ptrnt (:,:,j) = ptrnt(:,:,j) - mean(mean(ptrnt(:,:,j)));
    end;
    for k = 1:szd
        sizes = uint16(std(k).BoundingBox);
        curr = imresize(im(sizes(2):sizes(2) + sizes(4), sizes(1):sizes(1) + sizes(3)), [25, 15]);
        curr = imadjust(curr);
        curr = int16(curr);
        curr = curr - mean(curr(:));
        psnr = zeros(10, 1);
        for j = 1:10
            psnr(j) = sum(sum((curr - ptrnt(:,:,j)).^2));
        end;
        min_psnr = min(psnr);
        for j = 1:10
            if (psnr(j) == min_psnr)
                three_numbers(k) = j - 1;
                break;
            end;
        end;
    end;
end;
end

