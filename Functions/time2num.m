function num = time2num( in )
num=0;
str=in{1};
C=textscan(str, '%f%f', 1,'Delimiter', ':');
num=C{1}*60+C{2};
end

