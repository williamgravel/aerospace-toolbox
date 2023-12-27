function [sv] = selectSVCode(PRN)
arguments
    PRN (:,1) double
end

% Build table
data = table('VariableNames',{'PS 1','PS 2','Code Delay','First 10 Chips'},...
    'VariableTypes',{'double','double','double','string'},...
    'RowNames',string((1:37)'),'Size',[37,4]);

% Populate table
data('1',:) = {2, 6, 5, oct2hex('1440')};
data('2',:) = {3, 7, 6, oct2hex('1620')};
data('3',:) = {4, 8, 7, oct2hex('1710')};
data('4',:) = {5, 9, 8, oct2hex('1744')};
data('5',:) = {1, 9, 17, oct2hex('1133')};
data('6',:) = {2, 10, 18, oct2hex('1455')};
data('7',:) = {1, 8, 139, oct2hex('1131')};
data('8',:) = {2, 9, 140, oct2hex('1454')};
data('9',:) = {3, 10, 141, oct2hex('1626')};
data('10',:) = {2, 3, 251, oct2hex('1504')};
data('11',:) = {3, 4, 252, oct2hex('1642')};
data('12',:) = {5, 6, 254, oct2hex('1750')};
data('13',:) = {6, 7, 255, oct2hex('1764')};
data('14',:) = {7, 8, 256, oct2hex('1772')};
data('15',:) = {8, 9, 257, oct2hex('1775')};
data('16',:) = {9, 10, 258, oct2hex('1776')};
data('17',:) = {1, 4, 469, oct2hex('1156')};
data('18',:) = {2, 5, 470, oct2hex('1467')};
data('19',:) = {3, 6, 471, oct2hex('1633')};
data('20',:) = {4, 7, 472, oct2hex('1715')};
data('21',:) = {5, 8, 473, oct2hex('1746')};
data('22',:) = {6, 9, 474, oct2hex('1763')};
data('23',:) = {1, 3, 509, oct2hex('1063')};
data('24',:) = {4, 6, 512, oct2hex('1706')};
data('25',:) = {5, 7, 513, oct2hex('1743')};
data('26',:) = {6, 8, 514, oct2hex('1761')};
data('27',:) = {7, 9, 515, oct2hex('1770')};
data('28',:) = {8, 10, 516, oct2hex('1774')};
data('29',:) = {1, 6, 859, oct2hex('1127')};
data('30',:) = {2, 7, 860, oct2hex('1453')};
data('31',:) = {3, 8, 861, oct2hex('1625')};
data('32',:) = {4, 9, 862, oct2hex('1712')};
data('33',:) = {5, 10, 863, oct2hex('1745')};
data('34',:) = {4, 10, 950, oct2hex('1713')};
data('35',:) = {1, 7, 947, oct2hex('1134')};
data('36',:) = {2, 8, 948, oct2hex('1456')};
data('37',:) = {4, 10, 950, oct2hex('1713')};

% Select entries
sv = data(string(PRN),:);

end
