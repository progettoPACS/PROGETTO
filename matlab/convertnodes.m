clc

fileToRead1='GLL32.txt';
rawData1 = importdata(fileToRead1);
[~,name] = fileparts(fileToRead1);
newData1.(genvarname(name)) = rawData1;
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

GLL32(:,2)=GLL32(:,2)/2+1/2;
GLL32(:,3)=GLL32(:,3)/2;

for i=1:size(GLL32,1)
    fprintf(1,['const Real ql32ptx',num2str(i),' = ',num2str(GLL32(i,2), '%10.22e'),', ql32ptw',num2str(i),' = ',num2str(GLL32(i,3),'%10.22e'),';\n']);
end