clc
%Sito dove prendere i nodi
%http://keisan.casio.com/has10/SpecExec.cgi?id=system/2006/1280801905
fileToRead1='GLL64.txt';
rawData1 = importdata(fileToRead1);
[~,name] = fileparts(fileToRead1);
newData1.(genvarname(name)) = rawData1;
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

GLL64(:,2)=GLL64(:,2)/2+1/2;
GLL64(:,3)=GLL64(:,3)/2;

for i=1:size(GLL64,1)
    fprintf(1,['const Real ql64ptx',num2str(i),' = ',num2str(GLL64(i,2), '%10.22e'),', ql64ptw',num2str(i),' = ',num2str(GLL64(i,3),'%10.22e'),';\n']);
end

for i=1:64
    fprintf(1,['QuadraturePoint ( ql64ptx',num2str(i),', ql64ptw',num2str(i),'),\n']);
end