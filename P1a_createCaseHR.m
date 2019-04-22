clc;
clear variables;
close all;

folder='Data/IN/';
files=[];
files=dir([folder '*.csv']);
startcat=5;
endcat=5;
subratio=1;
casetitle=['_single_day_1_per_1_minute_in'];
i=1;
HF=[];
for i=startcat:endcat
file=files(i).name;
filepath=[folder '/' file];
fid1 = fopen(filepath, 'r');
Header = fgetl(fid1);
fclose(fid1);

Header = regexp(Header, '([^,]*)', 'tokens');
Header = cat(2, Header{:});
fid1 = fopen(filepath, 'r');
Data=[];

while ~feof(fid1)
 line=fgetl(fid1);
 Cline  = textscan(line, '%s%s%f', 1,'Delimiter', ',');
 Data=[Data;Cline];
end


HeartRateDataFrame=Data(2:end,3);
H=cell2mat(HeartRateDataFrame);
H(isnan(H))=[];
HF=[HF ;H];
end

%% Subsample

HF=HF(1:subratio:end);
HF=HF(1:(floor(sqrt(length(HF))))^2);


divider=optimalDivider(HF);
Hr=reshape(HF,divider,[]);
figure(2)
plot(HF);



figure(4)
[U,S,V]=svd(Hr);
SD=diag(S);
plot(SD/sum(SD));
grid on

figure(5)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8.3 2.5])
plot(Hr(:),'LineWidth',2)
set(findall(gcf,'-property','FontSize'),'FontSize',10)
grid on
print(['Figures/hr_plot.jpg'],'-djpeg','-r100')


csvwrite(['Data/' casetitle '.csv'],Hr);
