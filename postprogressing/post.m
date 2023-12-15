% data = readtable("data");
% x= data.Var2
% y= data.Var3
% trapz(x,y)
clear all
data.p = readtable("cell_P.csv");
data.v = readtable("cell_V.csv");
data.rho = readtable("cell_rho.csv");
time = data.p.Var1
yp= data.p(:,2:size(data.p,2))
yv= data.v(:,2:size(data.v,2))
yrho= data.rho(:,2:size(data.rho,2))
yp = table2array(yp)
yv = table2array(yv)
yrho = table2array(yrho)

x=1:size(yv,2)
% 
% data = readtable("LF.csv");
% y= data(:,2:size(data,2))
% y2 = table2array(y)
% 
% data = readtable("upwind.csv");
% y= data(:,2:size(data,2))
% y3 = table2array(y)
% 
% data = readtable("Lax-Wendroff.csv");
% y= data(:,2:size(data,2))
% y4 = table2array(y)





% 创建 VideoWriter 对象
outputVideo = VideoWriter('data.mp4', 'MPEG-4');
outputVideo.FrameRate = 100;
open(outputVideo);
% 遍历时间步长
for t = 1:size(yv,1)
    % 对网格上的点进行插值
    plot(x,yp(t,:),x,yv(t,:),x,yrho(t,:))
%     plot(x,y1(t,:),x,y2(t,:),x,y3(t,:),x,y4(t,:))
%     legend('FTCS','LF','upwind','Lax')   
%     xmin=0
%     ymin=min(min(y1))
%     xmax=size(y,2)
%     ymax =max(max(y1))
%     axis([xmin xmax ymin ymax]); 
    legend('pressure','velocity','density')
    % 更新图形
    drawnow;
    % 写入当前图形帧到视频文件
    writeVideo(outputVideo, getframe(gcf));
end

close(outputVideo);

