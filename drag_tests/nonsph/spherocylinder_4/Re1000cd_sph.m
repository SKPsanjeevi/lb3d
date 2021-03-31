close all;
deg = [0 10 30 45 60 80 90];
v = 0.09;
R = 26.4776125149;
string='Re2000ang';
result=zeros(size(deg,2),3);

ang={'0','10','30','45','60','80','90'};


hFig=figure('units','inches','position',[4 4 6 5]);
colormap(hFig, jet);
set(groot,'DefaultAxesColorOrder',[0 0 1; 0 .5 0; 1 0 0; 0 .75 .75; ...
                                   .75 0 .75; .75 .75 0; .25 .25 .25])
hold on;
for i=1:size(deg,2)
    dispchar=strcat('$',char(ang(i)),'^\circ$');
    target=strcat(string,num2str(deg(i)),'/Production');
    cd(target);
    delete('*.xdr','*.h5');
    fname=dir('*.asc');
    if(size(fname,1)>1)
        fname=dir('md-cfg_out_p00000001*.asc');
    end
    fid=fopen(fname.name);
    data=textscan(fid,'%*d%*f%f%f%f%*f%*f');
    fclose(fid);
    
    fid=fopen(fname.name);
    time=textscan(fid,'%d%*f%*f%*f%*f%*f%*f');
    fclose(fid);
    
    num=cell2mat(data);
    tdata=cell2mat(time);
    cl=num(:,1)/(0.5*pi*R^2*v^2);
    CD=num(:,2)/(0.5*pi*R^2*v^2);
    ct=num(:,3)/(0.5*pi*R^3*v^2);
    t = tdata(:,1);
    tds=double(t)*v/(2*R);
    plot(tds,CD,'DisplayName',dispchar);
    cd ../../
end

xlim([0 85]);
ylim([0.1 1.2]);
xlabel('$t|\textbf{u}_\infty|/d_p$');
ylabel('$C_D$');
