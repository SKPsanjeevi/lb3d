% deg = [0];
% deg = [0 24 35 45 55 66 90];
% deg = [0 10 20 30 40 50 60 70 80 90];
deg = [0 10 30 45 60 80 90];
% deg = [''];
% deg = '10Dds_1_4th';
% v = 0.043679;
% R = 20.604856366;
v = 0.09;
%v = 0.02;
%v = 0.043679;
%v = 0.069886;
%R = 11.4471424255;
%R = 10.8576704664;
%R = 13.572088083;
R = 26.4776125149;
% v = 0.02;
% R = 6.5;
string='Re1000ang';
%string='Re90';
%delete(string);
result=zeros(size(deg,2),3);
for i=1:size(deg,2)
%for i=1:1%size(deg,1)
    target=strcat(string,num2str(deg(i)),'/Production');
    %target=strcat(string,deg(i,:),'/Production');    
    %target=strcat(string,'/Production');    
    cd(target);
    delete('*.xdr','*.h5');
    fname=dir('*.asc');
    if(size(fname,1)>1)
        fname=dir('md-cfg_out_p00000001*.asc');
    end
    fid=fopen(fname.name);
    %data=textscan(fid,'%*d%*f%*f%*f%*f%*f%*f%*f%f%f%*f%*f%*f%f%*f%*f');
    data=textscan(fid,'%*d%*f%f%f%f%*f%*f');
    fclose(fid);
    num=cell2mat(data);
    last=size(num,1);
    begin=last-1500;
    cl=mean(num(begin:last,1))/(0.5*pi*R^2*v^2);
    CD=mean(num(begin:last,2))/(0.5*pi*R^2*v^2);
    ct=mean(num(begin:last,3))/(0.5*pi*R^3*v^2);
    result(i,:)=[CD cl ct];
    cd ../../
end
string=strcat(string,'_coeffs');
dlmwrite(string, result,'delimiter','\t');
