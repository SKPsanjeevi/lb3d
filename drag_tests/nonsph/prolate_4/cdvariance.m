% deg = [0];
% deg = [0 24 35 45 55 66 90];
% deg = [0 10 20 30 40 50 60 70 80 90];
deg = [0 10 30 45 60 80 90];
% deg = [90];
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
R = 20.3581321245;
%R = 27.1441761659;
% v = 0.02;
% R = 6.5;
string='Re300ang';
%string='Re90';
%delete(string);
result=zeros(size(deg,2),3);
sd=zeros(size(deg,2),2);
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
    data=textscan(fid,'%*d%*f%f%f%f%*f%*f%*f%f%f%*f%f%f');
%    data=textscan(fid,'%*d%*f%f%f%f%*f%*f');
    fclose(fid);
    num=cell2mat(data);
    last=size(num,1);
    begin=last-1000;
    
    cl=mean(num(begin:last,1))/(0.5*pi*R^2*v^2);
    CD=mean(num(begin:last,2))/(0.5*pi*R^2*v^2);
    ct=mean(num(begin:last,3))/(0.5*pi*R^3*v^2);

%    clp=mean(num(begin:last,4))/(0.5*pi*R^2*v^2);
%    cdp=mean(num(begin:last,5))/(0.5*pi*R^2*v^2);
%
%    clv=mean(num(begin:last,6))/(0.5*pi*R^2*v^2);
%    cdv=mean(num(begin:last,7))/(0.5*pi*R^2*v^2);

    CD_t=num(begin:last,2)/(0.5*pi*R^2*v^2);
    CD_sd=sqrt(mean((CD_t-CD).^2));
    CD_rel_sd=CD_sd/CD*100;


    result(i,:)=[CD cl ct];
    sd(i,:)=[CD_sd CD_rel_sd];
%    result(i,:)=[CD cl ct clp cdp clv cdv];

    cd ../../
end
string_coeffs=strcat(string,'_coeffs');
string_sd=strcat(string,'_sd');
dlmwrite(string_coeffs, result,'delimiter','\t');
dlmwrite(string_sd, sd,'delimiter','\t');
