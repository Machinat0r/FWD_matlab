clear;clc
fid=fopen('CSES_SCM_0123950_20200427_155704_20200427_163356_L0_0000044976.dat','rb');
data0=fread(fid,10000);%��ȡ����
array=dec2hex(data0);%��ʮ����ת��Ϊʮ������
array1=array(:,1);
array2=array(:,2);
data=[array1';array2'];
data=data(:)';%�����о��󽻲�ϲ�Ϊһ����ת��
load_location=strfind(data,'146F0102');%�����ַ����������غɴ����λ��
disp(load_location);

size1=load_location+20;
size2=load_location+23;
n=numel(load_location);

workmode1=[];workmode2=[];workmode3=[];
for i = 1:n
    size0=4*hex2dec(data(size1(i):size2(i)));
    switch size0
        case 3744
            workmode1(end+1) = load_location(i);
        case 20128
            workmode2(end+1) = load_location(i);
        case 20144
            workmode3(end+1) = load_location(i);
    end
end
j = length(workmode1)
k = length(workmode2)
l = length(workmode3)
   
time_sec1=load_location+24;%���ݰ�ʱ���������ݳ�ʼλ��
time_sec2=load_location+31;%���ݰ�ʱ����������ĩλ��
time_millisec1=load_location+32;%���ݰ�ʱ����������ݳ�ʼλ��
time_millisec2=load_location+35;%���ݰ�ʱ�����������ĩλ��
k=datenum(2009,01,01,00,00,00);

for i=1:n
    time_sec_dec=hex2dec(data(time_sec1(i):time_sec2(i)));%����ת��Ϊʮ����
    time_millisec_dec=hex2dec(data(time_millisec1(i):time_millisec2(i)));%������ת��Ϊʮ����
    time_sec(i)=time_sec_dec;
    time_millisec(i)=time_millisec_dec;
    time_day=time_sec_dec/86400+time_millisec_dec/86400000+k;
    SCM_time(i)=string(datestr(time_day,'yyyy-mm-dd HH-MM-SS'));
end  
 
mode1=load_location+36;%���ݰ�����ģʽ���ݳ�ʼλ��
mode2=load_location+39;%���ݰ�����ģʽ����ĩλ��
for i=1:n
    mode_hex(i)=hex2dec(data(mode1(i):mode2(i)));%��ȡ���ݰ�����ģʽ����
end

satpos_x1=load_location+68;
satpos_x2=load_location+75;
satpos_y1=load_location+76;
satpos_y2=load_location+83;
satpos_z1=load_location+84;
satpos_z2=load_location+91;
for i=1:n
    satpos_x(i)=hex2dec(data(satpos_x1(i):satpos_x2(i)));%����λ��X
    satpos_y(i)=hex2dec(data(satpos_y1(i):satpos_y2(i)));%����λ��Y
    satpos_z(i)=hex2dec(data(satpos_z1(i):satpos_z2(i)));%����λ��Z
end

GeoLonLat1=load_location+92;
GeoLonLat2=load_location+99;
GeoLonLat3=load_location+100;
GeoLonLat4=load_location+107;

for i=1:n
    GeoLonLat_hex1=hex2dec(data(GeoLonLat1(i):GeoLonLat2(i)));
    if GeoLonLat_hex1/1e7>180
        Longitude(i)=(GeoLonLat_hex1-16^8)/1e7;
    else
        Longitude(i)=GeoLonLat_hex1/1e7;
    end
    
    GeoLonLat_hex2=hex2dec(data(GeoLonLat3(i):GeoLonLat4(i)));
    if GeoLonLat_hex2/1e7>90
        Latitude(i)=(GeoLonLat_hex2-16^8)/1e7;
    else
        Latitude(i)=GeoLonLat_hex2/1e7;
    end
end



















