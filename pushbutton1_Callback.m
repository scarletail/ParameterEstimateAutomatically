% �����ť����txt��xls�ļ�-----------------------------------  
function pushbutton1_Callback(hObject, eventdata, handles)  
[filename filepath fileindex]=uigetfile({'*.xls','Excel�ļ�(*.xls)';'*.txt','�ı��ļ�(*.txt)';'*.*','�����ļ�(*.*)'},'ѡ���ļ�');  %���ļ��Ի���  
if fileindex~=0   %���û�е��ȡ��  
    l=length(filename);   %filename������׺��  
    if l<=4  
        errordlg('�����ļ�','�ļ��򿪴���');  %����Ի���  
        return;  
    end  
    test=filename(1,l-3:l);   %�ļ�������ȡ��׺��  
    switch test  
        case '.xls'   %����Ǳ���ļ�  
            str=[filepath filename];   %ƴ�Ӿ���·��  
            set(handles.edit2,'String',str);  %  
            h=waitbar(0,'���ڶ�ȡ�ļ�....');   %������(���ļ��Ƚ���)  
            [chengji xingming]=xlsread(str);  %���ݾ���·������xls�ļ���chengji��һ����ά��������xingming��һ��һά����������һ������chengjiֻ���洿���֣��ڶ�������xingmingֻ�����ַ�����  
            waitbar(1,h,'���');  
            delete(h);      
            set(handles.listbox1,'String',xingming(:,1));   %����listbox1����ʾ����  
            handles.chengji=chengji;     %��chengji���뵽ȫ�ֱ�����  
            guidata(hObject, handles);   %����gui����  
        case '.txt'  
            strfile=[filepath filename];   %ƴ�Ӿ���·��  
            set(handles.edit2,'String',strfile);    
            fin=fopen('chengji.txt','r');   %��txt�ļ�  
            str=fgetl(fin);     %��ȡtxt�ļ��ĵ�һ��  
            %   A=fscanf(fin,'%d','HeaderLines',1);  ����ǰһ�еı�������    
            %   [A B]=textscanf(fin,'%d %d');   ����һ��cell���͵�����[A B]��һ��һ�еĶ�    
            %   ftell();    �õ���ȡ��λ��  
              
            [str1 str2 str3 str4]=strread(str,'%s %s %s %s','delimiter',' '); %�Կո�ָ��ȡÿһ��  
            xingming(1)=str1;  
            counter=2;  %txt�ļ��ĵ�һ����(name yuwen shuxue yuwen)�����Դӵڶ��в�����Ҫ�����ݡ�  
            while feof(fin)==0   %����ܶ���txt�ļ������ݣ�(û�е�txt�ļ��Ľ�β)  
                str=fgetl(fin);  %��ȡtxt�ļ���һ������  
                [name yuwen shuxue yingyu]=strread(str,'%s %d %d %d','delimiter',' ');  
                xingming(counter)=name;  
                chengji(counter-1,:)=[yuwen shuxue yingyu];  
                counter=counter+1;  
            end  
            set(handles.listbox1,'String',xingming(1,:));  
            handles.chengji=chengji;  
            guidata(hObject, handles);  
            fclose(fin);  %�ر��ļ���  
              
        otherwise  
            errordlg('�ļ����ʹ���','�ļ�����');    
            return;  
    end        
end  
  
  
% listboxѡ��ı�Ļص�����-----------------------------------  
function listbox1_Callback(hObject, eventdata, handles)  
  
value=get(hObject,'Value')-1;  %value��ʾѡ����listbox�ĵڼ���(��1��ʼ)����1����Ϊ��һ���Ǳ���(name yuwen shuxue yuwen)  
if value>=1  
set(handles.edit1,'String',num2str(handles.chengji(value,:)));  
end  
  
  
% �˳���ť-------------------------------------------------  
function pushbutton2_Callback(hObject, eventdata, handles)  
  
clc;  
clear all;  
close(gcf);  