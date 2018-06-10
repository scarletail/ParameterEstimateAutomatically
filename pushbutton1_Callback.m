% 浏览按钮，打开txt或xls文件-----------------------------------  
function pushbutton1_Callback(hObject, eventdata, handles)  
[filename filepath fileindex]=uigetfile({'*.xls','Excel文件(*.xls)';'*.txt','文本文件(*.txt)';'*.*','所有文件(*.*)'},'选择文件');  %打开文件对话框  
if fileindex~=0   %如果没有点击取消  
    l=length(filename);   %filename包含后缀名  
    if l<=4  
        errordlg('错误文件','文件打开错误');  %错误对话框  
        return;  
    end  
    test=filename(1,l-3:l);   %文件名，截取后缀名  
    switch test  
        case '.xls'   %如果是表格文件  
            str=[filepath filename];   %拼接绝对路径  
            set(handles.edit2,'String',str);  %  
            h=waitbar(0,'正在读取文件....');   %进度条(打开文件比较慢)  
            [chengji xingming]=xlsread(str);  %根据绝对路径，打开xls文件。chengji是一个三维列向量，xingming是一个一维列向量。第一个参数chengji只保存纯数字，第二个参数xingming只保存字符串。  
            waitbar(1,h,'完成');  
            delete(h);      
            set(handles.listbox1,'String',xingming(:,1));   %设置listbox1的显示内容  
            handles.chengji=chengji;     %将chengji放入到全局变量中  
            guidata(hObject, handles);   %更新gui数据  
        case '.txt'  
            strfile=[filepath filename];   %拼接绝对路径  
            set(handles.edit2,'String',strfile);    
            fin=fopen('chengji.txt','r');   %打开txt文件  
            str=fgetl(fin);     %读取txt文件的第一行  
            %   A=fscanf(fin,'%d','HeaderLines',1);  忽略前一行的标题内容    
            %   [A B]=textscanf(fin,'%d %d');   返回一个cell类型的数据[A B]，一列一列的读    
            %   ftell();    得到读取的位置  
              
            [str1 str2 str3 str4]=strread(str,'%s %s %s %s','delimiter',' '); %以空格分割，截取每一行  
            xingming(1)=str1;  
            counter=2;  %txt文件的第一行是(name yuwen shuxue yuwen)，所以从第二行才是需要的内容。  
            while feof(fin)==0   %如果能读到txt文件的内容，(没有到txt文件的结尾)  
                str=fgetl(fin);  %读取txt文件的一行内容  
                [name yuwen shuxue yingyu]=strread(str,'%s %d %d %d','delimiter',' ');  
                xingming(counter)=name;  
                chengji(counter-1,:)=[yuwen shuxue yingyu];  
                counter=counter+1;  
            end  
            set(handles.listbox1,'String',xingming(1,:));  
            handles.chengji=chengji;  
            guidata(hObject, handles);  
            fclose(fin);  %关闭文件流  
              
        otherwise  
            errordlg('文件类型错误','文件错误');    
            return;  
    end        
end  
  
  
% listbox选项改变的回调函数-----------------------------------  
function listbox1_Callback(hObject, eventdata, handles)  
  
value=get(hObject,'Value')-1;  %value表示选择了listbox的第几项(从1开始)。减1是因为第一行是标题(name yuwen shuxue yuwen)  
if value>=1  
set(handles.edit1,'String',num2str(handles.chengji(value,:)));  
end  
  
  
% 退出按钮-------------------------------------------------  
function pushbutton2_Callback(hObject, eventdata, handles)  
  
clc;  
clear all;  
close(gcf);  