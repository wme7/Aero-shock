%function Modify_Grid()
% function [TotNum_DM,PolyDegN_max,PolyDegN_DM,DM_Vertex,DM_Connect,TotNum_EG,...
%           EG_glo_Connect,TotNum_in_EG,in1_EG,in2_EG,EG_in_Connect,TotNum_out_EG,...
%           out1_EG,out2_EG,out3_EG,out4_EG,EG_out_Connect] = Plot_Grid()
close all
clear all
BD_type=0;
NumData_in = input('Number of domain =','s') ; 

NumDM = sprintf('%sdom.in',NumData_in) ;
NumEG = sprintf('%sedg.in',NumData_in) ;

[fid,message] = fopen ( NumDM , 'r' ) ;

head{1} = fgetl(fid) ; 
[TotNum_DM] = fscanf ( fid , '%f' , 1 ) ;
head{1} = fgetl(fid) ;
[PolyDegN_max(1)] = fscanf ( fid , '%f' , 1 ) ;
head{1} = fgetl(fid) ;
[PolyDegN_max(2)] = fscanf ( fid , '%f' , 1 ) ;
head{1} = fgetl(fid) ;

for DDK = 1:TotNum_DM
    for j = 1:3
        head{1} = fgetl(fid) ;
    end
    [PolyDegN_DM(:,DDK)] = fscanf ( fid , '%f' , 2 ) ;
    for j = 1:2
        head{1} = fgetl(fid) ;
    end
    for j = 1:4
        [DM_Vertex(j,1:2,DDK)] = fscanf ( fid , '%f' ,  2 ) ;
    end
    for j = 1:2
        head{1} = fgetl(fid) ;
    end
    [DM_Connect(1:2,1:4,DDK)] = fscanf ( fid , '%f' , [2 4] ) ;
    head{1} = fgetl(fid) ;
end
fclose(fid) ;
Inx=[1 2 3 4 1];
figure(1)
hold on
axis([-1.5 1.5 -1.5 1.5]);
dx=2/sqrt(TotNum_DM)
Nx=round(sqrt(TotNum_DM))
for DDK = 1:TotNum_DM
plot(DM_Vertex(Inx,1,DDK),DM_Vertex(Inx,2,DDK));
cx=sum(DM_Vertex(1:4,1,DDK))/4-0.05*dx;
cy=sum(DM_Vertex(1:4,2,DDK))/4-0.05*dx;
text(cx,cy,num2str(DDK),'FontSize',10)
for j=1:4
ex=(DM_Vertex(Inx(j),1,DDK)+DM_Vertex(Inx(j+1),1,DDK))/2;
ey=(DM_Vertex(Inx(j),2,DDK)+DM_Vertex(Inx(j+1),2,DDK))/2;
ex=(ex*3+cx)/4;
ey=(ey*3+cy)/4;
DK=DM_Connect(1,j,DDK);
text(ex,ey,num2str(DK),'FontSize',6)
end
%DM_Connect(1:2,1:4,DDK)

 
end


[fid,message] = fopen ( NumEG , 'r' ) ;

head{1} = fgetl(fid) ; 
[TotNum_EG] = fscanf ( fid , '%f' , 1 ) ;
head{1} = fgetl(fid) ;

head{1} = fgetl(fid) ;
[EG_glo_Connect(1:5,1:TotNum_EG)] = fscanf ( fid , '%f' , [5 TotNum_EG] ); 
head{1} = fgetl(fid) ;

head{1} = fgetl(fid) ;
[TotNum_in_EG] = fscanf ( fid , '%f' , 1 ) ;
head{1} = fgetl(fid) ;


[in1_EG] = fscanf ( fid , '%f' , 1 ) ;
head{1} = fgetl(fid) ;
[in2_EG] = fscanf ( fid , '%f' , 1 ) ;
head{1} = fgetl(fid) ;

head{1} = fgetl(fid) ;
[EG_in_Connect(1:5,1:TotNum_in_EG)] = fscanf ( fid , '%f' , [5 TotNum_in_EG] ); 
if TotNum_in_EG ~= 0
    head{1} = fgetl(fid) ;
end

figure(2) %global edge index
hold on
axis([-1.5 1.5 -1.5 1.5]);
dx=2/sqrt(TotNum_DM)
for DDK = 1:TotNum_DM
plot(DM_Vertex(Inx,1,DDK),DM_Vertex(Inx,2,DDK));
cx=sum(DM_Vertex(1:4,1,DDK))/4-0.05*dx;
cy=sum(DM_Vertex(1:4,2,DDK))/4-0.05*dx;
text(cx,cy,num2str(DDK),'FontSize',10)
end
for j=1:TotNum_EG
    
  DDK=EG_glo_Connect(2,j);  
  inx_ls=EG_glo_Connect(3,j);    
ex=(DM_Vertex(Inx(inx_ls),1,DDK)+DM_Vertex(Inx(inx_ls+1),1,DDK))/2;
ey=(DM_Vertex(Inx(inx_ls),2,DDK)+DM_Vertex(Inx(inx_ls+1),2,DDK))/2;
cx=sum(DM_Vertex(1:4,1,DDK))/4-0.05*dx;
cy=sum(DM_Vertex(1:4,2,DDK))/4-0.05*dx;
text(ex,ey,num2str(j),'FontSize',8)
ex=(ex*3+cx)/4;
ey=(ey*3+cy)/4;
text(ex,ey,num2str(DDK),'FontSize',6)

  DDK=EG_glo_Connect(4,j);  
  inx_ls=EG_glo_Connect(5,j);    
ex=(DM_Vertex(Inx(inx_ls),1,DDK)+DM_Vertex(Inx(inx_ls+1),1,DDK))/2;
ey=(DM_Vertex(Inx(inx_ls),2,DDK)+DM_Vertex(Inx(inx_ls+1),2,DDK))/2;
cx=sum(DM_Vertex(1:4,1,DDK))/4-0.05*dx;
cy=sum(DM_Vertex(1:4,2,DDK))/4-0.05*dx;
ex=(ex*3+cx)/4;
ey=(ey*3+cy)/4;
text(ex,ey,num2str(DDK),'FontSize',6)
end 


figure(3)
hold on
axis([-1.5 1.5 -1.5 1.5]);
dx=2/sqrt(TotNum_DM)
for DDK = 1:TotNum_DM
plot(DM_Vertex(Inx,1,DDK),DM_Vertex(Inx,2,DDK));
cx=sum(DM_Vertex(1:4,1,DDK))/4-0.05*dx;
cy=sum(DM_Vertex(1:4,2,DDK))/4-0.05*dx;
text(cx,cy,num2str(DDK),'FontSize',10)
end
for j=1:TotNum_EG
    
  DDK=EG_in_Connect(2,j);  
  inx_ls=EG_in_Connect(3,j);    
ex=(DM_Vertex(Inx(inx_ls),1,DDK)+DM_Vertex(Inx(inx_ls+1),1,DDK))/2;
ey=(DM_Vertex(Inx(inx_ls),2,DDK)+DM_Vertex(Inx(inx_ls+1),2,DDK))/2;
cx=sum(DM_Vertex(1:4,1,DDK))/4-0.05*dx;
cy=sum(DM_Vertex(1:4,2,DDK))/4-0.05*dx;
text(ex,ey,num2str(j),'FontSize',8)
ex=(ex*3+cx)/4;
ey=(ey*3+cy)/4;
text(ex,ey,num2str(DDK),'FontSize',6)

  DDK=EG_in_Connect(4,j);  
  inx_ls=EG_in_Connect(5,j);    
ex=(DM_Vertex(Inx(inx_ls),1,DDK)+DM_Vertex(Inx(inx_ls+1),1,DDK))/2;
ey=(DM_Vertex(Inx(inx_ls),2,DDK)+DM_Vertex(Inx(inx_ls+1),2,DDK))/2;
cx=sum(DM_Vertex(1:4,1,DDK))/4-0.05*dx;
cy=sum(DM_Vertex(1:4,2,DDK))/4-0.05*dx;
ex=(ex*3+cx)/4;
ey=(ey*3+cy)/4;
text(ex,ey,num2str(DDK),'FontSize',6)
end 


head{1} = fgetl(fid) ;
[TotNum_out_EG] = fscanf ( fid , '%f' , 1 ) ;
head{1} = fgetl(fid) ;
[out1_EG] = fscanf ( fid , '%f' , 1 ) ;
head{1} = fgetl(fid) ;
[out2_EG] = fscanf ( fid , '%f' , 1 ) ;
head{1} = fgetl(fid) ;
[out3_EG] = fscanf ( fid , '%f' , 1 ) ;
head{1} = fgetl(fid) ;
[out4_EG] = fscanf ( fid , '%f' , 1 ) ;
head{1} = fgetl(fid) ;

head{1} = fgetl(fid) ;
[EG_out_Connect(1:5,1:TotNum_out_EG)] = fscanf ( fid , '%f' , [5 TotNum_out_EG] ) ;
head{1} = fgetl(fid) ;

fclose(fid) ;

figure(4)
hold on
axis([-1.5 1.5 -1.5 1.5]);

DM_Connect_new=DM_Connect;
for DDK = 1:TotNum_DM
for j=1:4
    DK=DM_Connect_new(1,j,DDK);
    Diff_DK=DK-DDK;
    switch Diff_DK
        case Nx*(Nx-1) % bottom boundary
            DM_Connect_new(1,j,DDK)=BD_type;           
        case -Nx*(Nx-1) % top boundary
            DM_Connect_new(1,j,DDK)=BD_type; 
        case Nx-1 %left boundary
            DM_Connect_new(1,j,DDK)=BD_type; 
        case 1-Nx % right boundary
            DM_Connect_new(1,j,DDK)=BD_type; 
        otherwise
    end

end
end
for DDK = 1:TotNum_DM
plot(DM_Vertex(Inx,1,DDK),DM_Vertex(Inx,2,DDK));
cx=sum(DM_Vertex(1:4,1,DDK))/4-0.05*dx;
cy=sum(DM_Vertex(1:4,2,DDK))/4-0.05*dx;
text(cx,cy,num2str(DDK),'FontSize',10)
for j=1:4
ex=(DM_Vertex(Inx(j),1,DDK)+DM_Vertex(Inx(j+1),1,DDK))/2;
ey=(DM_Vertex(Inx(j),2,DDK)+DM_Vertex(Inx(j+1),2,DDK))/2;
ex=(ex*3+cx)/4;
ey=(ey*3+cy)/4;
DK=DM_Connect_new(1,j,DDK);
text(ex,ey,num2str(DK),'FontSize',8)
end
end

figure(5) %global edge index
EG_glo_Connect_new=EG_glo_Connect;
hold on
axis([-1.5 1.5 -1.5 1.5]);

for DDK = 1:TotNum_DM
plot(DM_Vertex(Inx,1,DDK),DM_Vertex(Inx,2,DDK));
cx=sum(DM_Vertex(1:4,1,DDK))/4-0.05*dx;
cy=sum(DM_Vertex(1:4,2,DDK))/4-0.05*dx;
text(cx,cy,num2str(DDK),'FontSize',10)
end

EG_in_Connect_new=zeros(1:5,2*Nx*(Nx-1));
TotNum_in_EG_new=0;
in1_EG=0;
in2_EG=0;
EG_out_Connect_new=zeros(1:5,4*Nx);
out1_Connect=zeros(1:5,Nx);
out2_Connect=zeros(1:5,Nx);
out3_Connect=zeros(1:5,Nx);
out4_Connect=zeros(1:5,Nx);
TotNum_out_EG_new=0;
out1_EG_new=0;
out2_EG_new=0;
out3_EG_new=0;
out4_EG_new=0;

for j=1:TotNum_in_EG   
DDK=EG_in_Connect(2,j);  
inx_ls=EG_in_Connect(3,j);    
DDK_adj=DM_Connect_new(1,inx_ls,DDK);

DK=EG_in_Connect(4,j);  
inx_ls2=EG_in_Connect(5,j);  
DK_adj=DM_Connect_new(1,inx_ls2,DK);

if (DDK_adj==0)    
    % outer edges    
        switch inx_ls            
            case 1
                out1_EG_new=out1_EG_new+1;
                out3_EG_new=out3_EG_new+1;
                out1_Connect(1:5,out1_EG_new)=EG_in_Connect(1:5,j);
                out1_Connect(4,out1_EG_new)=0;
               
                out3_Connect(2:3,out3_EG_new)=EG_in_Connect(4:5,j);
                out3_Connect(4,out3_EG_new)=0;
                out3_Connect(5,out3_EG_new)=EG_in_Connect(3,j);                               
            case 2
                out2_EG_new=out2_EG_new+1;
                out4_EG_new=out4_EG_new+1;
                out2_Connect(1:5,out2_EG_new)=EG_in_Connect(1:5,j);
                out2_Connect(4,out2_EG_new)=0;
               
                out4_Connect(2:3,out4_EG_new)=EG_in_Connect(4:5,j);
                out4_Connect(4,out4_EG_new)=0;
                out4_Connect(5,out4_EG_new)=EG_in_Connect(3,j);    
            case 3
                out1_EG_new=out1_EG_new+1;
                out3_EG_new=out3_EG_new+1;
                out3_Connect(1:5,out3_EG_new)=EG_in_Connect(1:5,j);
                out3_Connect(4,out3_EG_new)=0;
               
                out1_Connect(2:3,out1_EG_new)=EG_in_Connect(4:5,j);
                out1_Connect(4,out1_EG_new)=0;
                out1_Connect(5,out1_EG_new)=EG_in_Connect(3,j);                               
            case 4
                out2_EG_new=out2_EG_new+1;
                out4_EG_new=out4_EG_new+1;
                out4_Connect(1:5,out4_EG_new)=EG_in_Connect(1:5,j);
                out4_Connect(4,out4_EG_new)=0;
               
                out2_Connect(2:3,out2_EG_new)=EG_in_Connect(4:5,j);
                out2_Connect(4,out2_EG_new)=0;
                out2_Connect(5,out2_EG_new)=EG_in_Connect(3,j);                    
            otherwise
        end
        
elseif (DK_adj==0)
    % outer edges
         switch inx_ls2 
             case 1
                out1_EG_new=out1_EG_new+1;
                out3_EG_new=out3_EG_new+1;
                out3_Connect(1:5,out3_EG_new)=EG_in_Connect(1:5,j);
                out3_Connect(4,out3_EG_new)=0;
               
                out1_Connect(2:3,out1_EG_new)=EG_in_Connect(4:5,j);
                out1_Connect(4,out1_EG_new)=0;
                out1_Connect(5,out1_EG_new)=EG_in_Connect(3,j);                               
            case 2
                out2_EG_new=out2_EG_new+1;
                out4_EG_new=out4_EG_new+1;
                out4_Connect(1:5,out4_EG_new)=EG_in_Connect(1:5,j);
                out4_Connect(4,out4_EG_new)=0;
               
                out2_Connect(2:3,out2_EG_new)=EG_in_Connect(4:5,j);
                out2_Connect(4,out2_EG_new)=0;
                out2_Connect(5,out2_EG_new)=EG_in_Connect(3,j);               
            case 3
                out1_EG_new=out1_EG_new+1;
                out3_EG_new=out3_EG_new+1;
                out1_Connect(1:5,out1_EG_new)=EG_in_Connect(1:5,j);
                out1_Connect(4,out1_EG_new)=0;
               
                out3_Connect(2:3,out3_EG_new)=EG_in_Connect(4:5,j);
                out3_Connect(4,out3_EG_new)=0;
                out3_Connect(5,out3_EG_new)=EG_in_Connect(3,j);                               
            case 4
                out2_EG_new=out2_EG_new+1;
                out4_EG_new=out4_EG_new+1;
                out2_Connect(1:5,out2_EG_new)=EG_in_Connect(1:5,j);
                out2_Connect(4,out2_EG_new)=0;
               
                out4_Connect(2:3,out4_EG_new)=EG_in_Connect(4:5,j);
                out4_Connect(4,out4_EG_new)=0;
                out4_Connect(5,out4_EG_new)=EG_in_Connect(3,j);    
 
            otherwise
         end 
else
    TotNum_in_EG_new=TotNum_in_EG_new+1;
    EG_in_Connect_new(1:5,TotNum_in_EG_new)=EG_in_Connect(1:5,j);
    if mod(inx_ls,2)==1
        in1_EG=in1_EG+1;
    else
        in2_EG=in2_EG+1;
    end  
end 

end

        EG_out_Connect_new(1:5,1:out1_EG_new)=out1_Connect;
        TotNum_out_EG_new=out1_EG_new;
        EG_out_Connect_new(1:5,TotNum_out_EG_new+1:TotNum_out_EG_new+out3_EG_new)=out3_Connect;
        TotNum_out_EG_new=TotNum_out_EG_new + out3_EG_new;
        EG_out_Connect_new(1:5,TotNum_out_EG_new+1:TotNum_out_EG_new+out2_EG_new)=out2_Connect;
        TotNum_out_EG_new=TotNum_out_EG_new + out2_EG_new;
        EG_out_Connect_new(1:5,TotNum_out_EG_new+1:TotNum_out_EG_new+out4_EG_new)=out4_Connect;
        TotNum_out_EG_new=TotNum_out_EG_new + out4_EG_new
        
        
 for j=1:TotNum_out_EG_new
    
  DDK=EG_out_Connect_new(2,j);  
  inx_ls=EG_out_Connect_new(3,j);
  if (DDK ==0)
  else
ex=(DM_Vertex(Inx(inx_ls),1,DDK)+DM_Vertex(Inx(inx_ls+1),1,DDK))/2;
ey=(DM_Vertex(Inx(inx_ls),2,DDK)+DM_Vertex(Inx(inx_ls+1),2,DDK))/2;
cx=sum(DM_Vertex(1:4,1,DDK))/4-0.05*dx;
cy=sum(DM_Vertex(1:4,2,DDK))/4-0.05*dx;
text(ex,ey,num2str(j),'FontSize',8)
ex=(ex*3+cx)/4;
ey=(ey*3+cy)/4;
text(ex,ey,num2str(DDK),'FontSize',6)
  end
  DDK=EG_out_Connect_new(4,j);  
  inx_ls=EG_out_Connect_new(5,j);
    if (DDK ==0)
  else
ex=(DM_Vertex(Inx(inx_ls),1,DDK)+DM_Vertex(Inx(inx_ls+1),1,DDK))/2;
ey=(DM_Vertex(Inx(inx_ls),2,DDK)+DM_Vertex(Inx(inx_ls+1),2,DDK))/2;
cx=sum(DM_Vertex(1:4,1,DDK))/4-0.05*dx;
cy=sum(DM_Vertex(1:4,2,DDK))/4-0.05*dx;
ex=(ex*3+cx)/4;
ey=(ey*3+cy)/4;
text(ex,ey,num2str(DDK),'FontSize',6)
  end
end    


figure(6)
hold on
axis([-1.5 1.5 -1.5 1.5]);

for DDK = 1:TotNum_DM
plot(DM_Vertex(Inx,1,DDK),DM_Vertex(Inx,2,DDK));
cx=sum(DM_Vertex(1:4,1,DDK))/4-0.05*dx;
cy=sum(DM_Vertex(1:4,2,DDK))/4-0.05*dx;
text(cx,cy,num2str(DDK),'FontSize',10)
end
for j=1:TotNum_in_EG_new
    
  DDK=EG_in_Connect_new(2,j);  
  inx_ls=EG_in_Connect_new(3,j);    
ex=(DM_Vertex(Inx(inx_ls),1,DDK)+DM_Vertex(Inx(inx_ls+1),1,DDK))/2;
ey=(DM_Vertex(Inx(inx_ls),2,DDK)+DM_Vertex(Inx(inx_ls+1),2,DDK))/2;
cx=sum(DM_Vertex(1:4,1,DDK))/4-0.05*dx;
cy=sum(DM_Vertex(1:4,2,DDK))/4-0.05*dx;
%text(ex,ey,num2str(j),'FontSize',8)
ex=(ex*3+cx)/4;
ey=(ey*3+cy)/4;
text(ex,ey,num2str(DDK),'FontSize',8)

  DDK=EG_in_Connect_new(4,j);  
  inx_ls=EG_in_Connect_new(5,j);    
ex=(DM_Vertex(Inx(inx_ls),1,DDK)+DM_Vertex(Inx(inx_ls+1),1,DDK))/2;
ey=(DM_Vertex(Inx(inx_ls),2,DDK)+DM_Vertex(Inx(inx_ls+1),2,DDK))/2;
cx=sum(DM_Vertex(1:4,1,DDK))/4-0.05*dx;
cy=sum(DM_Vertex(1:4,2,DDK))/4-0.05*dx;
ex=(ex*3+cx)/4;
ey=(ey*3+cy)/4;
text(ex,ey,num2str(DDK),'FontSize',8)
end 

EG_in_Connect_new(1,1:TotNum_in_EG_new)=1:TotNum_in_EG_new;
EG_glo_Connect_new(1:5,1:TotNum_in_EG_new)=EG_in_Connect_new;
TotNum_EG_new=TotNum_in_EG_new+TotNum_out_EG_new;
EG_out_Connect_new(1,1:TotNum_out_EG_new)=TotNum_in_EG_new+1:TotNum_EG_new;
EG_glo_Connect_new(1:5,TotNum_in_EG_new+1:TotNum_EG_new)=EG_out_Connect_new;

% Output
NumDM_new = sprintf('%sNdom.in',NumData_in) ;
NumEG_new = sprintf('%sNedg.in',NumData_in) ;

[fid,message] = fopen ( NumDM_new , 'w' ) ;
fprintf ( fid , '========= Input Parameters ========= \n') ;
fprintf ( fid , '%i   !Total Number of Domain \n' , TotNum_DM) ;
fprintf ( fid , '%i   !PolyDegN_max(1) \n' , PolyDegN_max(1) ) ;
fprintf ( fid , '%i   !PolyDegN_max(2) \n' , PolyDegN_max(2) ) ;
fprintf ( fid , '------------------------------------ \n') ;

for DDK = 1:TotNum_DM

fprintf ( fid , '---------- < Domain      %i  > ---------- \n', DDK) ;
fprintf ( fid , 'x12345678x12345678--------------------------\n') ;
fprintf ( fid , ' %8i %8i \n' , PolyDegN_DM(:,DDK));
fprintf ( fid , 'xxx-----.--------xxx-----.------------ \n') ;
    for j = 1:4
fprintf ( fid , '   %14.8f   %14.8f \n' ,DM_Vertex(j,1:2,DDK)) ;
    end
fprintf ( fid , 'x12345678x12345678--------------------------\n') ;
    for j = 1:4
fprintf ( fid , ' %8i %8i \n' ,DM_Connect_new(1:2,j,DDK)) ;
    end
    fprintf ( fid , '------------------------------------ \n') ;
end
fclose(fid) ;
 
[fid,message] = fopen ( NumEG_new , 'w' ) ;
fprintf ( fid , '========= Input Parameters ========= \n') ;
fprintf ( fid , '%i   !Total Number of global edge  \n' , TotNum_EG_new) ;
fprintf ( fid , '------------------------------------ \n') ;
for DDK = 1:TotNum_EG_new
fprintf ( fid , '   %8i   %8i   %8i   %8i   %8i \n' , EG_glo_Connect_new(1:5,DDK));
end
fprintf ( fid , '------------------------------------ \n') ;
fprintf ( fid , '%i   !Total Number of interior edge \n' , TotNum_in_EG_new); 
fprintf ( fid , '%i   !Total Number of interior edge 1 \n' , in1_EG); 
fprintf ( fid , '%i   !Total Number of interior edge 2 \n' , in2_EG); 
fprintf ( fid , '------------------------------------ \n') ;
for DDK = 1:TotNum_in_EG_new
fprintf ( fid , '   %8i   %8i   %8i   %8i   %8i \n' , EG_in_Connect_new(1:5,DDK));
end
fprintf ( fid , '------------------------------------ \n') ;
fprintf ( fid , '%i   !Total Number of outer edge \n' , TotNum_out_EG_new); 
fprintf ( fid , '%i   !Total Number of outer edge 1 \n' , out1_EG_new); 
fprintf ( fid , '%i   !Total Number of outer edge 2 \n' , out2_EG_new);
fprintf ( fid , '%i   !Total Number of outer edge 3 \n' , out3_EG_new); 
fprintf ( fid , '%i   !Total Number of outer edge 4 \n' , out4_EG_new); 
fprintf ( fid , '------------------------------------ \n') ;
for DDK = 1:TotNum_out_EG_new
fprintf ( fid , '   %8i   %8i   %8i   %8i   %8i \n' , EG_out_Connect_new(1:5,DDK));
end
fprintf ( fid , '------------------------------------ \n') ;
fprintf ( fid , '------------------------------------ \n') ;
fclose(fid) ;
% [TotNum_in_EG] = fscanf ( fid , '%f' , 1 ) ;
% head{1} = fgetl(fid) ;
% 
% 
% [in1_EG] = fscanf ( fid , '%f' , 1 ) ;
% head{1} = fgetl(fid) ;
% [in2_EG] = fscanf ( fid , '%f' , 1 ) ;
% head{1} = fgetl(fid) ;
% 
% head{1} = fgetl(fid) ;
% [EG_in_Connect(1:5,1:TotNum_in_EG)] = fscanf ( fid , '%f' , [5 TotNum_in_EG] ); 
% if TotNum_in_EG ~= 0
%     head{1} = fgetl(fid) ;
% end

% head{1} = fgetl(fid) ;
% [TotNum_out_EG] = fscanf ( fid , '%f' , 1 ) ;
% head{1} = fgetl(fid) ;
% [out1_EG] = fscanf ( fid , '%f' , 1 ) ;
% head{1} = fgetl(fid) ;
% [out2_EG] = fscanf ( fid , '%f' , 1 ) ;
% head{1} = fgetl(fid) ;
% [out3_EG] = fscanf ( fid , '%f' , 1 ) ;
% head{1} = fgetl(fid) ;
% [out4_EG] = fscanf ( fid , '%f' , 1 ) ;
% head{1} = fgetl(fid) ;
% 
% head{1} = fgetl(fid) ;
% [EG_out_Connect(1:5,1:TotNum_out_EG)] = fscanf ( fid , '%f' , [5 TotNum_out_EG] ) ;
% head{1} = fgetl(fid) ;
% 
% fclose(fid) ;

% [name]    :: PolyDegN_DM(Index1,Index2)
% [size]    :: (1:2,1:TotNum_DM)
% [Purpose] :: Store Polynomial Degree used in each domain
% [Detail]  :: Index2 is used to represent domain number
%              Index1 is used to distinguish the number of Polynomial
%              Degree used in x1 and x2 directions
%    
% Example   :: PolyDegN_DM(1,DDK) stores the following :
%              The degree of the polynomial used in the x1 direction of
%              Domain DDK
%
%              PolyDegN_DM(2,DDK) stores the following :
%              The degree of the polynomial used in the x2 direction of
%              Domain DDK
%              
% [name]    :: DM_Vertex(Index1,Index2)
% [size]    :: (1:5,1:2,1:TotNum_DM)
% [Purpose] :: To store the 4 physical coordinates of each domain
% [Detail]  :: The first index indicate the four vertex the fifth one
%              is a copy of the first one for coding convience
%
%              The second index stands for x and y (1 for x ,and 2 for y)
%              
%              The third index is the number of the Domain
%
% [name]    :: DM_Connect(Index1,Index2,Index3)
% [size]    :: (1:3,1:4,1:TotNum_DM)
% [Purpose] :: This variable is used to store the required info for
%              patching B.C. between two connect domains
% [Detail]  :: Index3 is used to denoted the domain number
%
%              Index2 is used to denoted the 4 edges :
%                1. xi_2 = -1
%                2. xi_1 =  1
%                3. xi_2 =  1
%                4. xi_1 = -1
%
%              Index1 is used to denoted the store info :
%                1 for the number of the connected domain
%                2 for the connecting edge
%                3 for patching direction
% 
% Example   :: DM_Connect(1,1,DDK) stores the following :
%              The number of the connected domain on the
%              1st edge of domain DDK
% 
%              DM_Connect(2,1,DDK) stores the following :
%              The edge number of the connected domain on the
%              1st edge of domain DDK
% 
%              DM_Connect(3,1,DDK) stores the following :
%              The patching direction of the connected domain on the
%              1st edge of domain DDK. (1 for the same ,2 for the oppsite)
%
%  Domain example :
%
%      ______________________________________
%      |           |           |            |
%      |     7     |     8     |     9      |
%      |           |           |            |
%      _________________(3)__________________
%      |           |           |            |
%      |     4    (4)    5    (2)    6      |
%      |           |           |            |
%      _________________(1)__________________
%      |           |           |            |
%      |     1     |     2     |     3      |
%      |           |           |            |
%      ______________________________________





