function dose2bed(ID)

%ID: Nombre de la carpeta donde se encuentran los archivos de Imagen,
%Dosis, Estucturas y Plan.

%%%%%%%%%%%%%%%%%%%%%%%%% DATOS EDITABLES %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Ruta a la carpeta raiz donde está la carpeta ID:
    %Ejemplo: Dosis, Estructuras plan e imagenes en: C:\Pacientes\Paciente1
    %ID='Paciente1'
    %Ruta='C:\Pacientes\'

Ruta='W:\RF\10 Personales\Jaime\CNN\Prueba\Dosis BED\Pacientes\';

%Dos posibilidades:
    %Interpolado de la dosis: Reduce generalmente el tamaño de la malla de cálculo:
    %P00=1
    %Interpolado de las estructuras: Mantiene la malla de cálculo, el histograma
    %presentará efectos de borde. P00=0

P00=1;

%Valor del cociente alfa/beta utilizado por defecto

AlfaB=3;

%Tipo de método de interpolación deseado:
    %nearest
    %cubic
    %linear

Metodo='linear';


%Tabla de equivalencias entre el nombre de la ROI y el valor alpha/beta
%asignado 
%SEGM{*,1}='Nombre de la estructura etiquetada en el TPS'
%SEGM{*,2}=Valor del cociente alfa/beta asignado
%El sistema permite añadir el número deseado de estructuras:

%SEGM{n+1,1}='organo';
%SEGM{n+1,2}=4;

SEGM{1,1}='GTV';
SEGM{1,2}=10;
SEGM{2,1}='CTV';
SEGM{2,2}=10;
SEGM{3,1}='FemoralHead (Left)';
SEGM{3,2}=4;
SEGM{4,1}='Bladder';
SEGM{4,2}=4;
SEGM{5,1}='Rectum';
SEGM{5,2}=4;
SEGM{6,1}='SmallBowel';
SEGM{6,2}=4;
SEGM{7,1}='Colon';
SEGM{7,2}=4;
SEGM{8,1}='Kidney (Left)';
SEGM{8,2}=2;
SEGM{9,1}='Kidney (Right)';
SEGM{9,2}=2;
SEGM{10,1}='FemoralHead (Right)';
SEGM{10,2}=4;
SEGM{11,1}='Lungs';
SEGM{11,2}=4;
SEGM{12,1}='Lung_L';
SEGM{12,2}=4;
SEGM{13,1}='Lung_R';
SEGM{13,2}=4;
SEGM{14,1}='ContralateralBreast';
SEGM{14,2}=3.3;
SEGM{15,1}='Stomach';
SEGM{15,2}=3;
SEGM{16,1}='ThyroidGland';
SEGM{16,2}=3;
SEGM{17,1}='ParotidGland (Right)';
SEGM{17,2}=3;
SEGM{18,1}='SpinalCanal';
SEGM{18,2}=2;
SEGM{19,1}='SpinalCord';
SEGM{19,2}=2;
SEGM{20,1}='Prostate';
SEGM{20,2}=1.5;
SEGM{21,1}='Eye, Left';
SEGM{21,2}=1;
SEGM{22,1}='Eye, Right';
SEGM{22,2}=1;



%%%%%%%%%%%%%%%%%%%%%%FIN DATOS EDITABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


PXX=strcat(Ruta,ID,'\*.*');
Lista=dir(PXX);
[F,~]=size(Lista);

%Inicializamos algunas variables

Mesa=[];
List={};
t=1;

%Detectamos dentro de la carpeta indicada el archivo de dosis, el de
%estructuras, el plan y las imágenes, obteniendo de cada uno la información
%necesaria
for i=1:1:F
    Nombre=Lista(i,1).name;
    if contains(lower(Nombre),'dcm')==1
        dicom=dicominfo(strcat(Ruta,ID,'\',Nombre));
        Modalidad=dicom.Modality;
        if strcmpi(Modalidad,'RTSTRUCT')==1
            SR=dicom;     
        elseif strcmpi(Modalidad,'RTPLAN')==1
            PN=dicom;       
            fx=PN.FractionGroupSequence.Item_1.NumberOfFractionsPlanned;%Numero de fracciones del tratamiento        
        elseif strcmpi(Modalidad,'RTDOSE')==1
            DS=Nombre;        
            Dosis=dicomread(strcat(Ruta,ID,'\',DS));%Lectura del cubo de dosis
            [Xd,Yd,~,Zd]=size(Dosis);
            Dosis_M=zeros(Xd,Yd,Zd);
            Dosis_M(:,:,:)=Dosis(:,:,1,:); 
            Dosis_D=dicom;           
            Espesor=Dosis_D.SliceThickness;%Tamaño del grid de dosis en la dimensión longitudinal
            Grid=Dosis_D.PixelSpacing;%Tamaño del grid de dosis en las dimensiones lateral y vertical
            Escala=Dosis_D.DoseGridScaling; %Factor de escalado entre 16bits y la dosis del plan            
            if isempty(Espesor)==1
                Espesor=Grid(1,1);
            end            
            IPPds=Dosis_D.ImagePositionPatient; %Origen de posición del paciente
            IOPds=Dosis_D.ImageOrientationPatient; %Vector de orientación del paciente (HFS,HFP,FFS...)
            Wdt=double(Dosis_D.Width);
            Hgt=double(Dosis_D.Height);
            Lgt=double(Dosis_D.NumberOfFrames);
        elseif strcmpi(Modalidad,'CT')==1            
            IM=dicom;            
            IPP=IM.ImagePositionPatient;
            Mesa=[Mesa;IPP(3,1) t];
            List{end+1}=IM.SOPInstanceUID;
            UH1(:,:,t)=dicomread(strcat(Ruta,ID,'\',Lista(i,1).name));
            R=IM.Rows;
            C=IM.Columns;
            Ps=IM.PixelSpacing;
            IOP=IM.ImageOrientationPatient; %Vector de orientación del paciente (HFS,HFP,FFS...)
            St=IM.SliceThickness;
            t=t+1;
        end
    end
end    

%Ordenamos las imagenes por posición de mesa
Mesa=sortrows(Mesa,'ascend');
[F]=size(Mesa,1);
Listado={};
UH=[];
for i=1:F
   Listado{end+1}=List{1,Mesa(i,2)};
   UH(:,:,i)=UH1(:,:,Mesa(i,2));
end

IPP(3,1)=Mesa(1,1);


Diccionario={};
for j=1:length(SEGM);
    Diccionario{1,j}=SEGM{j,1};
end


Estructuras=fieldnames(SR.StructureSetROISequence);%Numero de ROIs presentes en el estudio
[F]=size(Estructuras,1);
OARs=[];
%Matriz para cambio de base entre sistemas de referencia (sala-paciente)

%%%%ROTACION DE LA IMAGEN SI LA HUBIERA
V1=IOP(1:3,1);
V2=IOP(4:6,1); 
V3=cross(V1,V2);
ROT=[V1 V2 V3];
ImaP=(ROT^-1)*IPP;

%%%%ROTACION EL CUBO DE DOSIS SI LA HUBIERA
V1=IOPds(1:3,1);
V2=IOPds(4:6,1); 
V3=cross(V1,V2);
ROTds=[V1 V2 V3];
ImaPds=(ROTds^-1)*IPPds;




IPP=ImaP;%Posición del estudio de imagen referido a las coordenadas de la sala
IPPds=ImaPds;%posición del cubo de dosis referido a las coordenadas de la sala


%Estructuras presentes en el estudio y con alpha/betas distintos de 3

for i=1:1:F
    Campo=strcat('Item_',num2str(i));
    Campo=SR.StructureSetROISequence.(Campo);
        if isempty(find(strcmp(Diccionario,Campo.ROIName)==1))==0
            OARs=[OARs;i find(strcmp(Diccionario,Campo.ROIName)==1)];
        end
end

%Buscamos el número de imagenes del estudio
N_IM=numel(fieldnames(SR.ReferencedFrameOfReferenceSequence.Item_1.RTReferencedStudySequence.Item_1.RTReferencedSeriesSequence.Item_1.ContourImageSequence));

%almacenamos los ReferencedSOPInstanceUID para identificar las imagenes
SOP={};
for i=1:N_IM
    Imagen=strcat('Item_',num2str(i));
    Campo=SR.ReferencedFrameOfReferenceSequence.Item_1.RTReferencedStudySequence.Item_1.RTReferencedSeriesSequence.Item_1.ContourImageSequence.(Imagen);
    SOP{1,i}=Campo.ReferencedSOPInstanceUID;
end


%Generamos la matriz 3D de ceros donde iremos colocando los indices
M=zeros(R,C,N_IM);

%Recomponemos las estructuras del diccionario
for i=1:size(OARs,1)
    Campo=strcat('Item_',num2str(OARs(i,1)));
    Campo=SR.ROIContourSequence.(Campo);
    Campo=Campo.ContourSequence;
    CortesROI=fieldnames(Campo);
    [F]=size(CortesROI,1);
    for j=1:1:F
        Corte=strcat('Item_',num2str(j));
        Corte=Campo.(Corte);
        IMA=Corte.ContourImageSequence.Item_1.ReferencedSOPInstanceUID;
        Plano=find(strcmp(Listado,IMA)==1);
        MAT=Corte.ContourData; 
        if isempty(MAT)==0
            [fm]=size(MAT,1);
            mm=[];
            for k=1:fm/3
                VC=[MAT(3*k-2,1);MAT(3*k-1,1);MAT(3*k,1)];
                VC=(ROT^-1)*VC;
                mm=[mm; VC(1,1) VC(2,1)];            
            end     
            mm(:,2)=round(mm(:,2)./Ps(2,1))-round(ImaP(2,1)/Ps(2,1));
            mm(:,1)=round(mm(:,1)./Ps(1,1))-round(ImaP(1,1)/Ps(1,1));
            RR=poly2mask(mm(:,1),mm(:,2),double(R),double(C));                        
            m=double(RR);
        end
        Plano_N=M(:,:,Plano);
        Plano_N(Plano_N==0)=m(Plano_N==0)*double(SEGM{OARs(i,2),2}); % Le damos a cada pixel el valor del alpha/beta según el órgano
        M(:,:,Plano)=Plano_N;
    end
end
M(M==0)=AlfaB;  %Alpha beta = 3 para todo el tejido no etiquetado

[Xs,Ys,Zs]=size(M);

if P00==0
    %Resampleamos la matriz de estructuras a las dimensiones de la matriz de
    %dosis
    Mrs=imresize3(M,[round(Xs*Ps(1,1)/Grid(1,1)) round(Ys*Ps(2,1)/Grid(2,1)) round(Zs*St/Espesor)],Metodo);
    Mrs1=Mrs;
elseif P00==1
    %Resampleamos la matriz de dosis a las dimensiones de la matriz de
    %de estructuras
    
    [Xd,Yd,Zd]=size(Dosis_M);   
    xq=(0:Ps(1,1)/Grid(1,1):Hgt);    
    yq=(0:Ps(1,1)/Grid(2,1):Wdt);
    zq=(0:St/Espesor:Lgt);    
    F=griddedInterpolant(Dosis_M,Metodo);
    Dosis_sm=F({xq,yq,zq});
    [He,Wd,Fr]=size(Dosis_sm);
    Dosis_M=Dosis_sm;
    Mrs=M;
    
else
    disp('El valor de P debe ser 1 o 0.')
    return
end



%Reorientamos ambos cubos en caso de que hubiera rotaciones respecto al
%sistema de la sala
Rsr=maketform('affine',ROT);
Rds=maketform('affine',ROTds);
Mrs=imtransform(Mrs,Rsr);
Dosis_M=imtransform(Dosis_M,Rds);

[Xs,Ys,Zs]=size(Mrs);
[Xd,Yd,Zd]=size(Dosis_M);

%Alineamos los dos cubos, buscando las posiciones relativas de los origenes

%Longitudinal

%%%INICIOS EN LONGITUDINAL
if P00==1
    Zln=St;
else
    Zln=Espesor;
end
   
DeltaL=ceil((IPP(3,1)-IPPds(3,1))/Zln);
restoL=(ceil((IPP(3,1)-IPPds(3,1))/Zln)-((IPP(3,1)-IPPds(3,1))/Zln))*Zln;

if DeltaL>0 %Cubo de dosis es más largo en la dirección Z negativo
    SR_IL=1;%Inicio del cubo de estructuras
    DS_IL=DeltaL+1;    
elseif DeltaL<0 %Cubo de dosis es más corto en la dirección Z negativo
    SR_IL=abs(DeltaL)+1;
    DS_IL=1;
else
    SR_IL=1;
    DS_IL=1;
end
%%%FINALES EN LONGITUDINAL
SR_FL=Zd-DeltaL;
if SR_FL>Zs
    SR_FL=Zs;
end



%%%INICIOS EN LATERAL
if P00==1
    Xln=Ps(1,1);
else
    Xln=Grid(1,1);
end
DeltaLt=ceil((IPP(1,1)-IPPds(1,1))/Xln);
restoLt=(ceil((IPP(1,1)-IPPds(1,1))/Xln)-((IPP(1,1)-IPPds(1,1))/Xln))*Xln;
if DeltaLt>0 %Cubo de dosis es más largo en la dirección X negativo
    SR_ILt=1;%Inicio del cubo de estructuras
    DS_ILt=DeltaLt+1;
elseif DeltaLt<0 %Cubo de dosis es más corto en la dirección X negativo
    SR_ILt=abs(DeltaLt)+1;
    DS_ILt=1;
else
    SR_ILt=1;
    DS_ILt=1;
end
%%%FINALES EN LATERAL
SR_FLt=Yd-DeltaLt;
if SR_FLt>Ys
    SR_FLt=Ys;
end



%%%INICIOS EN VERTICAL
if P00==1
    Yln=Ps(2,1);
else
    Yln=Grid(2,1);
end
DeltaV=ceil((IPP(2,1)-IPPds(2,1))/Yln);
restoV=(ceil((IPP(2,1)-IPPds(2,1))/Yln)-((IPP(2,1)-IPPds(2,1))/Yln))*Yln;
if DeltaV>0 %Cubo de dosis es más largo en la dirección Y negativo
    SR_IV=1;%Inicio del cubo de estructuras
    DS_IV=DeltaV+1;
elseif DeltaV<0 %Cubo de dosis es más corto en la dirección Y negativo
    SR_IV=abs(DeltaV)+1;
    DS_IV=1;
else
    SR_IV=1;
    DS_IV=1;
end
%%%FINALES EN VERTICAL
SR_FV=Xd-DeltaV;
if SR_FV>Xs
    SR_FV=Xs;
end



%Creamos el cubo con coeficientes alpha/beta de igual tamaño que la malla
%de dosis y con los coeficientes dependiendo de las estrcuturas
SR_d=3.*ones(Xd,Yd,Zd);
xi=0;
for i=SR_IV:SR_FV
    yi=0;
    for j=SR_ILt:SR_FLt
        zi=0;
        for k=SR_IL:SR_FL
            SR_d(DS_IV+xi,DS_ILt+yi,DS_IL+zi)=Mrs(i,j,k);
            zi=zi+1;
        end
        yi=yi+1;
    end
    xi=xi+1;
end
 


[f,c,p]=size(SR_d);
BED=double(Dosis);
EQD2=double(Dosis);
Dosis_M=double(Dosis_M)*Escala;

%Cálculo de la DBE y el EQD2
for i=1:f
    for j=1:c
        for k=1:p
            alpha=SR_d(i,j,k);  
            BED(i,j,1,k)=double((Dosis_M(i,j,k))*(1+((Dosis_M(i,j,k)/fx)./alpha)));
            EQD2(i,j,1,k)=double(BED(i,j,1,k)./(1+(2/alpha)));          
        end
    end
end

save('Mat1.mat','SR_d');


Escala1=max(max(max(BED)))/65535;
Escala2=max(max(max(EQD2)))/65535;
BED=uint16(BED./Escala1);
EQD2=uint16(EQD2./Escala2);


Nombre=strrep(DS,'.dcm','');
DS=strcat(Nombre,'_BED','.dcm');
EQ=strcat(Nombre,'_EQD2','.dcm');

Dosis_D.DoseGridScaling=Escala1;


if P00==1
    %Anotamos los nuevos valores del voxel de dosis en su cabecera
    Dosis_D.PixelSpacing=Ps;
    Dosis_D.SliceThickness=St;
    z=(0:St:Fr*St-1);   
    Dosis_D.Width=Wd;
    Dosis_D.Height=He;
    Dosis_D.NumberOfFrames=Fr;
    Dosis_D.GridFrameOffsetVector=z;
    Dosis_D.ImagePositionPatient=IPPds-[restoLt;restoV;restoL];
end


Ruta2=strcat(Ruta,ID,'\',DS);
dicomwrite(BED,Ruta2,Dosis_D,'CreateMode','copy');
Dosis_D.DoseGridScaling=Escala2;
Ruta2=strcat(Ruta,ID,'\',EQ);
dicomwrite(EQD2,Ruta2,Dosis_D,'CreateMode','copy');


%%% dose2bed v1.0. Hospital HM Sanchinarro
%%% Jaime Martí Asenjo





