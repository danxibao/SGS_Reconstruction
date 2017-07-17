// RSCalculation.cpp: implementation of the RSCalculation class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SGS_Reconstruction.h"
#include "RSCalculation.h"
#include <time.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

RSCalculation::RSCalculation()
{

}

RSCalculation::~RSCalculation()
{

}

//射线在废物桶中的衰减长度
double RSCalculation::AttenuationLength(CartesianCoordinate PA,double dYOfDrum,double HOfDrum,double ThicknessOfDrum) 
{
	CartesianCoordinate PB;
	PB.x=DecSetupData.GHPGDX;
	PB.y=-1*dYOfDrum;
	PB.z=HOfDrum;
	double R=DecSetupData.GWDD*0.5+ThicknessOfDrum;
	double a=pow(PB.x-PA.x,2)+pow(PB.y-PA.y,2);
	double b=2*((PB.x-PA.x)*PA.x+(PB.y-PA.y)*PA.y);
	double c=pow(PA.x,2)+pow(PA.y,2)-R*R;
    double d=b*b-4*a*c;
	double k1,k2;
    if(d<=0) return 0;
	if(a<0.00000001) return 0;
	k1=(-1*b-sqrt(d))/(2*a);
	k2=(-1*b+sqrt(d))/(2*a);    
	double L=sqrt(pow(PA.x-PB.x,2)+pow(PA.y-PB.y,2)+pow(PA.z-PB.z,2));
	if(k1<0) return L*k2;//表示点源A在废物桶内
	else return L*(k2-k1);//表示点源A在废物桶外
}
double RSCalculation::AttenuationLength2(CartesianCoordinate PA,double dYOfDrum,double HOfDrum,double D) 
{
	CartesianCoordinate PB;
	PB.x=DecSetupData.GHPGDX;
	PB.y=-1*dYOfDrum;
	PB.z=HOfDrum;
	double R=0.5*D;
	double a=pow(PB.x-PA.x,2)+pow(PB.y-PA.y,2);
	double b=2*((PB.x-PA.x)*PA.x+(PB.y-PA.y)*PA.y);
	double c=pow(PA.x,2)+pow(PA.y,2)-R*R;
    double d=b*b-4*a*c;
	double k1,k2;
    if(d<=0) return 0;
	if(a<0.00000001) return 0;
	k1=(-1*b-sqrt(d))/(2*a);
	k2=(-1*b+sqrt(d))/(2*a);    
	double L=sqrt(pow(PA.x-PB.x,2)+pow(PA.y-PB.y,2)+pow(PA.z-PB.z,2));
	if(k1<0) return L*k2;//表示点源A在废物桶内
	else return L*(k2-k1);//表示点源A在废物桶外
}
//射线在废物桶各段层中的衰减长度
void RSCalculation::AttenuationLength(CartesianCoordinate PA,double dYOfDrum,double HOfDrum,double ThicknessOfDrum,double LInSegments[5]) 
{
	for(int i=0;i<5;i++) 
		LInSegments[i]=0;
	CartesianCoordinate PB;
	PB.x=DecSetupData.GHPGDX;
	PB.y=-1*dYOfDrum;
	PB.z=HOfDrum;
	double L0=AttenuationLength(PA,dYOfDrum,HOfDrum,ThicknessOfDrum);//在桶内的长度
	double L=sqrt(pow(PA.x-PB.x,2)+pow(PA.y-PB.y,2)+pow(PA.z-PB.z,2));//AB两点间的距离
	double dH=DecSetupData.GWDH/DecSetupData.DSNumZ;//段层高度
    double k1,k2;
    if(fabs(PA.z-PB.z)<=0.5*dH) {LInSegments[2]=L0;}//如果在点B在当前层
    else if(PA.z>PB.z+0.5*dH&&PA.z<=PB.z+1.5*dH)//如果在点B在偏上一层
	{
		k1=0.5*dH/(PA.z-PB.z);
		LInSegments[3]=L*(1-k1);
		LInSegments[2]=L0-LInSegments[3];	
	}
    else if(PA.z<PB.z-0.5*dH&&PA.z>=PB.z-1.5*dH)//如果在点B在偏下一层
	{
		k1=-0.5*dH/(PA.z-PB.z);
		LInSegments[1]=L*(1-k1);
		LInSegments[2]=L0-LInSegments[1];	
	}
    else if(PA.z>PB.z+1.5*dH)//如果在点B在偏上二层
	{
		k1=0.5*dH/(PA.z-PB.z);
		k2=1.5*dH/(PA.z-PB.z);	
		LInSegments[4]=L*(1-k2);
		LInSegments[3]=L*(1-k1)-LInSegments[4];
		LInSegments[2]=L0-LInSegments[3]-LInSegments[4];	
	}
    else if(PA.z>PB.z+1.5*dH)//如果在点B在偏下二层
	{
		k1=-0.5*dH/(PA.z-PB.z);
		k2=-1.5*dH/(PA.z-PB.z);	
		LInSegments[0]=L*(1-k2);
		LInSegments[1]=L*(1-k1)-LInSegments[0];
		LInSegments[2]=L0-LInSegments[0]-LInSegments[1];	
	}
}


//铅的质量衰减系数（cm2/g），能量范围0.1MeV至5MeV，Engergy单位为MeV
double RSCalculation::CAttenuationPerMassOfPb(double Energy) 
{
	return 10*exp(2.2003*pow(log10(Energy),2.0)-2.1953*pow(log10(Energy),1.0)-4.9528);
}
//轻质材料的质量衰减系数（cm2/g），能量范围0.1MeV至5MeV，Engergy单位为MeV
double RSCalculation::CAttenuationPerMass(double Energy) 
{
	return 10*exp(-1.0643*log10(Energy)-5.0829);
}

//////////////////////////////////////////////////////////////////////////////
//重建计算//
//////////////////////////////////////////////////////////////////////////////
//初始化
void RSCalculation::INIReconCal() 
{
	m_DensityOfWD=NULL;
	m_ActivityOfWD=NULL;
	m_TDHead=NULL;
	m_EDHead=NULL;
	CoefArrayHead=NULL;
	CoefArrayHead2=NULL;
	CoefArrayHead3=NULL;
	RelativeEHead=NULL;
	RadiusOfSourceHead=NULL;
	RS_E_Head=NULL;
	RS_E_Head2=NULL;
	EMResidHead=NULL;
	XHead=NULL;
	XHead2=NULL;
}
//输出结果
void RSCalculation::OutputReconResults() 
{
	//TriDataB.outputdata_3d(CoefArrayHead,m_Folder+CoefArrayHead->Title,_T("xyz"));
	TriDataB.outputdata_2d2(CoefArrayHead,m_Folder+_T("A.txt"),_T("YZ"),0);
	TriDataB.outputdata_1d(EMResidHead,m_Folder+EMResidHead->Title,-1,0,0);
	TriDataB.outputdata_3d(m_DensityOfWD,m_Folder+m_DensityOfWD->Title,_T("xyz"));
	TriDataB.outputdata_3d(m_ActivityOfWD,m_Folder+m_ActivityOfWD->Title,_T("xyz"));
}
//清除重建过程中的所有指针
void RSCalculation::DeleteAllPointers() 
{
	TriDataB.DeleteNodes(m_DensityOfWD);m_DensityOfWD=NULL;
	TriDataB.DeleteNodes(m_ActivityOfWD);m_ActivityOfWD=NULL;
	TriDataB.DeleteNodes(m_EDHead);m_EDHead=NULL;
	TriDataB.DeleteNodes(m_EDHead);m_TDHead=NULL;
	TriDataB.DeleteNodes(CoefArrayHead);CoefArrayHead=NULL;
	TriDataB.DeleteNodes(EMResidHead);EMResidHead=NULL;
	TriDataB.DeleteNodes(CoefArrayHead2);CoefArrayHead2=NULL;
	TriDataB.DeleteNodes(RelativeEHead);RelativeEHead=NULL;
	TriDataB.DeleteNodes(RadiusOfSourceHead);RadiusOfSourceHead=NULL;
	TriDataB.DeleteNodes(RS_E_Head);RS_E_Head=NULL;
	TriDataB.DeleteNodes(RS_E_Head2);RS_E_Head2=NULL;
	TriDataB.DeleteNodes(XHead);XHead=NULL;
	TriDataB.DeleteNodes(XHead2);XHead2=NULL;


}
//////////////////////////////////////////////////////////读取探测设置
bool RSCalculation::ReadDetectionSetup(CString FullFileName)
{
	//初值
	DecSetupData.DetectionType=_T("SGS");
	DecSetupData.nDetectionType=0;//项目类型：0SGS,1TGS，2ISGS，3STGS
	DecSetupData.TypeOfWD=_T("Steel200L");	//Type of the drum
	DecSetupData.TypeOfMatrix=_T("Compressure Waste");	//Type of the drum
	DecSetupData.GWDD=56;	//Diameter of the drum
	DecSetupData.GWDH=90;	//Height of the drum
	DecSetupData.GHPGDX=53;	//Initial x-coordinate of the inelt of HPG dector
	DecSetupData.GHPGDY=0;	//Initial y-coordinate of the inelt of HPG dector
	DecSetupData.GHPGDZ=5;	//Initial z-coordinate of the inelt of HPG dector
	DecSetupData.DSNumY=1;	//Times of transtation
	DecSetupData.DSNumZ=9;	//Number of segments
	DecSetupData.DSNumA=1;	//Times of rotation
	DecSetupData.DSLY=0;	//Step length of transtation
	DecSetupData.DSLZ=10;	//Step length of rise
	DecSetupData.DSLA=0;	//Step length of rotation
	DecSetupData.DSigleT=10;	//Period of single detection
	DecSetupData.VORotation=10;	//Velocity of rotation,min-1
	DecSetupData.IsTransmission=false;//Is transmission detection carried out
	DecSetupData.TransmissionSource=_T("Eu-152");//transmission source
	DecSetupData.WeightOfDrum=222;//Weight of Waste drum
	DecSetupData.TypeOfIteration=1;
	DecSetupData.MaxIteration=500;
	DecSetupData.NuGridsX=1;
	DecSetupData.NuGridsY=1;
	DecSetupData.NuGridsZ=9;
	DecSetupData.IsEqualVolumeOfGrids=true;
	//读取设置
		int i,j;
		CString CTemp,CSpe;
		char buf[300];
		char buft[300];
		int length_buf=300;
		int key=0;
		ifstream read(FullFileName);//打开文件
		if(read)
		{
			while(read)
			{///
			for(i=0;i<length_buf;i++) {buf[i]='\0';buft[i]='\0';}
			read.getline(buf,length_buf);//挨行读取文件
			CSpe=buf;
			if(CSpe.Find(_T("END"))>=0||CSpe.Find(_T("End"))>=0||CSpe.Find(_T("end"))>=0) break;//文件读取结束
			for(i=0,j=0,key=0;i<length_buf;i++) 
			{//
			if(buf[i]==';') break;
			else if(buf[i]==':'||buf[i]=='=') {key=1;}
			else if(key==1&&buf[i]!=' ') buft[j++]=buf[i];
			}//
			CTemp=buft;
			if(CSpe.Find(_T("FileName"))>=0) DecSetupData.FileName=CTemp;
			else if(CSpe.Find(_T("DetectionType"))>=0) DecSetupData.DetectionType=CTemp;
			else if(CSpe.Find(_T("NoOfWD"))>=0) DecSetupData.NoOfWD=CTemp;
			else if(CSpe.Find(_T("TypeOfWD"))>=0) DecSetupData.TypeOfWD=CTemp;
			else if(CSpe.Find(_T("TypeOfMatrix"))>=0) DecSetupData.TypeOfMatrix=CTemp;
			else if(CSpe.Find(_T("GWDD"))>=0) DecSetupData.GWDD=atof(CTemp);
			else if(CSpe.Find(_T("GWDH"))>=0) DecSetupData.GWDH=atof(CTemp);
			else if(CSpe.Find(_T("GHPGDX"))>=0) DecSetupData.GHPGDX=atof(CTemp);
			else if(CSpe.Find(_T("GHPGDY"))>=0) DecSetupData.GHPGDY=atof(CTemp);
			else if(CSpe.Find(_T("GHPGDZ"))>=0) DecSetupData.GHPGDZ=atof(CTemp);
			else if(CSpe.Find(_T("DSNumY"))>=0) DecSetupData.DSNumY=atoi(CTemp);
			else if(CSpe.Find(_T("DSNumZ"))>=0) DecSetupData.DSNumZ=atoi(CTemp);
			else if(CSpe.Find(_T("DSNumA"))>=0) DecSetupData.DSNumA=atoi(CTemp);
			else if(CSpe.Find(_T("DSLY"))>=0) DecSetupData.DSLY=atof(CTemp);
			else if(CSpe.Find(_T("DSLZ"))>=0) DecSetupData.DSLZ=atof(CTemp);
			else if(CSpe.Find(_T("DSLA"))>=0) DecSetupData.DSLA=atof(CTemp);
			else if(CSpe.Find(_T("DSigleT"))>=0) DecSetupData.DSigleT=atof(CTemp);
			else if(CSpe.Find(_T("VORotation"))>=0) DecSetupData.VORotation=atof(CTemp);
			else if(CSpe.Find(_T("IsTransmission"))>=0) {if(atoi(CTemp)>0) DecSetupData.IsTransmission=true; else DecSetupData.IsTransmission=false;}
			else if(CSpe.Find(_T("TransmissionSource"))>=0) DecSetupData.TransmissionSource=CTemp;
			else if(CSpe.Find(_T("WeightOfDrum"))>=0) DecSetupData.WeightOfDrum=atof(CTemp);
			else if(CSpe.Find(_T("TypeOfIteration"))>=0) DecSetupData.TypeOfIteration=atoi(CTemp);
			else if(CSpe.Find(_T("MaxIteration"))>=0) DecSetupData.MaxIteration=atoi(CTemp);
			else if(CSpe.Find(_T("NuGridsX"))>=0) DecSetupData.NuGridsX=atoi(CTemp);
			else if(CSpe.Find(_T("NuGridsY"))>=0) DecSetupData.NuGridsY=atoi(CTemp);
			else if(CSpe.Find(_T("NuGridsZ"))>=0) DecSetupData.NuGridsZ=atoi(CTemp);
			else if(CSpe.Find(_T("IsEqualVolumeOfGrids"))>=0) {if(atoi(CTemp)>0) DecSetupData.IsEqualVolumeOfGrids=true; else DecSetupData.IsEqualVolumeOfGrids=false;}
			}///
			read.close();
			if(DecSetupData.DetectionType.Find(_T("STGS"))>=0) DecSetupData.nDetectionType=3;
			else if(DecSetupData.DetectionType.Find(_T("ISGS"))>=0) DecSetupData.nDetectionType=2;
			else if(DecSetupData.DetectionType.Find(_T("TGS"))>=0) DecSetupData.nDetectionType=1;
			else if(DecSetupData.DetectionType.Find(_T("SGS"))>=0) DecSetupData.nDetectionType=0;
		}
		else return false;
return true;
}

//验证探测的有效性
bool RSCalculation::VerifyDetection(CString Folder) 
{
	//初始化
	INIReconCal();
	//读取探测设置
	if(!ReadDetectionSetup(Folder+_T("_DetectionSetup.txt"))) return false;
	//读取并验证透射测量数据
	if(DecSetupData.IsTransmission)//如果需要进行透射重建，读取探测文件
	{
		m_TDHead=TriDataB.ReadData(m_Folder+_T("_TransDect.dat"));//读取透射测量文件
		if(!m_TDHead) return false;//如果未能读取，返回失败
		//如果探测数据的IJK与设置不同，返回失败，测量文件的IJK依次未，转动、平动、升降
		if(m_TDHead->NuNodesInZ!=DecSetupData.DSNumZ) return false;
		//能量个数存在问题，返回失败
		if(m_TDHead->NuDataUsed<=0||m_TDHead->NuDataUsed>MaxNuDataInTriDNode) return false;
	}
	//读取并验证发射测量数据
		m_EDHead=TriDataB.ReadData(m_Folder+_T("_EmisDect.dat"));
		if(!m_EDHead) return false;//如果未能读取，返回失败
		//如果探测数据的IJK与设置不同，返回失败，测量文件的IJK依次未，转动、平动、升降
		if(m_EDHead->NuNodesInX!=DecSetupData.DSNumA) return false;
		if(DecSetupData.DetectionType!=2&&m_EDHead->NuNodesInY!=DecSetupData.DSNumY) return false;
		if(m_EDHead->NuNodesInZ!=DecSetupData.DSNumZ) return false;
		//能量个数存在问题，返回失败
		if(m_EDHead->NuDataUsed<=0||m_EDHead->NuDataUsed>MaxNuDataInTriDNode) return false;
		//SGS如果y方向测量次数不是1，返回失败
		if(DecSetupData.nDetectionType==0&&m_EDHead->NuNodesInY!=1) return false;//SGS
		//ISGS如果y方向测量次数不是2，返回失败
		if(DecSetupData.nDetectionType==2&&m_EDHead->NuNodesInY!=2) return false;//ISGS
		//STGS如果没有y方向测量次数不大于1，返回失败
		if(DecSetupData.nDetectionType==3&&m_EDHead->NuNodesInY<=1) return false;//STGS

	//创建密度结果指针，如果是SGS则段层内密度一样，TGS则分网格
		if(DecSetupData.nDetectionType==1) //TGS
		{
			 if(!CreatDensityOfWD(DecSetupData.NuGridsX,DecSetupData.NuGridsY,DecSetupData.NuGridsZ)) return false;
		}
		else 
		{	
			if(!CreatDensityOfWD(1,1,DecSetupData.NuGridsZ)) return false;//SGS
		}
	//创建活度结果指针
		if(DecSetupData.nDetectionType==2&&!CreatActivityHead(1,1,DecSetupData.NuGridsZ)) return false;//ISGS
		else if(!CreatActivityHead(DecSetupData.NuGridsX,DecSetupData.NuGridsY,DecSetupData.NuGridsZ)) return false;

	return true;
}
//对废物桶创建密度矩阵
bool RSCalculation::CreatDensityOfWD(int NuX,int NuY,int NuZ)
{
	m_DensityOfWD=TriDataB.CreatNodes(NuX,NuY,NuZ);
	if(!m_DensityOfWD) return false;
	m_DensityOfWD->Title=_T("Results_Density.dat");//数据空间的名称，可修改;
	int i,j,k,ndata,NData=1;
	double r1,r2;
	//如果透射测量文件存在，按照透射测量进行重建，赋相应的能量值，如果不存在，按照称重平均，无能量值
	if(m_TDHead) NData=m_TDHead->NuDataUsed;
	for(i=0;i<MaxNuDataInTriDNode&&i<NData;i++)
	{
		if(m_TDHead)
		{
			m_DensityOfWD->VARIABLESData[i]=m_TDHead->VARIABLESData[i];
			m_DensityOfWD->VARIABLES[i]=m_TDHead->VARIABLES[i];
		}
		else
		{
			m_DensityOfWD->VARIABLESData[i]=0;
			m_DensityOfWD->VARIABLES[i]=_T("Weighing");
		}
	}
	m_DensityOfWD->TitleXYZ[0]=_T("Angle(°)");//x坐标名称;
	m_DensityOfWD->TitleXYZ[1]=_T("Radius(cm)");//y坐标名称;
	m_DensityOfWD->TitleXYZ[2]=_T("Height(cm)");//z坐标名称;
	m_DensityOfWD->NuDataUsed=NData;//数据空间使用数据的个数;	
	m_DensityOfWD->IsScientific=true;
	double MeanDensity=DecSetupData.WeightOfDrum*1000/(PI*pow(DecSetupData.GWDD*0.5,2.0)*DecSetupData.GWDH);//单位：g/cm3
	TriDNodeDefine *dataz=m_DensityOfWD->FirstNode,*datax,*datay;
	//赋坐标值
	for(i=0;i<NuZ&&dataz;i++,dataz=dataz->NextZ)
	{
		for(j=0,datay=dataz;j<NuY&&datay;j++,datay=datay->NextY)
		{
			for(k=0,datax=datay;k<NuX&&datax;k++,datax=datax->NextX)
			{
				datax->x=360.0*k/NuX;
				datax->y=DecSetupData.GWDD*(j+0.5)/NuY;
				datax->z=DecSetupData.GWDH*(i+0.5)/NuZ;
				for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<NData;ndata++)
					datax->data[ndata]=MeanDensity;	//平均密度
				r2=DecSetupData.GWDD*(j+1)/NuY;
				r1=DecSetupData.GWDD*j/NuY;
				datax->s=1.0/NuX*PI*(pow(r2,2)-pow(r1,2))*DecSetupData.GWDH/NuZ;
			}		
		}
	}
	return true;
}
//创建存储重建活度结果的指针空间
bool RSCalculation::CreatActivityHead(int NuX,int NuY,int NuZ) 
{
	m_ActivityOfWD=TriDataB.CreatNodes(NuX,NuY,NuZ);
	if(!m_ActivityOfWD) return false;
	m_ActivityOfWD->Title=_T("Results_Activity.dat");//数据空间的名称，可修改;
	int i,j,k;
	double r2,r1;
	for(i=0;i<MaxNuDataInTriDNode&&i<m_EDHead->NuDataUsed;i++)
	{
		m_ActivityOfWD->VARIABLESData[i]=m_EDHead->VARIABLESData[i];//节点数据空间的横坐标，如果读数据说明是数据，自动识别成double类型;
		m_ActivityOfWD->VARIABLES[i]=m_EDHead->VARIABLES[i];//节点数据说明，随数据输出;
	}
	m_ActivityOfWD->TitleXYZ[0]=_T("Angle(°)");//x坐标名称;
	m_ActivityOfWD->TitleXYZ[1]=_T("Radius(cm)");//y坐标名称;
	m_ActivityOfWD->TitleXYZ[2]=_T("Height(cm)");//z坐标名称;
	m_ActivityOfWD->NuDataUsed=m_EDHead->NuDataUsed;//数据空间使用数据的个数;	
	m_ActivityOfWD->IsScientific=false;
	m_ActivityOfWD->NuDecimal=2;
	m_ActivityOfWD->IsOutputInf=false;
	TriDNodeDefine *dataz=m_ActivityOfWD->FirstNode,*datax,*datay;
	//赋坐标值
	for(i=0;i<NuZ&&dataz;i++,dataz=dataz->NextZ)
		for(j=0,datay=dataz;j<NuY&&datay;j++,datay=datay->NextY)
			for(k=0,datax=datay;k<NuX&&datax;k++,datax=datax->NextX)
			{
				datax->x=360.0*k/NuX;
				datax->y=0.5*DecSetupData.GWDD*(j+0.5)/NuY;
				datax->z=DecSetupData.GWDH*(i+0.5)/NuZ;
				datax->Inf=_T("//Bq");
				r2=DecSetupData.GWDD*(j+1)/NuY;
				r1=DecSetupData.GWDD*j/NuY;
				datax->s=1.0/NuX*PI*(pow(r2,2)-pow(r1,2))*DecSetupData.GWDH/NuZ;
			}		
	return true;
}

//SGS透射重建
bool RSCalculation::TransmissionReconstruction() 
{
	if(DecSetupData.nDetectionType==3) 
	{
		TransmissionReconstructionSTGS();
		return true;
	}
	//创建透射重建的密度矩阵
	if(!m_DensityOfWD) return false;
	//计算投射长度
	CartesianCoordinate PA;
	PA.x=-300;
	PA.y=DecSetupData.GHPGDY;
	PA.z=DecSetupData.GHPGDZ;
	double L=AttenuationLength(PA,DecSetupData.GHPGDY,DecSetupData.GHPGDZ,0.12);
	//进行重建
	TriDNodeDefine *dataTD=m_TDHead->FirstNode;//透射测量文件
	TriDNodeDefine *data=m_DensityOfWD->FirstNode;//密度
	double u,udd;
	int i,ndata;
	for(i=0;i<m_TDHead->NuNodesInZ&&data&&dataTD;i++)
	{
		for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<m_TDHead->NuDataUsed;ndata++)
		{
			u=-1*log(dataTD->data[ndata])/L;
			udd=CAttenuationPerMass(m_TDHead->VARIABLESData[ndata]);
			data->data[ndata]=u/udd;
		}	
		data=data->NextZ;
		dataTD=dataTD->NextZ;
	}
	//输出
	return true;
}
//STGS透射重建,求各环密度，再计算平均密度
bool RSCalculation::TransmissionReconstructionSTGS() 
{
	//创建透射重建的密度矩阵
	if(!m_DensityOfWD) return false;
	//计算投射长度
	CartesianCoordinate PA;
	PA.x=-300;
	PA.y=DecSetupData.GHPGDY;
	PA.z=DecSetupData.GHPGDZ;

	//进行重建
	int i,j,k,ndata;
	TriDNodeHeadDefine *HeadA=TriDataB.CreatNodes(DecSetupData.NuGridsY,m_TDHead->NuNodesInY,1);
	TriDNodeHeadDefine *HeadB=TriDataB.CreatNodes(1,m_TDHead->NuNodesInY,1);
	TriDNodeHeadDefine *HeadX=TriDataB.CreatNodes(1,DecSetupData.NuGridsY,1);
	if(!HeadA||!HeadB||!HeadX) return false;
	HeadA->Title=_T("A.dat");//数据空间的名称，可修改;
	HeadB->Title=_T("B.dat");//数据空间的名称，可修改;
	HeadX->Title=_T("X.dat");//数据空间的名称，可修改;
	for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<m_TDHead->NuDataUsed;ndata++)
	{
		HeadA->VARIABLESData[ndata]=m_TDHead->VARIABLESData[ndata];//节点数据空间的横坐标，如果读数据说明是数据，自动识别成double类型;
		HeadA->VARIABLES[ndata]=m_TDHead->VARIABLES[ndata];//节点数据说明，随数据输出;
		HeadB->VARIABLESData[ndata]=m_TDHead->VARIABLESData[ndata];//节点数据空间的横坐标，如果读数据说明是数据，自动识别成double类型;
		HeadB->VARIABLES[ndata]=m_TDHead->VARIABLES[ndata];//节点数据说明，随数据输出;
		HeadX->VARIABLESData[ndata]=m_TDHead->VARIABLESData[ndata];//节点数据空间的横坐标，如果读数据说明是数据，自动识别成double类型;
		HeadX->VARIABLES[ndata]=m_TDHead->VARIABLES[ndata];//节点数据说明，随数据输出;
	}
	HeadA->NuDataUsed=m_TDHead->NuDataUsed;//数据空间使用数据的个数;
	HeadB->NuDataUsed=m_TDHead->NuDataUsed;//数据空间使用数据的个数;
	HeadX->NuDataUsed=m_TDHead->NuDataUsed;//数据空间使用数据的个数;	


	TriDNodeDefine *dataTD=m_TDHead->FirstNode,*dataTDY;//透射测量文件
	TriDNodeDefine *data=m_DensityOfWD->FirstNode;//密度
	TriDNodeDefine *Bdata,*Xdata;
	TriDNodeDefine *AdataX,*AdataY;
	double DCircle,L,LOld,S[10];
	double udd,density;
	for(k=0;k<DecSetupData.NuGridsY;k++)
	{
		if(!DecSetupData.IsEqualVolumeOfGrids)
			S[k]=(pow(k+1.0,2)-pow(k*1.0,2))/pow(DecSetupData.NuGridsY*1.0,2.0);
		else S[k]=1.0/DecSetupData.NuGridsY;
	}
	for(i=0;i<m_TDHead->NuNodesInZ;i++)
	{
		dataTDY=dataTD;
		Bdata=HeadB->FirstNode;
		AdataY=HeadA->FirstNode;
		for(j=0;j<m_TDHead->NuNodesInY;j++)
		{ //j
			AdataX=AdataY;
			LOld=0;
			//系数矩阵
			for(k=0;k<DecSetupData.NuGridsY;k++)
			{
				PA.y=DecSetupData.GHPGDY+DecSetupData.DSLY*j;
				if(!DecSetupData.IsEqualVolumeOfGrids) DCircle=DecSetupData.GWDD*(k+1.0)/DecSetupData.NuGridsY;//均匀半径;
				else DCircle=DecSetupData.GWDD*sqrt((k+1.0)/DecSetupData.NuGridsY);//均匀面积;
				L=AttenuationLength2(PA,-1*PA.y,PA.z,DCircle);	
				AdataX->data[0]=L-LOld;				
				LOld=L;
				for(ndata=1;ndata<MaxNuDataInTriDNode&&ndata<m_TDHead->NuDataUsed;ndata++)
					AdataX->data[ndata]=AdataX->data[0];
				AdataX=AdataX->NextX;			
			}
			//B值的确定
			for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<m_TDHead->NuDataUsed;ndata++)
			{
				if(dataTDY->data[ndata]>0.00000001) 
				Bdata->data[ndata]=-1*log(dataTDY->data[ndata]);
				else Bdata->data[ndata]=18.42;
			}
			Bdata=Bdata->NextY;
			AdataY=AdataY->NextY;
			dataTDY=dataTDY->NextY;
		}//j
		//x清零
		Xdata=HeadX->FirstNode;
		for(j=0;j<DecSetupData.NuGridsY;j++)
		{
			for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<m_TDHead->NuDataUsed;ndata++)
				Xdata->data[ndata]=0;
			Xdata=Xdata->NextY;
		}
		//透射迭代
		double Cancha[10];
		TriDataB.InitializationX(HeadA,_T("YX"),HeadB,_T("Y"),HeadX,_T("Y"));
		for(int noIteration=0;noIteration<50;noIteration++)
		TriDataB.MLEM_Iteration(HeadA,_T("YX"),HeadB,_T("Y"),HeadX,_T("Y"),Cancha);
		//求平均密度
		for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<m_TDHead->NuDataUsed;ndata++)
		{
			udd=CAttenuationPerMass(m_TDHead->VARIABLESData[ndata]);
			Xdata=HeadX->FirstNode;
			for(j=0;j<DecSetupData.NuGridsY;j++)
			{
				density=Xdata->data[ndata]/udd;
				if(j==0) data->data[ndata]=0;
				data->data[ndata]+=density*S[j];
				Xdata=Xdata->NextY;
			}
		}
		data=data->NextZ;
		dataTD=dataTD->NextZ;
	}
	TriDataB.outputdata_3d(HeadA,m_Folder+HeadA->Title,_T("xyz"));
	TriDataB.outputdata_3d(HeadX,m_Folder+HeadX->Title,_T("xyz"));
	TriDataB.DeleteNodes(HeadA);
	TriDataB.DeleteNodes(HeadB);
	TriDataB.DeleteNodes(HeadX);
	//输出
	return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
//SGS、STGS发射重建
////////////////////////////////////////////////////////////////////////////////////////////////////
bool RSCalculation::EmissionReconstruction() 
{
	if(DecSetupData.nDetectionType==2)
		return (EmissionReconstructionISGS());//ISGS

	//发射探测文件
	if(!m_EDHead) return false;
	//创建发射重建的结果存储指针
	int NuZ=DecSetupData.NuGridsX*DecSetupData.NuGridsY*DecSetupData.NuGridsZ;
	XHead=TriDataB.CreatNodes(1,1,NuZ);
	if(!XHead) return false;
	XHead->Title=_T("X.dat");//数据空间的名称，可修改;
	int i,j,k,ndata;
	for(i=0;i<MaxNuDataInTriDNode&&i<m_EDHead->NuDataUsed;i++)
	{
		XHead->VARIABLESData[i]=m_EDHead->VARIABLESData[i];//节点数据空间的横坐标，如果读数据说明是数据，自动识别成double类型;
		XHead->VARIABLES[i]=m_EDHead->VARIABLES[i];//节点数据说明，随数据输出;
	}
	XHead->TitleXYZ[0]=_T("");//x坐标名称;
	XHead->TitleXYZ[1]=_T("");//y坐标名称;
	XHead->TitleXYZ[2]=_T("Height(cm)");//z坐标名称;
	XHead->NuDataUsed=m_EDHead->NuDataUsed;//数据空间使用数据的个数;	
	XHead->IsScientific=false;
	XHead->NuDecimal=0;
	XHead->IsOutputInf=true;
	TriDNodeDefine *data=XHead->FirstNode;
	//创建系数矩阵
	CoefArrayHead=CreatCoefArrayHead();
	if(!CoefArrayHead) return false;
	//创建残差存储指针
	EMResidHead=TriDataB.CreatNodes(DecSetupData.MaxIteration,1,1);
	if(!EMResidHead) return false;
	EMResidHead->Title=_T("迭代残差.txt");
	EMResidHead->TitleXYZ[0]=_T("NoIteration");
	EMResidHead->TitleXYZ[1]=_T("None");
	EMResidHead->TitleXYZ[2]=_T("None");
	for(i=0;i<m_EDHead->NuDataUsed&&i<MaxNuDataInTriDNode;i++)
	{
		EMResidHead->VARIABLES[i]=m_EDHead->VARIABLES[i];
		EMResidHead->VARIABLESData[i]=m_EDHead->VARIABLESData[i];
	}
	EMResidHead->NuDataUsed=m_EDHead->NuDataUsed;
	EMResidHead->IsScientific=true;
	//进行重建
	TriDNodeHeadDefine *m_EDHeadTemp=TriDataB.CreatNodes(1,1,m_EDHead->NuNodesInX*m_EDHead->NuNodesInY*m_EDHead->NuNodesInZ);
	TriDNodeDefine *datax,*datay,*dataz,*datatemp=m_EDHeadTemp->FirstNode;

	for(j=0,datay=m_EDHead->FirstNode;j<m_EDHead->NuNodesInY&&datay;j++,datay=datay->NextY)
		for(k=0,dataz=datay;k<m_EDHead->NuNodesInZ&&dataz;k++,dataz=dataz->NextZ)

			for(i=0,datax=dataz;i<m_EDHead->NuNodesInX&&datax;i++,datax=datax->NextX)
			{
				for(ndata=0;ndata<m_EDHead->NuDataUsed&&ndata<MaxNuDataInTriDNode;ndata++)
					datatemp->data[ndata]=datax->data[ndata];
				datatemp=datatemp->NextZ;
			}
			m_EDHeadTemp->NuDataUsed=m_EDHead->NuDataUsed;
	TriDataB.outputdata_3d(m_EDHeadTemp,m_Folder+"B.dat",_T("xyz"));///////////////////////////////////////////////////////调使用
	//初始化
	TriDataB.MultiplyNodesData(CoefArrayHead,100000);
	TriDataB.MultiplyNodesData(m_EDHeadTemp,100000);

	if(!TriDataB.InitializationX(CoefArrayHead,_T("ZY"),m_EDHeadTemp,_T("Z"),XHead,_T("Z"))) return false;
	bool key=true;
	//	AfxMessageBox(_T("1"));
	TriDNodeDefine *EMResidNodeOld=NULL,*EMResidNode=EMResidHead->FirstNode;

//TriDataB.LST_Changed(CoefArrayHead,_T("ZY"),m_EDHeadTemp,_T("Z"));

	for(int noIteration=0;noIteration<DecSetupData.MaxIteration;noIteration++)
	{
		EMResidNode->x=noIteration+1;
		if(DecSetupData.TypeOfIteration==1) key=TriDataB.MLEM_Iteration(CoefArrayHead,_T("ZY"),m_EDHeadTemp,_T("Z"),XHead,_T("Z"),EMResidNode->data);
		else if(DecSetupData.TypeOfIteration==2) key=TriDataB.ART_Iteration(CoefArrayHead,_T("ZY"),m_EDHeadTemp,_T("Z"),XHead,_T("Z"),EMResidNode->data);
		else if(DecSetupData.TypeOfIteration==3) key=TriDataB.Gauss_Iteration(CoefArrayHead,_T("ZY"),m_EDHeadTemp,_T("Z"),XHead,_T("Z"),EMResidNode->data);
		else if(DecSetupData.TypeOfIteration==4) key=TriDataB.IterationSGS(CoefArrayHead,_T("ZY"),m_EDHeadTemp,_T("Z"),XHead,_T("Z"),EMResidNode->data);	
		if(!key) return false;
		EMResidNode=EMResidNode->NextX;
	}
	//输出
	SetdataActivityHead(XHead);
	return true;
}

//创建系数矩阵,2维矩阵，I为1，J为列数，即未知数个数，K为行数，即探测次数
TriDNodeHeadDefine *RSCalculation::CreatCoefArrayHead()
{
	if(!m_EDHead) return NULL;
	if(!m_DensityOfWD) return NULL;
	if(!m_ActivityOfWD) return NULL;
	//生成系数矩阵的指针
	int NuZ=m_EDHead->NuNodesInX*m_EDHead->NuNodesInY*m_EDHead->NuNodesInZ;//行数，即探测数
	int NuY=DecSetupData.NuGridsX*DecSetupData.NuGridsY*DecSetupData.NuGridsZ;//列数，即未知数个数
	CoefArrayHead=TriDataB.CreatNodes(1,NuY,NuZ);
	if(!CoefArrayHead) return NULL;
	int i,j,ndata;
	for(i=0;i<m_EDHead->NuDataUsed&&i<MaxNuDataInTriDNode;i++)
	{
		CoefArrayHead->VARIABLES[i]=m_EDHead->VARIABLES[i];
		CoefArrayHead->VARIABLESData[i]=m_EDHead->VARIABLESData[i];
	}
	CoefArrayHead->IsScientific=true;
	CoefArrayHead->NuDecimal=4;
	CoefArrayHead->NuDataUsed=m_EDHead->NuDataUsed;
	CoefArrayHead->Title=_T("A_SGS.txt");
	CoefArrayHead->TitleXYZ[0]=_T("None");
	CoefArrayHead->TitleXYZ[1]=_T("Row");
	CoefArrayHead->TitleXYZ[2]=_T("Line");
	//创建计算
	TriDNodeHeadDefine *DEBody;
	int NoPos=0;
	int nx,ny,nz,nl;
	double x1[20],y1[20],x2[1],y2[1];
	double DEBVersusR[10][100];
	double DEBofOneEquation[2000];
	TriDNodeDefine *dataDEBx,*dataDEBy,*dataDEBz;//体源（网格）
	TriDNodeDefine *dataDEDy=m_EDHead->FirstNode,*dataDEDz;//探测
	TriDNodeDefine *dataDensity;//密度
	TriDNodeDefine *dataCoefArrayy,*dataCoefArrayz=CoefArrayHead->FirstNode;
	CString Ctemp;


	for(i=0;i<m_EDHead->NuNodesInY&&dataDEDy;i++,dataDEDy=dataDEDy->NextY)//测量的偏心位置
	{
		NoPos=int(floor(dataDEDy->y/3.5+0.5));//偏心位置
		//Ctemp.Format(_T("%d"),NoPos);
		//if(key++<4)AfxMessageBox(Ctemp);
		//生成体源效率指针，按照网格计算（IJK分别为：密度、半径、层高），已经根据探测结果，插值生成了对应能量的效率。对于不同偏心位置，要计算偏心序号，重新生成一次
		DEBody=CreatDEFilesB.DE_ProcessingEInCBody(NoPos,DecSetupData.NuGridsY,DecSetupData.IsEqualVolumeOfGrids,m_EDHead->VARIABLESData,m_EDHead->NuDataUsed);
		if(!DEBody) return NULL;
		dataDensity=m_DensityOfWD->FirstNode;//各层密度
		for(j=0,dataDEDz=dataDEDy;j<m_EDHead->NuNodesInZ&&dataDEDz;j++,dataDEDz=dataDEDz->NextZ)//测量的升降位置
		{
			for(ndata=0;ndata<CoefArrayHead->NuDataUsed&&ndata<MaxNuDataInTriDNode;ndata++)//不同能量
			{
				dataDEBz=DEBody->FirstNode;//划分的体源效率		
				for(nz=0;nz<DEBody->NuNodesInZ&&dataDEBz;nz++,dataDEBz=dataDEBz->NextZ)//偏层的网格位置，分：正对层，偏一层，偏二层
					for(ny=0,dataDEBy=dataDEBz;ny<DEBody->NuNodesInY&&dataDEBy;ny++,dataDEBy=dataDEBy->NextY)//不同半径的体源（网格），半径从小到大排列
					{///
						for(nx=0,dataDEBx=dataDEBy;nx<DEBody->NuNodesInX&&dataDEBx;nx++,dataDEBx=dataDEBx->NextX)//不同密度	
						{
							x1[nx]=dataDEBx->x;
							y1[nx]=dataDEBx->data[ndata];
						}
						x2[0]=dataDensity->data[0];
						//插值
						TriDataB.SplineInterpolation(x1, y1, DEBody->NuNodesInX, x2, y2,1,true);
						DEBVersusR[nz][ny]=y2[0];//插值后的一个体源的效率暂存
					}///
				for(nl=0;nl<2000;nl++)
					DEBofOneEquation[nl]=0;//所有值都先定为零
				for(nz=0;nz<DEBody->NuNodesInZ;nz++)//偏层的网格位置，分：正对层，偏一层，偏二层
					for(ny=0;ny<DEBody->NuNodesInY;ny++)//不同半径的体源（网格），半径从小到大排列
					{
						nl=DecSetupData.NuGridsY*j;
						DEBofOneEquation[nl+nz*DecSetupData.NuGridsY+ny]=DEBVersusR[nz][ny];
						if(nz>0&&(nl-nz*DecSetupData.NuGridsY+ny)>=0)
						DEBofOneEquation[nl-nz*DecSetupData.NuGridsY+ny]=DEBVersusR[nz][ny];				
					}
				dataCoefArrayy=dataCoefArrayz;
				for(nl=0;nl<2000&&dataCoefArrayy;nl++,dataCoefArrayy=dataCoefArrayy->NextY)
					dataCoefArrayy->data[ndata]=DEBofOneEquation[nl];
			}//能量循环结束			
			dataDensity=dataDensity->NextZ;//下一层高度的探测位置
			dataCoefArrayz=dataCoefArrayz->NextZ;//下一个方程	
		}//下一个升降位置结束	
			TriDataB.DeleteNodes(DEBody);//删除一个偏心位置的体源效率，计算下一个偏心位置
	}//下一个偏心位置结束
//TriDataB.outputdata_3d(CoefArrayHead,m_Folder+CoefArrayHead->Title,_T("xyz"));
	return CoefArrayHead;
}

//创建系数矩阵,2维矩阵，I为1，J为列数，即未知数个数，K为行数，即探测次数
TriDNodeHeadDefine *RSCalculation::CreatCoefArrayHeadZXWG()
{
	if(!m_EDHead) return NULL;
	if(!m_DensityOfWD) return NULL;
	if(!m_ActivityOfWD) return NULL;
	//生成系数矩阵的指针
	int NuZXWGInR=8;
	int NuZ=m_EDHead->NuNodesInX*m_EDHead->NuNodesInY*m_EDHead->NuNodesInZ;//行数，即探测数
	int NuY=DecSetupData.NuGridsX*NuZXWGInR*DecSetupData.NuGridsZ;//列数，即未知数个数
	CoefArrayHead3=TriDataB.CreatNodes(1,NuY,NuZ);
	if(!CoefArrayHead3) return NULL;
	int i,j,ndata;
	for(i=0;i<m_EDHead->NuDataUsed&&i<MaxNuDataInTriDNode;i++)
	{
		CoefArrayHead3->VARIABLES[i]=m_EDHead->VARIABLES[i];
		CoefArrayHead3->VARIABLESData[i]=m_EDHead->VARIABLESData[i];
	}
	CoefArrayHead3->IsScientific=true;
	CoefArrayHead3->NuDecimal=4;
	CoefArrayHead3->NuDataUsed=m_EDHead->NuDataUsed;
	CoefArrayHead3->Title=_T("A_SGS.txt");
	CoefArrayHead3->TitleXYZ[0]=_T("None");
	CoefArrayHead3->TitleXYZ[1]=_T("Row");
	CoefArrayHead3->TitleXYZ[2]=_T("Line");
	//创建计算
	TriDNodeHeadDefine *DEBody;
	int NoPos=0;
	int nx,ny,nz,nl;
	double x1[20],y1[20],x2[1],y2[1];
	double DEBVersusR[10][100];
	double DEBofOneEquation[2000];
	TriDNodeDefine *dataDEBx,*dataDEBy,*dataDEBz;//体源（网格）
	TriDNodeDefine *dataDEDy=m_EDHead->FirstNode,*dataDEDz;//探测
	TriDNodeDefine *dataDensity;//密度
	TriDNodeDefine *dataCoefArrayy,*dataCoefArrayz=CoefArrayHead3->FirstNode;
	CString Ctemp;

	for(i=0;i<m_EDHead->NuNodesInY&&dataDEDy;i++,dataDEDy=dataDEDy->NextY)//测量的偏心位置
	{
		NoPos=int(floor(dataDEDy->y/3.5+0.5));//偏心位置
		//Ctemp.Format(_T("%d"),NoPos);
		//if(key++<4)AfxMessageBox(Ctemp);
		//生成体源效率指针，按照网格计算（IJK分别为：密度、半径、层高），已经根据探测结果，插值生成了对应能量的效率。对于不同偏心位置，要计算偏心序号，重新生成一次
		DEBody=CreatDEFilesB.DE_ProcessingEInCBody(NoPos,NuZXWGInR,DecSetupData.IsEqualVolumeOfGrids,m_EDHead->VARIABLESData,m_EDHead->NuDataUsed);
		if(!DEBody) return NULL;
		dataDensity=m_DensityOfWD->FirstNode;//各层密度
		for(j=0,dataDEDz=dataDEDy;j<m_EDHead->NuNodesInZ&&dataDEDz;j++,dataDEDz=dataDEDz->NextZ)//测量的升降位置
		{
			for(ndata=0;ndata<CoefArrayHead3->NuDataUsed&&ndata<MaxNuDataInTriDNode;ndata++)//不同能量
			{
				dataDEBz=DEBody->FirstNode;//划分的体源效率		
				for(nz=0;nz<DEBody->NuNodesInZ&&dataDEBz;nz++,dataDEBz=dataDEBz->NextZ)//偏层的网格位置，分：正对层，偏一层，偏二层
					for(ny=0,dataDEBy=dataDEBz;ny<DEBody->NuNodesInY&&dataDEBy;ny++,dataDEBy=dataDEBy->NextY)//不同半径的体源（网格），半径从小到大排列
					{///
						for(nx=0,dataDEBx=dataDEBy;nx<DEBody->NuNodesInX&&dataDEBx;nx++,dataDEBx=dataDEBx->NextX)//不同密度	
						{
							x1[nx]=dataDEBx->x;
							y1[nx]=dataDEBx->data[ndata];
						}
						x2[0]=dataDensity->data[0];
						//插值
						TriDataB.SplineInterpolation(x1, y1, DEBody->NuNodesInX, x2, y2,1,true);
						DEBVersusR[nz][ny]=y2[0];//插值后的一个体源的效率暂存
					}///
				for(nl=0;nl<2000;nl++)
					DEBofOneEquation[nl]=0;//所有值都先定为零
				for(nz=0;nz<DEBody->NuNodesInZ;nz++)//偏层的网格位置，分：正对层，偏一层，偏二层
					for(ny=0;ny<DEBody->NuNodesInY;ny++)//不同半径的体源（网格），半径从小到大排列
					{
						nl=DecSetupData.NuGridsY*j;
						DEBofOneEquation[nl+nz*DecSetupData.NuGridsY+ny]=DEBVersusR[nz][ny];
						if(nz>0&&(nl-nz*DecSetupData.NuGridsY+ny)>=0)
						DEBofOneEquation[nl-nz*DecSetupData.NuGridsY+ny]=DEBVersusR[nz][ny];				
					}
				dataCoefArrayy=dataCoefArrayz;
				for(nl=0;nl<2000&&dataCoefArrayy;nl++,dataCoefArrayy=dataCoefArrayy->NextY)
					dataCoefArrayy->data[ndata]=DEBofOneEquation[nl];
			}//能量循环结束			
			dataDensity=dataDensity->NextZ;//下一层高度的探测位置
			dataCoefArrayz=dataCoefArrayz->NextZ;//下一个方程	
		}//下一个升降位置结束	
			TriDataB.DeleteNodes(DEBody);//删除一个偏心位置的体源效率，计算下一个偏心位置
	}//下一个偏心位置结束
//TriDataB.outputdata_3d(CoefArrayHead,m_Folder+CoefArrayHead->Title,_T("xyz"));
	return CoefArrayHead;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
//ISGS发射重建
////////////////////////////////////////////////////////////////////////////////////////////////////
bool RSCalculation::EmissionReconstructionISGS()//ISGS
{
	//发射探测文件
	if(!m_EDHead) return false;

	//创建ISGS用的指针空间
	if(!CreatCoefInISGS()) return false;

	//创建残差存储指针
	EMResidHead=TriDataB.CreatNodes(DecSetupData.MaxIteration,2,1);
	if(!EMResidHead) return false;
	EMResidHead->Title=_T("迭代残差.txt");
	EMResidHead->TitleXYZ[0]=_T("NoIteration");
	EMResidHead->TitleXYZ[1]=_T("dY");
	EMResidHead->TitleXYZ[2]=_T("None");
	int i,j,k,ndata;
	for(i=0;i<m_EDHead->NuDataUsed&&i<MaxNuDataInTriDNode;i++)
	{
		EMResidHead->VARIABLES[i]=m_EDHead->VARIABLES[i];
		EMResidHead->VARIABLESData[i]=m_EDHead->VARIABLESData[i];
	}
	EMResidHead->NuDataUsed=m_EDHead->NuDataUsed;
	EMResidHead->IsScientific=true;
	//进行重建
	TriDNodeHeadDefine *m_EDHeadTemp1=TriDataB.CreatNodes(1,1,m_EDHead->NuNodesInZ);
	TriDNodeHeadDefine *m_EDHeadTemp2=TriDataB.CreatNodes(1,1,m_EDHead->NuNodesInZ);
	TriDNodeDefine *data1=m_EDHeadTemp1->FirstNode;
	TriDNodeDefine *data2=m_EDHeadTemp2->FirstNode;
	TriDNodeDefine *datatemp=m_EDHead->FirstNode;

	for(k=0;k<m_EDHead->NuNodesInZ&&data1&&data2;k++)
	{
		for(ndata=0;ndata<m_EDHead->NuDataUsed&&ndata<MaxNuDataInTriDNode;ndata++)
		{
			data1->data[ndata]=datatemp->data[ndata];
			data2->data[ndata]=datatemp->NextY->data[ndata];		
		}
		data1=data1->NextZ;
		data2=data2->NextZ;	
		datatemp=datatemp->NextZ;	
	}
	m_EDHeadTemp1->NuDataUsed=m_EDHead->NuDataUsed;
	m_EDHeadTemp2->NuDataUsed=m_EDHead->NuDataUsed;
	bool key;
	int nstep=0;
CString CTemp;
	ofstream output_data;
	output_data.open(m_Folder+_T("123.txt"),ios::trunc);

	TriDNodeDefine *EMResidNodeOld=NULL,*EMResidNode=EMResidHead->FirstNode;
	//////
	{
	for(i=0;i<20;i++)
		for(ndata=0;ndata<20;ndata++)
		{
			rInTime[i][ndata]=0;
			for(j=0;j<28;j++)
				NurInTime[i][ndata][j]=0;
		}
	}
int NuCircle=200;
	for(int nCircle=0,nInTime=0;nCircle<NuCircle;nCircle++,nInTime++)
	{
		//根据等效半径生成系数矩阵，初始半径为零
		if(!CreatCoefArrayHeadISGS()) return false;
		//CTemp.Format(_T("%02d"),nCircle);
	//	TriDataB.outputdata_2d2(CoefArrayHead,m_Folder+CTemp+_T("A.txt"),_T("YZ"),0);
		//TriDataB.outputdata_2d2(CoefArrayHead2,m_Folder+CTemp+_T("A2.txt"),_T("YZ"),0);
		//初始化
		TriDataB.SetNodesData(XHead,0);
		TriDataB.SetNodesData(XHead2,0);
		if(!TriDataB.InitializationX(CoefArrayHead,_T("ZY"),m_EDHeadTemp1,_T("Z"),XHead,_T("Z"))) return false;
		if(!TriDataB.InitializationX(CoefArrayHead2,_T("ZY"),m_EDHeadTemp2,_T("Z"),XHead2,_T("Z"))) return false;
		for(int noIteration=0;noIteration<int(DecSetupData.MaxIteration*1.0/NuCircle);noIteration++,nstep++)
		{
			EMResidNode->x=nstep+1;
			if(DecSetupData.TypeOfIteration==1)
				key=TriDataB.MLEM_Iteration(CoefArrayHead,_T("ZY"),m_EDHeadTemp1,_T("Z"),XHead,_T("Z"),EMResidNode->data);
			else if(DecSetupData.TypeOfIteration==2)
				key=TriDataB.ART_Iteration(CoefArrayHead,_T("ZY"),m_EDHeadTemp1,_T("Z"),XHead,_T("Z"),EMResidNode->data);
		if(!key) return false;
			if(DecSetupData.TypeOfIteration==1)
				key=TriDataB.MLEM_Iteration(CoefArrayHead2,_T("ZY"),m_EDHeadTemp2,_T("Z"),XHead2,_T("Z"),EMResidNode->NextY->data);
			else if(DecSetupData.TypeOfIteration==2)
				key=TriDataB.ART_Iteration(CoefArrayHead2,_T("ZY"),m_EDHeadTemp2,_T("Z"),XHead2,_T("Z"),EMResidNode->NextY->data);
		if(!key) return false;
		EMResidNode=EMResidNode->NextX;
		}


double x1=XHead->FirstNode->NextZ->NextZ->NextZ->NextZ->data[1];
double x2=XHead2->FirstNode->NextZ->NextZ->NextZ->NextZ->data[1];
double a1=CoefArrayHead->FirstNode->NextZ->NextY->NextZ->NextY->NextZ->NextY->NextZ->NextY->data[1];
double a2=CoefArrayHead2->FirstNode->NextZ->NextY->NextZ->NextY->NextZ->NextY->NextZ->NextY->data[1];
double r1=RadiusOfSourceHead->FirstNode->NextZ->NextZ->NextZ->NextZ->data[0];
double r2=RadiusOfSourceHead->FirstNode->NextZ->NextZ->NextZ->NextZ->data[1];
double r3=RadiusOfSourceHead->FirstNode->NextZ->NextZ->NextZ->NextZ->data[2];
	//	output_data<<nstep<<" "<<x1<<" "<<x2<<" "<<a1<<" "<<a2<<" "<<x2*a2/(x1*a1)<<" "<<r<<endl;
		output_data<<nstep<<" "<<r1<<" "<<r2<<" "<<r3<<endl;

		if(!IdentifyRadiusInISGS()) return false;
	}
	output_data.close();
	//输出
	SetdataActivityHead(XHead);
	TriDataB.outputdata_3d(RadiusOfSourceHead,m_Folder+RadiusOfSourceHead->Title,_T("xyz"));
	return true;
}

//创建ISGS用的指针空间
bool RSCalculation::CreatCoefInISGS()
{
	if(!m_EDHead) return false;
	if(!m_DensityOfWD) return false;
//	if(!m_ActivityOfWD) return false;
	//生成系数矩阵的指针
	CoefArrayHead=TriDataB.CreatNodes(1,DecSetupData.NuGridsZ,m_EDHead->NuNodesInZ);//dy=0
	CoefArrayHead2=TriDataB.CreatNodes(1,DecSetupData.NuGridsZ,m_EDHead->NuNodesInZ);//dy>0
	if(!CoefArrayHead||!CoefArrayHead2) return false;
	int i,j,k,ndata;
	for(i=0;i<m_EDHead->NuDataUsed&&i<MaxNuDataInTriDNode;i++)
	{
		CoefArrayHead->VARIABLES[i]=m_EDHead->VARIABLES[i];
		CoefArrayHead->VARIABLESData[i]=m_EDHead->VARIABLESData[i];
		CoefArrayHead2->VARIABLES[i]=m_EDHead->VARIABLES[i];
		CoefArrayHead2->VARIABLESData[i]=m_EDHead->VARIABLESData[i];
	}
	CoefArrayHead->IsScientific=true;
	CoefArrayHead2->IsScientific=true;
	CoefArrayHead->NuDecimal=4;
	CoefArrayHead2->NuDecimal=4;
	CoefArrayHead->NuDataUsed=m_EDHead->NuDataUsed;
	CoefArrayHead2->NuDataUsed=m_EDHead->NuDataUsed;
	CoefArrayHead->Title=_T("A_SGS.txt");
	CoefArrayHead2->Title=_T("A2_SGS.txt");
	CoefArrayHead->TitleXYZ[0]=_T("None");
	CoefArrayHead->TitleXYZ[1]=_T("Row");
	CoefArrayHead->TitleXYZ[2]=_T("Line");
	CoefArrayHead2->TitleXYZ[0]=_T("None");
	CoefArrayHead2->TitleXYZ[1]=_T("Row");
	CoefArrayHead2->TitleXYZ[2]=_T("Line");

	//创建发射重建的结果存储指针
	XHead=TriDataB.CreatNodes(1,1,DecSetupData.NuGridsZ);
	XHead2=TriDataB.CreatNodes(1,1,DecSetupData.NuGridsZ);
	if(!XHead||!XHead2) return false;
	XHead->Title=_T("X.dat");//数据空间的名称，可修改;
	XHead2->Title=_T("X2.dat");//数据空间的名称，可修改;
	for(i=0;i<MaxNuDataInTriDNode&&i<m_EDHead->NuDataUsed;i++)
	{
		XHead->VARIABLESData[i]=m_EDHead->VARIABLESData[i];//节点数据空间的横坐标，如果读数据说明是数据，自动识别成double类型;
		XHead->VARIABLES[i]=m_EDHead->VARIABLES[i];//节点数据说明，随数据输出;
		XHead2->VARIABLESData[i]=m_EDHead->VARIABLESData[i];//节点数据空间的横坐标，如果读数据说明是数据，自动识别成double类型;
		XHead2->VARIABLES[i]=m_EDHead->VARIABLES[i];//节点数据说明，随数据输出;
	}
	XHead->TitleXYZ[0]=_T("None");//x坐标名称;
	XHead->TitleXYZ[1]=_T("None");//y坐标名称;
	XHead->TitleXYZ[2]=_T("Height(cm)");//z坐标名称;
	XHead2->TitleXYZ[0]=_T("None");//x坐标名称;
	XHead2->TitleXYZ[1]=_T("None");//y坐标名称;
	XHead2->TitleXYZ[2]=_T("Height(cm)");//z坐标名称;
	XHead->NuDataUsed=m_EDHead->NuDataUsed;//数据空间使用数据的个数;
	XHead2->NuDataUsed=m_EDHead->NuDataUsed;//数据空间使用数据的个数;	

	//读取两个探测器位置的环面探测效率，并按照实际的探测能量，插值形成指针
	RS_E_Head=NULL;	RS_E_Head2=NULL;
	TriDNodeHeadDefine *RS_E_HeadTemp;
	TriDNodeDefine *datax,*datay,*dataz,*datax2,*datay2,*dataz2,*DensityData;

	//dy=0的效率插值
	RS_E_HeadTemp=CreatDEFilesB.ReadRSE_ISGS(DecSetupData.GHPGDY);//dy=0;
	if(!RS_E_HeadTemp) return false;
	if(RS_E_HeadTemp->NuNodesInZ!=3) return false;
	RS_E_Head=TriDataB.CreatNodes(m_EDHead->NuNodesInZ,RS_E_HeadTemp->NuNodesInY,RS_E_HeadTemp->NuNodesInZ);
	if(!RS_E_Head) {TriDataB.DeleteNodes(RS_E_HeadTemp);return false;}
	RS_E_Head->Title=_T("RS_E.dat");
	for(i=0;i<MaxNuDataInTriDNode&&i<m_EDHead->NuDataUsed;i++)
	{
		RS_E_Head->VARIABLESData[i]=m_EDHead->VARIABLESData[i];//节点数据空间的横坐标，如果读数据说明是数据，自动识别成double类型;
		RS_E_Head->VARIABLES[i]=m_EDHead->VARIABLES[i];//节点数据说明，随数据输出;
	}
	RS_E_Head->TitleXYZ[0]=_T("Segments");//段层;
	RS_E_Head->TitleXYZ[1]=_T("Radius");//半径;
	RS_E_Head->TitleXYZ[2]=_T("Height");//z坐标名称;
	RS_E_Head->NuDataUsed=m_EDHead->NuDataUsed;//当前层高度，偏一层高度，偏二层高度;
	double x1[100],y1[3][100],x2[1],y2[1],ydY0[20][MaxNuDataInTriDNode],ydY1[20][MaxNuDataInTriDNode],ydY2[20][MaxNuDataInTriDNode];
	int N1;
	//半径方向
	for(j=0,datay=RS_E_HeadTemp->FirstNode,datay2=RS_E_Head->FirstNode;j<RS_E_HeadTemp->NuNodesInY&&datay&&datay2;j++,datay=datay->NextY,datay2=datay2->NextY)
	{
		//能量方向
		for(ndata=0;ndata<RS_E_HeadTemp->NuDataUsed&&ndata<MaxNuDataInTriDNode;ndata++)
		{///
			//密度方向
			for(k=0,datax=datay;k<RS_E_HeadTemp->NuNodesInX&&datax;k++,datax=datax->NextX)	
			{//
				x1[k]=datax->x;
				y1[0][k]=datax->data[ndata];
				if(datax->NextZ) y1[1][k]=datax->NextZ->data[ndata];
				if(datax->NextZ&&datax->NextZ->NextZ) y1[2][k]=datax->NextZ->NextZ->data[ndata];
				N1=k;
			}//
			//段层方向
			DensityData=m_DensityOfWD->FirstNode;
			for(i=0;i<m_EDHead->NuNodesInZ&&DensityData;i++,DensityData=DensityData->NextZ)
			{//
				x2[0]=DensityData->data[0];//当前层的密度
				TriDataB.SplineInterpolation(x1, y1[0], N1, x2, y2, 1, true);
				ydY0[i][ndata]=y2[0];	
				TriDataB.SplineInterpolation(x1, y1[1], N1, x2, y2, 1, true);
				ydY1[i][ndata]=y2[0];				
				TriDataB.SplineInterpolation(x1, y1[2], N1, x2, y2, 1, true);
				ydY2[i][ndata]=y2[0];				
			}//	
		}///
		//能量插值
			for(i=0,datax2=datay2;i<m_EDHead->NuNodesInZ&&datax2;i++,datax2=datax2->NextX)//注意，RS_E_Head的x方向为段层方向
			{
				//当前层
				TriDataB.SplineInterpolation(RS_E_HeadTemp->VARIABLESData,ydY0[i], RS_E_HeadTemp->NuDataUsed, RS_E_Head->VARIABLESData,datax2->data,RS_E_Head->NuDataUsed,true);
				//偏一层	
				TriDataB.SplineInterpolation(RS_E_HeadTemp->VARIABLESData,ydY1[i], RS_E_HeadTemp->NuDataUsed, RS_E_Head->VARIABLESData,datax2->NextZ->data,RS_E_Head->NuDataUsed,true);
				//偏二层	
				TriDataB.SplineInterpolation(RS_E_HeadTemp->VARIABLESData,ydY2[i], RS_E_HeadTemp->NuDataUsed, RS_E_Head->VARIABLESData,datax2->NextZ->NextZ->data,RS_E_Head->NuDataUsed,true);			
			}
		//坐标值
			DensityData=m_DensityOfWD->FirstNode;
			for(i=0,datax2=datay2;i<m_EDHead->NuNodesInZ&&datax2;i++,datax2=datax2->NextX)//注意，RS_E_Head的x方向为段层方向
			{
				for(k=0,dataz=datay,dataz2=datax2;k<RS_E_HeadTemp->NuNodesInZ&&dataz&&dataz2;k++,dataz=dataz->NextZ,dataz2=dataz2->NextZ)
				{
					dataz2->x=DensityData->z;
					dataz2->y=dataz->y;
					dataz2->z=dataz->z;				
				}	
				DensityData=DensityData->NextZ;
			}
	}
	RS_E_Head=TriDataB.ExpendDatabufInY(RS_E_Head,281);
	TriDataB.DeleteNodes(RS_E_HeadTemp); RS_E_HeadTemp=NULL;

	//dy>0偏心位置的效率插值
	RS_E_HeadTemp=CreatDEFilesB.ReadRSE_ISGS(DecSetupData.GHPGDY+DecSetupData.DSLY);//dy>0;
	if(!RS_E_HeadTemp) return false;
	if(RS_E_HeadTemp->NuNodesInZ!=3) return false;
	RS_E_Head2=TriDataB.CreatNodes(m_EDHead->NuNodesInZ,RS_E_HeadTemp->NuNodesInY,RS_E_HeadTemp->NuNodesInZ);
	if(!RS_E_Head2) {TriDataB.DeleteNodes(RS_E_HeadTemp);return false;}
	RS_E_Head2->Title=_T("RS_E2.dat");
	for(i=0;i<MaxNuDataInTriDNode&&i<m_EDHead->NuDataUsed;i++)
	{
		RS_E_Head2->VARIABLESData[i]=m_EDHead->VARIABLESData[i];//节点数据空间的横坐标，如果读数据说明是数据，自动识别成double类型;
		RS_E_Head2->VARIABLES[i]=m_EDHead->VARIABLES[i];//节点数据说明，随数据输出;
	}
	RS_E_Head2->TitleXYZ[0]=_T("Segments");//段层;
	RS_E_Head2->TitleXYZ[1]=_T("Radius");//半径;
	RS_E_Head2->TitleXYZ[2]=_T("Height");//z坐标名称;
	RS_E_Head2->NuDataUsed=m_EDHead->NuDataUsed;//当前层高度，偏一层高度，偏二层高度;
	//半径方向
	for(j=0,datay=RS_E_HeadTemp->FirstNode,datay2=RS_E_Head2->FirstNode;j<RS_E_HeadTemp->NuNodesInY&&datay&&datay2;j++,datay=datay->NextY,datay2=datay2->NextY)
	{
		//能量方向
		for(ndata=0;ndata<RS_E_HeadTemp->NuDataUsed&&ndata<MaxNuDataInTriDNode;ndata++)
		{
			//密度方向
			for(k=0,datax=datay;k<RS_E_HeadTemp->NuNodesInX&&datax;k++,datax=datax->NextX)	
			{
				x1[k]=datax->x;
				y1[0][k]=datax->data[ndata];
				if(datax->NextZ) y1[1][k]=datax->NextZ->data[ndata];
				if(datax->NextZ&&datax->NextZ->NextZ) y1[2][k]=datax->NextZ->NextZ->data[ndata];
				N1=k;
			}
			//段层方向
			DensityData=m_DensityOfWD->FirstNode;
			for(i=0;i<m_EDHead->NuNodesInZ&&DensityData;i++,DensityData=DensityData->NextZ)
			{
				x2[0]=DensityData->data[0];//当前层的密度
				TriDataB.SplineInterpolation(x1, y1[0], N1, x2, y2, 1,true);
				ydY0[i][ndata]=y2[0];	
				TriDataB.SplineInterpolation(x1, y1[1], N1, x2, y2, 1,true);
				ydY1[i][ndata]=y2[0];				
				TriDataB.SplineInterpolation(x1, y1[2], N1, x2, y2, 1,true);
				ydY2[i][ndata]=y2[0];				
			}	
		}
		//能量插值
			for(i=0,datax2=datay2;i<m_EDHead->NuNodesInZ&&datax2;i++,datax2=datax2->NextX)//注意，RS_E_Head的x方向为段层方向
			{
				//当前层
				TriDataB.SplineInterpolation(RS_E_HeadTemp->VARIABLESData,ydY0[i], RS_E_HeadTemp->NuDataUsed, RS_E_Head2->VARIABLESData,datax2->data,RS_E_Head2->NuDataUsed,true);
				//偏一层	
				TriDataB.SplineInterpolation(RS_E_HeadTemp->VARIABLESData,ydY1[i], RS_E_HeadTemp->NuDataUsed, RS_E_Head2->VARIABLESData,datax2->NextZ->data,RS_E_Head2->NuDataUsed,true);
				//偏二层	
				TriDataB.SplineInterpolation(RS_E_HeadTemp->VARIABLESData,ydY2[i], RS_E_HeadTemp->NuDataUsed, RS_E_Head2->VARIABLESData,datax2->NextZ->NextZ->data,RS_E_Head2->NuDataUsed,true);			
			}
			//坐标值
			DensityData=m_DensityOfWD->FirstNode;
			for(i=0,datax2=datay2;i<m_EDHead->NuNodesInZ&&datax2;i++,datax2=datax2->NextX)//注意，RS_E_Head的x方向为段层方向
			{
				for(k=0,dataz=datay,dataz2=datax2;k<RS_E_HeadTemp->NuNodesInZ&&dataz&&dataz2;k++,dataz=dataz->NextZ,dataz2=dataz2->NextZ)
				{
					dataz2->x=DensityData->z;
					dataz2->y=dataz->y;
					dataz2->z=dataz->z;				
				}	
				DensityData=DensityData->NextZ;
			}
	}
	TriDataB.DeleteNodes(RS_E_HeadTemp); RS_E_HeadTemp=NULL;
	RS_E_Head2=TriDataB.ExpendDatabufInY(RS_E_Head2,281);

	//探测效率的比值，每层每个能量都有比值
	RelativeEHead=TriDataB.CreatNodes(1,RS_E_Head->NuNodesInY,RS_E_Head->NuNodesInX);
	if(!RelativeEHead) return false;
	for(i=0;i<RS_E_Head->NuDataUsed&&i<MaxNuDataInTriDNode;i++)
	{
		RelativeEHead->VARIABLES[i]=RS_E_Head->VARIABLES[i];
		RelativeEHead->VARIABLESData[i]=RS_E_Head->VARIABLESData[i];
	}

	RelativeEHead->IsScientific=true;
	RelativeEHead->NuDecimal=8;
	RelativeEHead->NuDataUsed=RS_E_Head->NuDataUsed;
	RelativeEHead->Title=_T("RelativeE_ISGS.txt");
	RelativeEHead->TitleXYZ[0]=_T("None");
	RelativeEHead->TitleXYZ[1]=_T("Radius");//各半径,注意半径是根据指针坐标由大变小的
	RelativeEHead->TitleXYZ[2]=_T("Height");//各断层

	TriDNodeDefine *REdataz,*REdatay;
	dataz=RS_E_Head->FirstNode;
	dataz2=RS_E_Head2->FirstNode;
	REdataz=RelativeEHead->FirstNode;
	for(i=0;i<RS_E_Head->NuNodesInX&&dataz&&dataz2&&REdataz;i++,dataz=dataz->NextX,dataz2=dataz2->NextX,REdataz=REdataz->NextZ)
	for(j=0,datay=dataz,datay2=dataz2,REdatay=REdataz;j<RS_E_Head->NuNodesInY&&datay&&datay2&&REdatay;j++,datay=datay->NextY,datay2=datay2->NextY,REdatay=REdatay->NextY)
	{
		REdatay->y=datay->y;
		REdatay->z=datay->x;
		for(ndata=0;ndata<RS_E_Head->NuDataUsed&&ndata<MaxNuDataInTriDNode;ndata++)		
				REdatay->data[ndata]=datay2->data[ndata]/datay->data[ndata];
	}

	//等效半径
	RadiusOfSourceHead=TriDataB.CreatNodes(1,1,m_EDHead->NuNodesInZ);
	if(!RadiusOfSourceHead) return false;
	for(i=0;i<m_EDHead->NuDataUsed&&i<MaxNuDataInTriDNode;i++)
	{
		RadiusOfSourceHead->VARIABLES[i]=m_EDHead->VARIABLES[i];
		RadiusOfSourceHead->VARIABLESData[i]=m_EDHead->VARIABLESData[i];
	}
	RadiusOfSourceHead->IsScientific=true;
	RadiusOfSourceHead->NuDecimal=8;
	RadiusOfSourceHead->NuDataUsed=m_EDHead->NuDataUsed;
	RadiusOfSourceHead->Title=_T("RadiusOfSource_ISGS.txt");
	RadiusOfSourceHead->TitleXYZ[0]=_T("None");
	RadiusOfSourceHead->TitleXYZ[1]=_T("None");
	RadiusOfSourceHead->TitleXYZ[2]=_T("Radius");
	//TriDataB.SetNodesData(RadiusOfSourceHead,14);

	TriDataB.outputdata_3d(RS_E_Head,m_Folder+RS_E_Head->Title,_T("xyz"));
	TriDataB.outputdata_3d(RS_E_Head2,m_Folder+RS_E_Head2->Title,_T("xyz"));
	TriDataB.outputdata_3d(RelativeEHead,m_Folder+RelativeEHead->Title,_T("xyz"));

	return true;
}

//判断等效半径
bool RSCalculation::IdentifyRadiusInISGS()
{
	if(!XHead||!XHead2||!RadiusOfSourceHead||!RelativeEHead) return false;
	if(!CoefArrayHead||!CoefArrayHead2) return false;
	TriDNodeDefine *datax1=XHead->FirstNode;
	TriDNodeDefine *datax2=XHead2->FirstNode;
	TriDNodeDefine *datar=RadiusOfSourceHead->FirstNode;
	TriDNodeDefine *dataRE=RelativeEHead->FirstNode,*dataREy,*dataREy0,*dataREyEnd;
	TriDNodeDefine *CoefA1=CoefArrayHead->FirstNode;
	TriDNodeDefine *CoefA2=CoefArrayHead2->FirstNode;

	int i,j,ndata;
	double RatioOfX;
	double rtemp=0;
//	int MaxrNo,MaxJ;
	for(i=0;i<XHead->NuNodesInZ&&i<XHead2->NuNodesInZ&&i<RadiusOfSourceHead->NuNodesInZ;i++)
	{////
		for(ndata=0;ndata<XHead->NuDataUsed&&ndata<RadiusOfSourceHead->NuDataUsed&&ndata<MaxNuDataInTriDNode;ndata++)	
		{///
			dataREy0=dataRE;
			if(fabs(datax1->data[ndata])>1e-8)
			{//
				RatioOfX=datax2->data[ndata]*CoefA2->data[ndata]/(datax1->data[ndata]*CoefA1->data[ndata]);
				//原坐标是从28到0，将指针变到末尾，最小值处
				dataREyEnd=dataREy0;
				for(j=0;dataREy0&&dataREy0->NextY;j++)
					dataREy0=dataREy0->NextY;
				//如果r=0的比值小于r=28的比值
				dataREy=dataREy0;
				rtemp=0;
				while(dataREy)
				{
					if(dataREy->LastY==NULL) {rtemp=dataREy->y;break;}
					else if(dataREy->data[ndata]>=dataREy->LastY->data[ndata])//沿LastY递减
					{
						if(RatioOfX>=dataREy->data[ndata]) {rtemp=dataREy->y;break;}
						else if(RatioOfX>=dataREy->LastY->data[ndata]) {rtemp=0.5*(dataREy->y+dataREy->LastY->y);break;}					
					}
					else//沿LastY递增
					{
						if(RatioOfX<=dataREy->data[ndata]) {rtemp=dataREy->y;break;}
						else if(RatioOfX<=dataREy->LastY->data[ndata]) {rtemp=0.5*(dataREy->y+dataREy->LastY->y);break;}				
					}
						dataREy=dataREy->LastY;
				}
			}//
			else rtemp=0;
		//	datar->data[ndata]=rtemp;
		double relax=0.9;
			datar->data[ndata]=datar->data[ndata]*relax+rtemp*(1-relax);
			
		//	rInTime[i][ndata]+=rtemp;		
		//	datar->data[ndata]=rInTime[i][ndata]/(nInTime+1);//等效半径取迭代过程中的平均值
			/*
			NurInTime[i][ndata][int(rtemp)]++;
			for(j=0,MaxrNo=0,MaxJ=0;j<28;j++)
				if(MaxrNo<NurInTime[i][ndata][j]) {MaxrNo=NurInTime[i][ndata][j];MaxJ=j;}
					datar->data[ndata]=MaxJ;
					*/


		}///
		datax1=datax1->NextZ;
		datax2=datax2->NextZ;
		datar=datar->NextZ;
		dataRE=dataRE->NextZ;
		if(CoefA1->NextY) CoefA1=CoefA1->NextY->NextZ;
		if(CoefA2->NextY) CoefA2=CoefA2->NextY->NextZ;
	}////
	return true;
}
//创建ISGS的系数矩阵
bool RSCalculation::CreatCoefArrayHeadISGS()
{
	if(!m_EDHead) return false;
	if(!CoefArrayHead||!CoefArrayHead2) return false;
	if(!RadiusOfSourceHead) return false;
	if(!RS_E_Head||!RS_E_Head2) return false;
	TriDNodeDefine *datal1=CoefArrayHead->FirstNode;
	TriDNodeDefine *datal2=CoefArrayHead2->FirstNode;
	TriDNodeDefine *dataE1=RS_E_Head->FirstNode,*dataE21;
	TriDNodeDefine *dataE2=RS_E_Head2->FirstNode,*dataE22;
	TriDNodeDefine *datar=RadiusOfSourceHead->FirstNode;
	int i,j,line;
	double E11,E12,E10,E21,E22,E20;
	for(line=0;line<m_EDHead->NuNodesInZ&&datal1&&datal2;line++)
	{
		for(i=0;i<m_EDHead->NuDataUsed&&i<MaxNuDataInTriDNode;i++)
		{
			dataE21=dataE1;
			dataE22=dataE2;
			for(j=0;j<RS_E_Head->NuNodesInY&&dataE21&&dataE22;j++,dataE21=dataE21->NextY,dataE22=dataE22->NextY)
			{
				if(datar->data[i]>=dataE21->y) 
				{
					E10=dataE21->data[i];E11=dataE21->NextZ->data[i];E12=dataE21->NextZ->NextZ->data[i];
					E20=dataE22->data[i];E21=dataE22->NextZ->data[i];E22=dataE22->NextZ->NextZ->data[i];					
					break;
				}
				else if(dataE21->NextY==NULL) 
				{
					E10=dataE21->data[i];E11=dataE21->NextZ->data[i];E12=dataE21->NextZ->NextZ->data[i];
					E20=dataE22->data[i];E21=dataE22->NextZ->data[i];E22=dataE22->NextZ->NextZ->data[i];					
					break;
				}
				else if(datar->data[i]<=dataE21->y&&datar->data[i]>=dataE21->NextY->y)
				{
					E10=0.5*(dataE21->data[i]+dataE21->NextY->data[i]);
					E11=0.5*(dataE21->NextZ->data[i]+dataE21->NextY->NextZ->data[i]);
					E12=0.5*(dataE21->NextZ->NextZ->data[i]+dataE21->NextY->NextZ->NextZ->data[i]);
					E20=0.5*(dataE22->data[i]+dataE22->NextY->data[i]);
					E21=0.5*(dataE22->NextZ->data[i]+dataE22->NextY->NextZ->data[i]);
					E22=0.5*(dataE22->NextZ->NextZ->data[i]+dataE22->NextY->NextZ->NextZ->data[i]);
					break;
				}		
			}
			/*
			if(dataE21->NextY==NULL)
			{
				E10=dataE21->data[i];E11=dataE21->NextZ->data[i];E12=dataE21->NextZ->NextZ->data[i];
				E20=dataE22->data[i];E21=dataE22->NextZ->data[i];E22=dataE22->NextZ->NextZ->data[i];			
			}
			else
			{
				k=(datar->data[i]-dataE21->NextY->y)/(dataE21->y-dataE21->NextY->y);			
				E10=dataE21->NextY->data[i]+k*(dataE21->data[i]-dataE21->NextY->data[i]);
				E11=dataE21->NextZ->NextY->data[i]+k*(dataE21->NextZ->data[i]-dataE21->NextZ->NextY->data[i]);
				E12=dataE21->NextZ->NextZ->NextY->data[i]+k*(dataE21->NextZ->NextZ->data[i]-dataE21->NextZ->NextZ->NextY->data[i]);
				E20=dataE22->NextY->data[i]+k*(dataE22->data[i]-dataE22->NextY->data[i]);
				E21=dataE22->NextZ->NextY->data[i]+k*(dataE22->NextZ->data[i]-dataE22->NextZ->NextY->data[i]);
				E22=dataE22->NextZ->NextZ->NextY->data[i]+k*(dataE22->NextZ->NextZ->data[i]-dataE22->NextZ->NextZ->NextY->data[i]);			
			}
			*/
			datal1->data[i]=E10;
			if(datal1->LastY) datal1->LastY->data[i]=E11;
			if(datal1->LastY&&datal1->LastY->LastY) datal1->LastY->LastY->data[i]=E12;		
			if(datal1->NextY) datal1->NextY->data[i]=E11;
			if(datal1->NextY&&datal1->NextY->NextY) datal1->NextY->NextY->data[i]=E12;			
			datal2->data[i]=E20;
			if(datal2->LastY) datal2->LastY->data[i]=E21;
			if(datal2->LastY&&datal2->LastY->LastY) datal1->LastY->LastY->data[i]=E22;		
			if(datal2->NextY) datal2->NextY->data[i]=E21;
			if(datal2->NextY&&datal2->NextY->NextY) datal1->NextY->NextY->data[i]=E22;	

		}
		if(datal1->NextZ) datal1=datal1->NextZ->NextY;
		if(datal2->NextZ) datal2=datal2->NextZ->NextY;	
		datar=datar->NextZ;
		dataE1=dataE1->NextX;
		dataE2=dataE2->NextX;
	}
	return true;
}

//创建存储重建活度结果的指针空间
//注意：
bool RSCalculation::SetdataActivityHead(TriDNodeHeadDefine *headX) 
{
	int i,j,k,ndata;
	TriDNodeDefine *dataz=m_ActivityOfWD->FirstNode,*datay,*datax;
	TriDNodeDefine *dataX=headX->FirstNode;
	//赋坐标值
	for(i=0;i<m_ActivityOfWD->NuNodesInZ&&dataz;i++,dataz=dataz->NextZ)
		for(j=0,datay=dataz;j<m_ActivityOfWD->NuNodesInY&&datay;j++,datay=datay->NextY)
			for(k=0,datax=datay;k<m_ActivityOfWD->NuNodesInX&&datax;k++,datax=datax->NextX)
			{
				for(ndata=0;ndata<m_ActivityOfWD->NuDataUsed&&ndata<MaxNuDataInTriDNode;ndata++)//不同能量
					datax->data[ndata]=dataX->data[ndata];				
				if(dataX) dataX=dataX->NextZ;
			}		
	return true;
}

//输出重建结果到文件
double  RSCalculation::OutputResults()
{
	if(!m_DensityOfWD) return -999;
	if(!m_ActivityOfWD) return -999;
    int i,j;
	int NLine=0;
	CString CTemp,Ctemp;

	///////////////////////////////////////////////////////////////建立结果报告文件
	ofstream output;
	output.open(m_Folder+_T("#DetectionReport.txt"),ios::trunc);
	//获得创建时间
	CTime time;
    time = CTime::GetCurrentTime();
	output<<(LPCTSTR)_T("Report for TEST of ")<<(LPCTSTR)DecSetupData.NoOfWD;
	output<<time.Format(_T("    %Y-%m-%d %H:%M:%S"))<<endl;NLine++;
	CTemp.Format(_T("%s %s"),_T("Developer:"),_T("Shanghai Jiao Tong University")); output<<(LPCTSTR)CTemp<<endl<<endl;NLine++;

	//项目信息
	ClassfyingFoldersInf(m_Folder);//读取时间节点信息
	output<<(LPCTSTR)_T("========== Item Information ==========")<<endl;NLine++;
	CTemp.Format(_T("%-25s%s"),_T("Item ID:"),DecSetupData.NoOfWD); 	output<<(LPCTSTR)CTemp<<endl;NLine++;
	CTemp.Format(_T("%-25s%s"),_T("Detection Time:"),DetectionTime); output<<(LPCTSTR)CTemp<<endl;NLine++;	
	CTemp.Format(_T("%-25s%s"),_T("Reconstruction Time:"),ReconstructionTime); output<<(LPCTSTR)CTemp<<endl;NLine++;	
	CTemp.Format(_T("%-25s%s"),_T("Matrix Type:"),DecSetupData.TypeOfMatrix); output<<(LPCTSTR)CTemp<<endl;	NLine++;	
	CTemp.Format(_T("%-25s%s"),_T("Container Type:"),DecSetupData.TypeOfWD); output<<(LPCTSTR)CTemp<<endl;NLine++;	
	CTemp.Format(_T("%-25s%g"),_T("Weight (kg):"),DecSetupData.WeightOfDrum); output<<(LPCTSTR)CTemp<<endl;NLine++;
	double volume=PI*pow(0.5*DecSetupData.GWDD,2)*DecSetupData.GWDH/1000;
	CTemp.Format(_T("%-25s%.2f"),_T("Container Volume (l):"),volume); output<<(LPCTSTR)CTemp<<endl;NLine++;
	output<<endl;NLine++;
	//探测信息
	output<<(LPCTSTR)_T("========== Detection Information ==========")<<endl;NLine++;
	if(DecSetupData.DetectionType.Find(_T("TGS"))>=0)	CTemp.Format(_T("%-25s%s"),_T("Detection Type:"),_T("Tomographic Gamma Scanning"));	
	else if(DecSetupData.DetectionType.Find(_T("STGS"))>=0)	CTemp.Format(_T("%-25s%s"),_T("Detection Type:"),_T("Semi-Tomographic Gamma Scanning"));
	else if(DecSetupData.DetectionType.Find(_T("ISGS"))>=0)	CTemp.Format(_T("%-25s%s"),_T("Detection Type:"),_T("Improved Segmented Gamma Scanning"));
	else	CTemp.Format(_T("%-25s%s"),_T("Detection Type:"),_T("Segmented Gamma Scanning"));
	output<<(LPCTSTR)CTemp<<endl;NLine++;
	CTemp.Format(_T("%-25s%d\t dH = %.2f cm"),_T("Segments:"),DecSetupData.DSNumZ,DecSetupData.DSLZ); output<<(LPCTSTR)CTemp<<endl;NLine++;
	if(DecSetupData.DSNumY>1)//如果半径划分数大于1
	{CTemp.Format(_T("%-25s%d"),_T("Division versus Radius:"),DecSetupData.DSNumY); output<<(LPCTSTR)CTemp<<endl;NLine++;}
	if(DecSetupData.DetectionType.Find(_T("TGS"))>=0)//如果TGS
	{CTemp.Format(_T("%-25s%d"),_T("Division versus Circle:"),DecSetupData.DSNumA); output<<(LPCTSTR)CTemp<<endl;NLine++;}
	CTemp.Format(_T("%-25sx=%g y=%g z=%g (cm)"),_T("Detector starts from:"),DecSetupData.GHPGDX,DecSetupData.GHPGDY,DecSetupData.GHPGDZ); output<<(LPCTSTR)CTemp<<endl;NLine++;
	if(DecSetupData.IsTransmission)
	{CTemp.Format(_T("%-25s%s (Source %s)"),_T("Transmission Reconstructon:"),_T("True"),DecSetupData.TransmissionSource); output<<(LPCTSTR)CTemp<<endl;NLine++;}	
	else 
	{CTemp.Format(_T("%-25s%s (Density is calculated by weight)"),_T("Transmission Reconstructon:"),_T("False")); output<<(LPCTSTR)CTemp<<endl;NLine++;}	
	if(DecSetupData.TypeOfIteration==1) 
	{CTemp.Format(_T("%-25s%s"),_T("Iteration type:"),_T("ML-EMA")); output<<(LPCTSTR)CTemp<<endl;NLine++;}
	else if(DecSetupData.TypeOfIteration==2)
	{CTemp.Format(_T("%-25s%s"),_T("Iteration type:"),_T("ART")); output<<(LPCTSTR)CTemp<<endl;NLine++;}
	CTemp.Format(_T("%-25s%d"),_T("Maximal iteration number:"),DecSetupData.MaxIteration); output<<(LPCTSTR)CTemp<<endl;NLine++;		
	output<<endl;NLine++;

	//重建结果信息
	output<<(LPCTSTR)_T("========== Reconstruction Results ==========")<<endl;NLine++;
	//统计计算
	Statistics();
	CTemp.Format(_T("%-25s%.2f"),_T("Mean Density (g/cm^3):"),m_MeanDensity); output<<(LPCTSTR)CTemp<<endl;NLine++;
	CTemp.Format(_T("%-25s%.2e"),_T("Total Activity (Bq):"),m_TotalActivity); output<<(LPCTSTR)CTemp<<endl;NLine++;
	CTemp.Format(_T("%-25s"),_T("Activity versus energy:")); output<<(LPCTSTR)CTemp<<endl;NLine++;
	CTemp.Format(_T("| %3s | %10s | %10s | %10s | %14s |"),_T("No."),_T("Nuclide"),_T("Engergy"),_T("Activity"),_T("Concentration"));
	output<<(LPCTSTR)CTemp<<endl;NLine++;
	CTemp.Format(_T("| %3s | %10s | %10s | %10s | %14s |"),_T(" "),_T(" "),_T("MeV"),_T("Bq"),_T("Bq/g"));
	output<<(LPCTSTR)CTemp<<endl;NLine++;
	CString Nucilde;
	for(j=0;j<m_ActivityOfWD->NuDataUsed;j++)
	{
		Nucilde=IdentifyNuclide(m_ActivityOfWD->VARIABLESData[j]);
		CTemp.Format(_T("| %03d | %10s | %10s | %10.2e | %14.2e |"),j+1,Nucilde,m_ActivityOfWD->VARIABLES[j],m_MeanActivity[j],m_MeanActivity[j]/DecSetupData.WeightOfDrum/1000);
		output<<(LPCTSTR)CTemp<<endl;NLine++;		
	}
	output<<endl;NLine++;

	//while(NLine++<32) output<<endl;

	//密度分布
	output<<(LPCTSTR)_T("========== Density Distribution Versus Height ==========")<<endl;
	output<<_T("Density (g/cm^3):")<<endl;
	//output<<_T("   _______________________________________________")<<endl;
	CTemp.Format(_T("|%5s|%8s |%8s |"),_T("Seg."),_T("H (cm)"),_T("Density"));
	output<<(LPCTSTR)CTemp<<endl;
	for(i=0;i<m_DensityOfWD->NuNodesInZ;i++)
	{
		CTemp.Format(_T("| %03d |%8.2f |%8.4f |"),i+1,m_H[i],m_Density[i]);
		output<<(LPCTSTR)CTemp<<endl;//输出各层的平均密度
	}
	//平均密度
	CTemp.Format(_T("| Mean Density  |%8.4f |"),m_MeanDensity);
	output<<(LPCTSTR)CTemp<<endl;
	output<<endl;

	//活度
	output<<(LPCTSTR)_T("========== Activity Distribution Versus Height ==========")<<endl;
	output<<_T("Activity (Bq):")<<endl;
	//发射重建结果输出
	CTemp.Format(_T("|%5s|%8s |"),_T("Seg."),_T("H (cm)"));
	for(j=0;j<m_ActivityOfWD->NuDataUsed;j++)
	{
		Ctemp.Format(_T("%7sMeV |"),m_ActivityOfWD->VARIABLES[j]);
		CTemp+=Ctemp;
	}
		Ctemp.Format(_T("%10s |"),_T("Total"));
		CTemp+=Ctemp;
	output<<(LPCTSTR)CTemp<<endl;

	double SumActi;
	for(i=0;i<m_ActivityOfWD->NuNodesInZ;i++)
	{
		CTemp.Format(_T("| %03d |%8.2f |"),i+1,m_H[i]);
		SumActi=0;
		for(j=0;j<m_ActivityOfWD->NuDataUsed;j++)
		{
			Ctemp.Format(_T("%10.0f |"),m_Activity[i][j]);
			SumActi+=m_Activity[i][j];//统计各层活度
			CTemp+=Ctemp;
		}
		Ctemp.Format(_T("%10.0f |"),SumActi);
		CTemp+=Ctemp;
		output<<(LPCTSTR)CTemp<<endl;
	}

	CTemp.Format(_T("|Total Activity |"));
	for(j=0;j<m_ActivityOfWD->NuDataUsed;j++)
		{
			Ctemp.Format(_T("%10.0f |"),m_MeanActivity[j]);
			CTemp+=Ctemp;
		}
	Ctemp.Format(_T("%10.0f |"),m_TotalActivity);
	CTemp+=Ctemp;
	output<<(LPCTSTR)CTemp<<endl;
	output.close();

	ofstream output2;
	output2.open(m_FolderDetection+_T("#ReconstructionReportList.txt"),ios::app);
	Ctemp.Format(_T("%25s %10.4f "),m_ProjectName,m_MeanDensity);
	output2<<(LPCTSTR)Ctemp;
	CTemp.Empty();
	for(j=0;j<m_ActivityOfWD->NuDataUsed;j++)
		{
			Ctemp.Format(_T("%10.0f "),m_MeanActivity[j]);
			CTemp+=Ctemp;
		}
	output2<<(LPCTSTR)CTemp;
	Ctemp.Format(_T("%10.0f\r\n"),m_TotalActivity);
	output2<<(LPCTSTR)Ctemp;
	output2.close();
	return m_TotalActivity;
}

//统计
void RSCalculation::Statistics()
{
	int i,j,k,ndata;
	//初值
	for(i=0;i<20;i++)
	{
		m_Density[i]=0;
		for(j=0;j<MaxNuDataInTriDNode;j++)
			m_Activity[i][j]=0;	
	}
	m_MeanDensity=0;
	for(j=0;j<MaxNuDataInTriDNode;j++)
		m_MeanActivity[j]=0;	
	m_TotalActivity=0;

	if(!m_DensityOfWD) return;
	if(!m_ActivityOfWD) return;


	TriDNodeDefine *dataz,*datay,*datax;
	//计算密度
	double S,density;
	//透射重建结果统计,将各层各能量各网格内的密度根据密度加权求和平均
	dataz=m_DensityOfWD->FirstNode;
	for(i=0;i<m_DensityOfWD->NuNodesInZ&&dataz;i++,dataz=dataz->NextZ)
	{
		S=0;
		for(j=0,datay=dataz;j<m_DensityOfWD->NuNodesInY&&datay;j++,datay=datay->NextY)
			for(k=0,datax=datay;k<m_DensityOfWD->NuNodesInX&&datax;k++,datax=datax->NextX)
			{
				density=0;
				for(ndata=0;ndata<m_DensityOfWD->NuDataUsed&&ndata<MaxNuDataInTriDNode;ndata++)//不同能量
					density+=datax->data[ndata];
				m_Density[i]+=datax->s*density/m_DensityOfWD->NuDataUsed;	
				S+=datax->s;
			}
		m_H[i]=dataz->z;
		m_Density[i]/=S;
		m_MeanDensity+=m_Density[i];
	}
	m_MeanDensity/=m_DensityOfWD->NuNodesInZ;

	//发射重建结果统计,直接求和
	dataz=m_ActivityOfWD->FirstNode;
	for(i=0;i<m_ActivityOfWD->NuNodesInZ&&dataz;i++,dataz=dataz->NextZ)
		for(j=0,datay=dataz;j<m_ActivityOfWD->NuNodesInY&&datay;j++,datay=datay->NextY)
			for(k=0,datax=datay;k<m_ActivityOfWD->NuNodesInX&&datax;k++,datax=datax->NextX)
				for(ndata=0;ndata<m_ActivityOfWD->NuDataUsed&&ndata<MaxNuDataInTriDNode;ndata++)//不同能量
				{
					m_Activity[i][ndata]+=datax->data[ndata];
					m_MeanActivity[ndata]+=datax->data[ndata];
				}
	m_TotalActivity=0;
	for(ndata=0;ndata<m_ActivityOfWD->NuDataUsed&&ndata<MaxNuDataInTriDNode;ndata++)//不同能量
		m_TotalActivity+=m_MeanActivity[ndata];
}



//读过程记录文件
void RSCalculation::ClassfyingFoldersInf(CString CurrentFolder)
{
	ifstream read(CurrentFolder+_T("_DetInformation.txt"));
	int NuLines=0;
	int i;
	char buf[100];
	CString NoOfWD;
	CString CTemp,CTemp2;
	KeyOfResult=-1;
	while(read)
	{///
		for(i=0;i<100;i++) buf[i]='\0';
		read.getline(buf,100);
		CTemp.Empty();CTemp=buf;
		if(CTemp.GetLength()>0) 
		{
			if(NuLines++==0) NoOfWD=CTemp;
			if(CTemp.Find(_T("=0="))>=0) {CreatingTime=ReadTime(CTemp);KeyOfResult=0;}
			if(CTemp.Find(_T("=1="))>=0) {DetectionTime=ReadTime(CTemp);KeyOfResult=1;}		
			if(CTemp.Find(_T("=2="))>=0) {ReconstructionTime=ReadTime(CTemp);KeyOfResult=2;}	
			if(CTemp.Find(_T("=-1="))>=0) {ReconstructionTime=ReadTime(CTemp);KeyOfResult=-1;}	
		}
	}	
	read.close();
}

//从测量记录里读取时间
CString RSCalculation::ReadTime(CString CTemp)
{
	CString Time;
	if(CTemp.GetLength()==0) return Time;
	if(CTemp.Find(':')<0||CTemp.Find(';')<0) return Time;
	char buf[100],buf2[30];
	int i,j,key=0;
	strcpy(buf,CTemp);
	for(i=0,j=0;i<100&&j<30;i++)
	{
		if(key==1)
		{
			if(buf[i]==';') break;
			buf2[j++]=buf[i];
			buf2[j]='\0';
		}
		if(buf[i]==':') key=1;		
	}	
	Time=buf2;	
	return Time;
}

//根据能量辨别核素
CString RSCalculation::IdentifyNuclide(double Energy)
{
	double EnergySameple[11];
	CString Nuclide[11];
	int No=0,MaxNo;
	//0.32 Cr-51
	EnergySameple[No]=0.32;
	Nuclide[No++]=_T("Cr-51");
	//0.5639 Sb-122
	EnergySameple[No]=0.5639;
	Nuclide[No++]=_T("Sb-122");
	//0.6027 Sb-124
	EnergySameple[No]=0.6027;
	Nuclide[No++]=_T("Sb-124");
	//0.6617 Cs-137
	EnergySameple[No]=0.6617;
	Nuclide[No++]=_T("Cs-137");
	//0.7242 Zr-95
	EnergySameple[No]=0.7242;
	Nuclide[No++]=_T("Zr-95");
	//0.7653 Nb-95
	EnergySameple[No]=0.7653;
	Nuclide[No++]=_T("Nb-95");
	//0.8108 Co-58
	EnergySameple[No]=0.8108;
	Nuclide[No++]=_T("Co-58");
	//0.8348 Mn-54
	EnergySameple[No]=0.8348;
	Nuclide[No++]=_T("Mn-54");
	//1.1732 Co-60
	EnergySameple[No]=1.1732;
	Nuclide[No++]=_T("Co-60");
	//1.3325 Co-60
	EnergySameple[No]=1.3325;
	Nuclide[No++]=_T("Co-60");
	//1.408 Eu-152  
	EnergySameple[No]=1.408;
	Nuclide[No++]=_T("Eu-152");

	MaxNo=No;
	double Ratio=0.01;
    for(No=0;No<MaxNo;No++)
	{
		if(fabs(Energy-EnergySameple[No])<=Ratio*EnergySameple[No]) return Nuclide[No];
	}
	return _T("Unknown");
}

