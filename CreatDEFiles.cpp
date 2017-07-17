// CreatDEFiles.cpp: implementation of the CCreatDEFiles class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SGS_Reconstruction.h"
#include "CreatDEFiles.h"
#include <math.h>
#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CCreatDEFiles::CCreatDEFiles()
{

}

CCreatDEFiles::~CCreatDEFiles()
{

}

//��ȡ�����ļ������ú���ʱ�����ļ���FullFileName
void CCreatDEFiles::DECalSetup()
{
	CString	FullFileName=m_FolderDE+_T("Detector.txt");
	//���ó�ֵ
	LOfCollimator=15;//׼ֱ�����ȣ���λcm;
	XOfCollimator=38;
	ROfWasteDrum=28;//����Ͱ�뾶����λcm;
	dHOfWasteDrum=10;//����Ͱ��ߣ���λcm;

	//��ȡ�����ļ�
		int i,j;
		CString CTemp,CSpe;
		char buf[300];
		char buft[300];
		int length_buf=300;
		int key=0;
		int keyStart=0;
		ifstream read(FullFileName);//���ļ�
		if(read)
		{
			while(read)
			{///
			for(i=0;i<length_buf;i++) {buf[i]='\0';buft[i]='\0';}
			read.getline(buf,length_buf);//���ж�ȡ�ļ�
			CSpe=buf;
			if(CSpe.Find(_T("END"))>=0||CSpe.Find(_T("End"))>=0||CSpe.Find(_T("end"))>=0) break;//�ļ���ȡ����
			for(i=0,j=0,key=0;i<length_buf;i++) 
			{//
			if(buf[i]==';') break;
			else if(buf[i]==':'||buf[i]=='=') {key=1;}
			else if(key==1&&buf[i]!=' ') buft[j++]=buf[i];
			}//
			CTemp=buft;
			//��ȡʶ������ֵ
			if(CSpe.Find(_T("LOfCollimator"))>=0) LOfCollimator=atof(CTemp);
			else if(CSpe.Find(_T("ROfWasteDrum"))>=0) ROfWasteDrum=atof(CTemp);
			else if(CSpe.Find(_T("XOfCollimator"))>=0) XOfCollimator=atof(CTemp);
			else if(CSpe.Find(_T("dHOfWasteDrum"))>=0) dHOfWasteDrum=atof(CTemp);			
			}///
			read.close();//�رն��ļ�
		}
	XOfCrystal=XOfCollimator+LOfCollimator;//����X���꣬�Է���ͰԲ��Ϊ���������ģ�׼ֱ����X�������򣬵�λcm;
}


//�����ڷ���Ͱ�е�˥������
double CCreatDEFiles::AttenuationLength(CartesianCoordinate PA,double dYOfDrum,double HOfDrum,double ThicknessOfDrum) 
{
	CartesianCoordinate PB;
	PB.x=XOfCrystal;
	PB.y=-1*dYOfDrum;
	PB.z=HOfDrum;
	double R=ROfWasteDrum+ThicknessOfDrum;
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
	if(k1<0) return L*k2;//��ʾ��ԴA�ڷ���Ͱ��
	else return L*(k2-k1);//��ʾ��ԴA�ڷ���Ͱ��
}

//Ǧ������˥��ϵ����cm2/g����������Χ0.1MeV��5MeV��Engergy��λΪMeV
double CCreatDEFiles::CAttenuationPerMassOfPb(double Energy) 
{
	return 10*exp(2.2003*pow(log10(Energy),2.0)-2.1953*pow(log10(Energy),1.0)-4.9528);
}
//���ʲ��ϵ�����˥��ϵ����cm2/g����������Χ0.1MeV��5MeV��Engergy��λΪMeV
double CCreatDEFiles::CAttenuationPerMass(double Energy) 
{
	return 10*exp(-1.0643*log10(Energy)-5.0829);
}

//һ��̽����λ��ʱ������Դ��̽��Ч��
bool CCreatDEFiles::DE_ProcessingEInCircle(int NoPos) 
{
	CString DEName[_NdY];
	double dYOfDrum[_NdY];
	CString CTemp;
	for(int ndY=0;ndY<_NdY;ndY++)
	{	
		CTemp.Format(_T("E_WD%03d.dat"),ndY*35);
		DEName[ndY]=CTemp;
		dYOfDrum[ndY]=3.5*ndY;
	}
	//��ȡһ��̽����λ�õĵ�Դ̽��Ч��
	TriDNodeHeadDefine *DEHead;
	DEHead=TriDataB.ReadData(m_FolderDE+DEName[NoPos]);
	if(!DEHead) return false;

	//������ʱ����
	double Density[11]={0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2,2.25,2.5};
	TriDNodeHeadDefine * HeadTemp=TriDataB.CreatNodes(11,DEHead->NuNodesInY,DEHead->NuNodesInZ);
	if(!HeadTemp) {TriDataB.DeleteNodes(DEHead);return false;}
	HeadTemp->DimensionOfNodes=3;
	HeadTemp->Title=_T("1RL_")+DEName[NoPos];//���ݿռ�����ƣ����޸�,Radius Line;
	int ndata;
	for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<DEHead->NuDataUsed;ndata++)
	{
		HeadTemp->VARIABLESData[ndata]=DEHead->VARIABLESData[ndata];//�ڵ����ݿռ�ĺ����꣬���������˵�������ݣ��Զ�ʶ���double����;
		HeadTemp->VARIABLES[ndata]=DEHead->VARIABLES[ndata];//�ڵ�����˵�������������;
	}
	HeadTemp->TitleXYZ[0]=_T("Density");//x��������;
	HeadTemp->TitleXYZ[1]=_T("Radius");//y��������;
	HeadTemp->TitleXYZ[2]=_T("Height");//z��������;
	HeadTemp->NuDataUsed=DEHead->NuDataUsed;//���ݿռ�ʹ�����ݵĸ���;
	HeadTemp->NuDecimal=8;
	HeadTemp->IsScientific=true;	

	TriDNodeDefine *dataDEx,*dataDEy,*dataDEz;	
	TriDNodeDefine *datax,*datay,*dataz;	
	int i;
	CartesianCoordinate PA;
	double LAtten,DEAtten,u,dA;
	dataDEz=DEHead->FirstNode;
	dataz=HeadTemp->FirstNode;
	while(dataDEz&&dataz)
	{
		dataDEy=dataDEz;
		datay=dataz;
		while(dataDEy&&datay)
		{
			datax=datay;
			for(i=0;i<HeadTemp->NuNodesInX&&datax;i++)
			{
				datax->x=Density[i];//�ܶ�
				datax->y=dataDEy->y;//�뾶
				datax->z=dataDEy->z;//�߶�
				dataDEx=dataDEy;
				while(dataDEx)
				{
					PA.x=dataDEx->y*cos(dataDEx->x*PI/180);
					PA.y=dataDEx->y*sin(dataDEx->x*PI/180);
					PA.z=dataDEy->z;
					LAtten=AttenuationLength(PA,dYOfDrum[NoPos],0,0);//����Ͱ��˥�����ȣ����Ʒ���Ͱ�ں�
					if(dataDEx->NextX&&dataDEx->LastX) dA=0.5*(dataDEx->NextX->x-dataDEx->LastX->x);
					else if(dataDEx->NextX&&!dataDEx->LastX) dA=(dataDEx->NextX->x-dataDEx->x);
					else if(!dataDEx->NextX&&dataDEx->LastX) dA=(dataDEy->x+360-dataDEx->LastX->x);	
					for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<DEHead->NuDataUsed;ndata++)//��ͬ����
					{
						DEAtten=dataDEx->data[ndata];//ԭʼ�ĵ�Դ̽��Ч��
						u=CAttenuationPerMass(HeadTemp->VARIABLESData[ndata])*Density[i];//˥��ϵ��u�����������u/d���ٳ����ܶ�
						DEAtten*=exp(-1*u*LAtten);//˥����ĵ�Դ̽��Ч��
						datax->data[ndata]+=DEAtten*dA;//ͳ���߻�Դ��˥��У����̽��Ч��
					}
					dataDEx=dataDEx->NextX;
				}
				datax=datax->NextX;
			}
			datax=datay;
			for(i=0;i<HeadTemp->NuNodesInX&&datax;i++)
			{
				for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<DEHead->NuDataUsed;ndata++)//��ͬ����
					datax->data[ndata]/=360.0;//ͳ���߻�Դ��˥��У����̽��Ч��
				datax=datax->NextX;
			}
			dataDEy=dataDEy->NextY;
			datay=datay->NextY;
		}	
		dataDEz=dataDEz->NextZ;
		dataz=dataz->NextZ;
	}

	TriDataB.outputdata_3d(HeadTemp,m_SaveEFolderDE+HeadTemp->Title,_T("yxz"));
	TriDataB.DeleteNodes(HeadTemp);
	TriDataB.DeleteNodes(DEHead);
	return true;
}
//һ��̽����λ��ʱ������Դ��̽��Ч��
bool CCreatDEFiles::DE_ProcessingEInCylinder(int NoPos) 
{
	CString DEName[_NdY];
	double dYOfDrum[_NdY];
	CString CTemp;
	for(int ndY=0;ndY<_NdY;ndY++)
	{	
		CTemp.Format(_T("E_WD%03d.dat"),ndY*35);
		DEName[ndY]=CTemp;
		dYOfDrum[ndY]=3.5*ndY;
	}
	//��ȡһ��̽����λ�õĵ�Դ̽��Ч��
	TriDNodeHeadDefine *DEHead;
	CString FileName=_T("1RL_")+DEName[NoPos];
	DEHead=TriDataB.ReadData(m_SaveEFolderDE+FileName);
	if(!DEHead) return false;

	//������ʱ����
	double HightBody[3]={0,10,20};
	TriDNodeHeadDefine * HeadTemp=TriDataB.CreatNodes(DEHead->NuNodesInX,DEHead->NuNodesInY,3);
	if(!HeadTemp) {TriDataB.DeleteNodes(DEHead);return false;}
	HeadTemp->DimensionOfNodes=3;
	HeadTemp->Title=_T("2RS_")+DEName[NoPos];//���ݿռ�����ƣ����޸�, Radius Sylinder;
	int ndata;
	for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<DEHead->NuDataUsed;ndata++)
	{
		HeadTemp->VARIABLESData[ndata]=DEHead->VARIABLESData[ndata];//�ڵ����ݿռ�ĺ����꣬���������˵�������ݣ��Զ�ʶ���double����;
		HeadTemp->VARIABLES[ndata]=DEHead->VARIABLES[ndata];//�ڵ�����˵�������������;
	}
	HeadTemp->TitleXYZ[0]=_T("Density");//x��������;
	HeadTemp->TitleXYZ[1]=_T("Radius");//y��������;
	HeadTemp->TitleXYZ[2]=_T("Height");//z��������;
	HeadTemp->NuDataUsed=DEHead->NuDataUsed;//���ݿռ�ʹ�����ݵĸ���;
	HeadTemp->NuDecimal=8;
	HeadTemp->IsScientific=true;

	TriDNodeDefine *dataDEx,*dataDEy,*dataDEz;	
	TriDNodeDefine *datax,*datay,*dataz;	
	double X1[20],Y1[20],X2[251],Y2[251];
	int N1,N2=251;
	double MaxH,dH,h;
	double SumEh1,SumEh2,SumEh3;
	int i,n1,n2,n3;
	dataDEx=DEHead->FirstNode;
	datax=HeadTemp->FirstNode;
	while(dataDEx&&datax)
	{
		dataDEy=dataDEx;
		datay=datax;
		while(dataDEy&&datay)
		{
			dataz=datay;
			//���긳ֵ
			for(i=0;i<3&&dataz;i++)
			{
				dataz->x=dataDEy->x;
				dataz->y=dataDEy->y;
				dataz->z=HightBody[i];	
				dataz=dataz->NextZ;
			}
			//�����������μ���
			for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<DEHead->NuDataUsed;ndata++)//��ͬ����
			{
				dataDEz=dataDEy;
				for(N1=0;N1<20&&N1<DEHead->NuNodesInZ&&dataDEz;N1++)
				{
					X1[N1]=dataDEz->z;
					MaxH=dataDEz->z;
					Y1[N1]=dataDEz->data[ndata];
					dataDEz=dataDEz->NextZ;
				}
			//2.5���ϲ�߶��ڷֳ�251��
				dH=MaxH/(N2-1);//�ֳɵĲ��
				for(i=0;i<N2;i++)//����ĸ߶�
					X2[i]=dH*i;
				TriDataB.SplineInterpolation(X1,Y1,N1,X2,Y2,N2,true);//������ֵ
				SumEh1=0;SumEh2=0;SumEh3=0;//��ֵ
				n1=0;n2=0;n3=0;
				for(i=0;i<N2-1;i++)
				{
					h=dH*i;
					if(h<MaxH*0.2) {SumEh1+=Y2[i];n1++;}//����ǵ�ǰ��
					else if(h<MaxH*0.6) {SumEh2+=Y2[i];n2++;}//�������һ��
					else {SumEh3+=Y2[i];n3++;}//�����������		
				}
				dataz=datay;
				for(i=0;i<3&&dataz;i++)//ƽ��ֵ
				{
					if(i==0) dataz->data[ndata]=SumEh1/n1;
					else if(i==1) dataz->data[ndata]=SumEh2/n2;
					else dataz->data[ndata]=SumEh3/n3;
					dataz=dataz->NextZ;
				}
			}	
			dataDEy=dataDEy->NextY;
			datay=datay->NextY;
		}	
		dataDEx=dataDEx->NextX;
		datax=datax->NextX;
	}

	TriDataB.outputdata_3d(HeadTemp,m_SaveEFolderDE+HeadTemp->Title,_T("yxz"));//���������
	TriDataB.DeleteNodes(HeadTemp);//ɾ���ڵ�
	TriDataB.DeleteNodes(DEHead);//ɾ���ڵ�
	return true;
}

//һ��̽����λ��ʱ�β��ڻ������̽��Ч��
bool CCreatDEFiles::DE_ProcessingEInCBody(int NoPos,bool IsEqualVolume) 
{
	CString DEName[_NdY];
	double dYOfDrum[_NdY];
	CString CTemp;
	for(int ndY=0;ndY<_NdY;ndY++)
	{	
		CTemp.Format(_T("E_WD%03d.dat"),ndY*35);
		DEName[ndY]=CTemp;
		dYOfDrum[ndY]=3.5*ndY;
	}
	//��ȡһ��̽����λ�õĵ�Դ̽��Ч��
	TriDNodeHeadDefine *DEHead;
	CString FileName=_T("2RS_")+DEName[NoPos];
	DEHead=TriDataB.ReadData(m_SaveEFolderDE+FileName);
	if(!DEHead) return false;

	//������ʱ����
	int i,j;
	int NC=8;
	double RB[8];
	double RBH[8];
	double SC[8];
	for(i=0;i<NC;i++)
	{
	
		if(!IsEqualVolume) {RB[i]=(i+1)*1.0/NC;RBH[i]=(i+0.5)*1.0/NC;}//���ȷֲ�
		else {RB[i]=sqrt((i+1.0)/NC);RBH[i]=sqrt((i+0.5)/NC);}//���ȷֲ�	
	}
	TriDNodeHeadDefine * HeadTemp=TriDataB.CreatNodes(DEHead->NuNodesInX,NC,DEHead->NuNodesInZ);
	if(!HeadTemp) {TriDataB.DeleteNodes(DEHead);return false;}
	HeadTemp->DimensionOfNodes=3;
	if(IsEqualVolume) HeadTemp->Title=_T("3SCB_EB_")+DEName[NoPos];//�����, Segment Cylinder Body, Equal Body;
	else HeadTemp->Title=_T("3SCB_ER_")+DEName[NoPos];//�Ȱ뾶,Segment Cylinder Body, Equal Radius;
	int ndata;
	for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<DEHead->NuDataUsed;ndata++)
	{
		HeadTemp->VARIABLESData[ndata]=DEHead->VARIABLESData[ndata];//�ڵ����ݿռ�ĺ����꣬���������˵�������ݣ��Զ�ʶ���double����;
		HeadTemp->VARIABLES[ndata]=DEHead->VARIABLES[ndata];//�ڵ�����˵�������������;
	}
	HeadTemp->TitleXYZ[0]=_T("Density");//x��������;
	HeadTemp->TitleXYZ[1]=_T("Radius");//y��������;
	HeadTemp->TitleXYZ[2]=_T("Height");//z��������;
	HeadTemp->NuDataUsed=DEHead->NuDataUsed;//���ݿռ�ʹ�����ݵĸ���;
	HeadTemp->NuDecimal=8;
	HeadTemp->IsScientific=true;

	TriDNodeDefine *dataDEx,*dataDEy,*dataDEz;	
	TriDNodeDefine *datax,*datay,*dataz;	
	double X1[50],Y1[50],X2[281],Y2[281];
	int N1,N2=281;
	double MaxR,dR,r;
	dataDEx=DEHead->FirstNode;
	datax=HeadTemp->FirstNode;
	while(dataDEx&&datax)
	{
		dataDEz=dataDEx;
		dataz=datax;
		while(dataDEz&&dataz)
		{
			//�����������μ���
			for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<DEHead->NuDataUsed;ndata++)//��ͬ����
			{
				dataDEy=dataDEz;MaxR=0;
				while(dataDEy&&dataDEy->NextY) dataDEy=dataDEy->NextY;//ע�⣺����ֵ����Ϊ�뾶�����Ǵ������ڽ��б��
				for(N1=0;N1<50&&N1<DEHead->NuNodesInY&&dataDEy;N1++)
				{
					X1[N1]=dataDEy->y;
					if(MaxR<dataDEy->y)MaxR=dataDEy->y;
					Y1[N1]=dataDEy->data[ndata];
					dataDEy=dataDEy->LastY;
				}
			// �뾶��Χ�ڷֳ�281��
				dR=MaxR/(N2-1);//�뾶�Ľڵ㳤��
				for(i=0;i<N2;i++)//���ڵ�뾶
					X2[i]=dR*i;
				TriDataB.SplineInterpolation(X1,Y1,N1,X2,Y2,N2,true);//������ֵ
				for(j=0;j<NC;j++)
					SC[j]=0;
				datay=dataz;
				for(i=1,j=0;i<N2&&j<NC&&datay;i++)
				{
					r=dR*i;
					if(r>RB[j]*MaxR&&datay) {datay=datay->NextY;datay->data[ndata]=0;j++;}
					if(r<=RB[j]*MaxR&&datay)
					{
						datay->data[ndata]+=(Y2[i]+Y2[i-1])*dR*dR*(i-0.5);
					    SC[j]+=2*dR*dR*(i-0.5);					
					}//����ǵ�ǰ��
				}
				//��ƽ��
				datay=dataz;
				for(j=0;j<NC&&datay;j++)
				{
					//���긳ֵ
					if(ndata==0)
					{
					datay->x=dataDEz->x;
					datay->y=RBH[j]*MaxR;
					datay->z=dataDEz->z;
					}
					datay->data[ndata]/=SC[j];
					datay=datay->NextY;					
				}
			}	
			dataz=dataz->NextZ;
			dataDEz=dataDEz->NextZ;
		}	
		dataDEx=dataDEx->NextX;
		datax=datax->NextX;
	}
	TriDataB.outputdata_3d(HeadTemp,m_SaveEFolderDE+HeadTemp->Title,_T("yxz"));//���������
	TriDataB.DeleteNodes(HeadTemp);//ɾ���ڵ�
	TriDataB.DeleteNodes(DEHead);//ɾ���ڵ�
	return true;
}

//һ��̽����λ��ʱ�β��ڻ������̽��Ч��
TriDNodeHeadDefine * CCreatDEFiles::DE_ProcessingEInCBody(int NoPos,int NuDivisionInR,bool IsEqualVolume) 
{
	//�ڰ뾶������������󻮷���
	if(NuDivisionInR>MaxNuDivisionInR) NuDivisionInR=MaxNuDivisionInR;
	//���ڵ�Բ����Ч�ʣ���10��ƫ��
	CString DEName[_NdY];
	double dYOfDrum[_NdY];
	CString CTemp;
	for(int ndY=0;ndY<_NdY;ndY++)
	{	
		CTemp.Format(_T("E_WD%03d.dat"),ndY*35);
		DEName[ndY]=CTemp;
		dYOfDrum[ndY]=3.5*ndY;
	}
	//��ȡһ��̽����λ�õĵ�Դ̽��Ч��
	TriDNodeHeadDefine *DEHead;
	CString FileName=_T("2RS_")+DEName[NoPos];
	DEHead=TriDataB.ReadData(m_SaveEFolderDE+FileName);
	if(!DEHead) return NULL;
	//������ʱ����
	int i,j;
	int NC=NuDivisionInR;

	double RB[MaxNuDivisionInR];
	double RBH[MaxNuDivisionInR];
	double SC[MaxNuDivisionInR];
	double BX2[251];
	//��ֵ��x����
	for(i=0;i<251;i++)
		BX2[i]=pow(1.0*i/250,2);
	//���񻮷ְ뾶����Ľڵ�����
	for(i=0;i<NC;i++)
	{
		if(!IsEqualVolume) {RB[i]=(i+1)*1.0/NC;	RBH[i]=(i+0.5)/NC;}//���Ⱦ��ȷֲ�
		else {RB[i]=sqrt((i+1.0)/NC);RBH[i]=sqrt((i+0.5)/NC);}//������ֲ�	
	}
	//�ڵ������滻�������
	for(j=0;j<NC;j++)
	{
		for(i=0;i<250;i++)
		{
			if(BX2[i]<=RB[j]&&RB[j]<=BX2[i+1])
			{
				if(fabs(RB[j]-BX2[i])<=fabs(BX2[i+1]-RB[j])) BX2[i]=RB[j];
				else BX2[i+1]=RB[j];
				break;
			}
		}
	}
	//����ָ��
	TriDNodeHeadDefine * HeadTemp=TriDataB.CreatNodes(DEHead->NuNodesInX,NC,DEHead->NuNodesInZ);
	if(!HeadTemp) {TriDataB.DeleteNodes(DEHead);return NULL;}
	HeadTemp->DimensionOfNodes=3;
	//DIR��Division in Radius
	if(IsEqualVolume) {CTemp.Format(_T("DIR%02d_SCB_EB_"),NC);HeadTemp->Title=CTemp+DEName[NoPos];}//�����, Segment Cylinder Body, Equal Body;
	else {CTemp.Format(_T("DIR%02d_SCB_ER_"),NC);HeadTemp->Title=CTemp+DEName[NoPos];}//�Ȱ뾶,Segment Cylinder Body, Equal Radius;
	int ndata;
	for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<DEHead->NuDataUsed;ndata++)
	{
		HeadTemp->VARIABLESData[ndata]=DEHead->VARIABLESData[ndata];//�ڵ����ݿռ�ĺ����꣬���������˵�������ݣ��Զ�ʶ���double����;
		HeadTemp->VARIABLES[ndata]=DEHead->VARIABLES[ndata];//�ڵ�����˵�������������;
	}
	HeadTemp->IsScientific=true;
	HeadTemp->TitleXYZ[0]=_T("Density");//x��������;
	HeadTemp->TitleXYZ[1]=_T("Radius");//y��������;
	HeadTemp->TitleXYZ[2]=_T("Height");//z��������;
	HeadTemp->NuDataUsed=DEHead->NuDataUsed;//���ݿռ�ʹ�����ݵĸ���;
	HeadTemp->NuDecimal=8;
	HeadTemp->IsScientific=true;

	TriDNodeDefine *dataDEx,*dataDEy,*dataDEz;	
	TriDNodeDefine *datax,*datay,*dataz;	
	double X1[50],Y1[50],X2[251],Y2[251];
	int N1,N2=251;
	double MaxR,dR,r;
	dataDEx=DEHead->FirstNode;
	datax=HeadTemp->FirstNode;
	while(dataDEx&&datax)
	{
		dataDEz=dataDEx;
		dataz=datax;
		while(dataDEz&&dataz)
		{
			//�����������μ���
			for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<DEHead->NuDataUsed;ndata++)//��ͬ����
			{
				dataDEy=dataDEz;MaxR=0;
				while(dataDEy&&dataDEy->NextY) dataDEy=dataDEy->NextY;//ע�⣺����ֵ����Ϊ�뾶�����Ǵ������ڽ��б��
				for(N1=0;N1<50&&N1<DEHead->NuNodesInY&&dataDEy;N1++)
				{
					X1[N1]=dataDEy->y;
					if(MaxR<dataDEy->y)MaxR=dataDEy->y;
					Y1[N1]=dataDEy->data[ndata];
					dataDEy=dataDEy->LastY;
				}
			// �뾶��Χ�ڷֳ�251��
				for(i=0;i<N2;i++)//���ڵ�뾶
					X2[i]=BX2[i]*MaxR;
				TriDataB.SplineInterpolation(X1,Y1,N1,X2,Y2,N2,true);//������ֵ
				for(j=0;j<NC;j++)
					SC[j]=0;
				datay=dataz;
				for(i=1,j=0;i<N2&&j<NC&&datay;i++)
				{
					r=X2[i];
					if(r>RB[j]*MaxR&&datay) {datay=datay->NextY;datay->data[ndata]=0;j++;}
					if(r<=RB[j]*MaxR&&datay)
					{
						dR=X2[i]-X2[i-1];
						datay->data[ndata]+=(Y2[i]+Y2[i-1])*dR*dR*(i-0.5);
					    SC[j]+=2*dR*dR*(i-0.5);					
					}
				}
				//��ƽ��
				datay=dataz;
				for(j=0;j<NC&&datay;j++)
				{
					//���긳ֵ
					if(ndata==0)
					{
						datay->x=dataDEz->x;
						datay->y=RBH[j]*MaxR;
						datay->z=dataDEz->z;
					}
					datay->data[ndata]/=SC[j];
					datay=datay->NextY;					
				}
			}	
			dataz=dataz->NextZ;
			dataDEz=dataDEz->NextZ;
		}	
		dataDEx=dataDEx->NextX;
		datax=datax->NextX;
	}
	//TriDataB.DeleteNodes(HeadTemp);//ɾ���ڵ�
	TriDataB.DeleteNodes(DEHead);//ɾ���ڵ�
	return HeadTemp;
}

//һ��̽����λ��ʱ�β��ڻ������̽��Ч��,���̶�������
TriDNodeHeadDefine * CCreatDEFiles::DE_ProcessingEInCBody(int NoPos,int NuDivisionInR,bool IsEqualVolume,double Energy[],int NuEnergy) 
{
	CString Ctemp;
	TriDNodeHeadDefine *HeadTemp=DE_ProcessingEInCBody(NoPos,NuDivisionInR,IsEqualVolume);
	if(!HeadTemp) return NULL;
	if(NuEnergy<=0) return NULL;
	TriDNodeDefine *datax,*datay,*dataz;
	int i,j,k,ndata;
	double xSameple[MaxNuDataInTriDNode],ySameple[MaxNuDataInTriDNode],y[MaxNuDataInTriDNode];
	for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<HeadTemp->NuDataUsed;ndata++)
		xSameple[ndata]=HeadTemp->VARIABLESData[ndata];//�̶ȵ�����
	dataz=HeadTemp->FirstNode;
	for(i=0;i<HeadTemp->NuNodesInZ&&dataz;i++,dataz=dataz->NextZ)
	{
		for(j=0,datay=dataz;j<HeadTemp->NuNodesInY&&datay;j++,datay=datay->NextY)
		{
			for(k=0,datax=datay;k<HeadTemp->NuNodesInX&&datax;k++,datax=datax->NextX)
			{
				for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<HeadTemp->NuDataUsed;ndata++)
					ySameple[ndata]=datax->data[ndata];//�̶�Ч��
				TriDataB.SplineInterpolation(xSameple, ySameple, HeadTemp->NuDataUsed, Energy, y,NuEnergy,true);//������ֵ���м���ָ���任
				for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<NuEnergy;ndata++)
					datax->data[ndata]=y[ndata];//�̶�Ч��
			}		
		}
	}
	//����������Ϣ
	HeadTemp->NuDataUsed=NuEnergy;
	CString CTemp;
	for(ndata=0;ndata<MaxNuDataInTriDNode;ndata++)
	{
		HeadTemp->VARIABLESData[ndata]=0;
		HeadTemp->VARIABLES[ndata].Empty();
		if(ndata<HeadTemp->NuDataUsed)
		{
			HeadTemp->VARIABLESData[ndata]=Energy[ndata];
			CTemp.Format(_T("%g"),Energy[ndata]);
			HeadTemp->VARIABLES[ndata]=CTemp;
		}	
	}
		xSameple[ndata]=HeadTemp->VARIABLESData[ndata];//�̶ȵ�����
	TriDataB.outputdata_3d(HeadTemp,m_SaveEFolderDE+HeadTemp->Title,_T("yxz"));//���������
	return HeadTemp;
}

//һ��̽����λ��ʱ�β��ڻ������̽��Ч��
bool CCreatDEFiles::DE_ProcessingEInBody() 
{
	CString DEName[_NdY];
	double dYOfDrum[_NdY];
	int NdY=_NdY;
	CString CTemp;
	for(int ndY=0;ndY<_NdY;ndY++)
	{	
		CTemp.Format(_T("E_WD%03d.dat"),ndY*35);
		DEName[ndY]=CTemp;
		dYOfDrum[ndY]=3.5*ndY;
	}

	//��ȡһ��̽����λ�õĵ�Դ̽��Ч��
	TriDNodeHeadDefine *DEHead=NULL;
	TriDNodeHeadDefine * HeadTemp=NULL;
	CString FileName;
	int ndata;
	TriDNodeDefine *dataDEx,*dataDEy,*dataDEz;	
	TriDNodeDefine *datax,*datay,*dataz;
	for(int i=0;i<NdY;i++)
	{//
		FileName=_T("3SCB_EB_")+DEName[i];
		DEHead=TriDataB.ReadData(m_SaveEFolderDE+FileName);
		if(!DEHead)
		{///	
			TriDataB.DeleteNodes(HeadTemp);//ɾ���ڵ�
			return false;
		}///
		if(i==0)
		{///
			HeadTemp=TriDataB.CreatNodes(DEHead->NuNodesInX,NdY,DEHead->NuNodesInZ);
			if(!HeadTemp) {TriDataB.DeleteNodes(DEHead);return false;}
			HeadTemp->DimensionOfNodes=3;
			HeadTemp->Title=_T("4SB_DE.dat");//Segment Body;
			for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<DEHead->NuDataUsed;ndata++)
			{
			HeadTemp->VARIABLESData[ndata]=DEHead->VARIABLESData[ndata];//�ڵ����ݿռ�ĺ����꣬���������˵�������ݣ��Զ�ʶ���double����;
			HeadTemp->VARIABLES[ndata]=DEHead->VARIABLES[ndata];//�ڵ�����˵�������������;
			}
			HeadTemp->TitleXYZ[0]=_T("Density");//x��������;
			HeadTemp->TitleXYZ[1]=_T("dY");//y��������;
			HeadTemp->TitleXYZ[2]=_T("Height");//z��������;
			HeadTemp->NuDataUsed=DEHead->NuDataUsed;//���ݿռ�ʹ�����ݵĸ���;
			datay=HeadTemp->FirstNode;
			HeadTemp->NuDecimal=8;
			HeadTemp->IsScientific=true;
		}///
		dataDEx=DEHead->FirstNode;
		datax=datay;
		while(dataDEx&&datax)
		{///
			dataDEz=dataDEx;dataz=datax;
			while(dataDEz&&dataz)
			{
				dataz->x=dataDEz->x;
				dataz->z=dataDEz->z;				
				dataz->y=dYOfDrum[i];
				dataDEy=dataDEz;
				while(dataDEy)
				{			
					for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<DEHead->NuDataUsed;ndata++)
						dataz->data[ndata]+=dataDEy->data[ndata];
					dataDEy=dataDEy->NextY;
				}
				for(ndata=0;ndata<MaxNuDataInTriDNode&&ndata<DEHead->NuDataUsed;ndata++)
					dataz->data[ndata]/=DEHead->NuNodesInY;
				dataz=dataz->NextZ;
				dataDEz=dataDEz->NextZ;
			}	
			dataDEx=dataDEx->NextX;
			datax=datax->NextX;
		}///	
		datay=datay->NextY;
		TriDataB.DeleteNodes(DEHead);//ɾ���ڵ�
	}
	TriDataB.outputdata_3d(HeadTemp,m_SaveEFolderDE+HeadTemp->Title,_T("yxz"));//���������
	TriDataB.DeleteNodes(HeadTemp);//ɾ���ڵ�
	return true;
}

//ISGS�ڶ�ȡdy=0��ƫ��λ�õĻ���Ч��
TriDNodeHeadDefine *CCreatDEFiles::ReadRSE_ISGS(double dy) 
{
	CString DEName[_NdY];
	CString DEName2[_NdY];
	double dYOfDrum[_NdY];
	int ndY;
	CString CTemp;
	for(ndY=0;ndY<_NdY;ndY++)
	{	
		CTemp.Format(_T("2RS_E_WD%03d.dat"),ndY*35);
		DEName[ndY]=CTemp;
		CTemp.Format(_T("5RL_E_WD%03d.dat"),ndY*35);
		DEName2[ndY]=CTemp;
		dYOfDrum[ndY]=3.5*ndY;
	}
	for(ndY=0;ndY<_NdY-1;ndY++)
	{
		if(dYOfDrum[ndY]<=dy&&dYOfDrum[ndY+1]>dy) break;
	}
	return TriDataB.ReadData(m_SaveEFolderDE+DEName[ndY]);
}