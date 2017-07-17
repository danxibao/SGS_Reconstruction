// RSCalculation.h: interface for the RSCalculation class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RSCALCULATION_H__A8C59F77_2E97_4472_AF98_1C077B94E5F8__INCLUDED_)
#define AFX_RSCALCULATION_H__A8C59F77_2E97_4472_AF98_1C077B94E5F8__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include <fstream.h>
#include <iostream.h>
#include <iomanip.h>
#include "math.h"
#include "TriData.h"
#include "CreatDEFiles.h"

class DecSetupDefine
{
public:
	CString FileName;//
	CString DetectionType;//̽������
	int nDetectionType;//̽������
	CString NoOfWD;//����Ͱ���
	CString TypeOfWD;//����Ͱ����
	CString TypeOfMatrix;//���������
	double GWDD;//����Ͱֱ��
	double GWDH;//����Ͱ�߶�
	double GHPGDX;//̽����������ʼλ��x���꣬cm
	double GHPGDY;//̽����������ʼλ��y���꣬cm
	double GHPGDZ;//̽����������ʼλ��z���꣬cm
	int DSNumY;//ƽ������
	int DSNumZ;//��������
	int DSNumA;//ת������
	double DSLY;//ƽ���������룬cm
	double DSLZ;//�����������룬cm
	double DSLA;//ת���������룬cm
	double DSigleT;//̽��ʱ�䣬min
	double VORotation;//ת��������min-1
	bool IsTransmission;//�Ƿ����͸�����
	CString TransmissionSource;//���������
	double WeightOfDrum;
	int TypeOfIteration;
	int MaxIteration;
	int NuGridsX;
	int NuGridsY;
	int NuGridsZ;
	bool IsEqualVolumeOfGrids;
};

class RSCalculation  
{
public:
	RSCalculation();
	virtual ~RSCalculation();


	TriData TriDataB;
	CCreatDEFiles CreatDEFilesB;
	//����
	CString m_Folder;
	CString m_FolderDetection;
	CString m_FolderDE;
	CString m_SaveEFolderDE;
	CString m_ProjectName;
	DecSetupDefine DecSetupData;

	//�ؽ�����
	int m_MaxIteration_SGS;
	int m_TypeOfIteration;
	TriDNodeHeadDefine *m_DensityOfWD;
	TriDNodeHeadDefine *m_ActivityOfWD;
	TriDNodeHeadDefine *m_TDHead;
	TriDNodeHeadDefine *m_EDHead;
	TriDNodeHeadDefine *CoefArrayHead;
	TriDNodeHeadDefine *CoefArrayHead2;//ISGS
	TriDNodeHeadDefine *CoefArrayHead3;//ISGS
	TriDNodeHeadDefine *RelativeEHead;//ISGS
	TriDNodeHeadDefine *RadiusOfSourceHead;//ISGS
	TriDNodeHeadDefine *RS_E_Head;//ISGS
	TriDNodeHeadDefine *RS_E_Head2;//ISGS
	TriDNodeHeadDefine *EMResidHead;
	TriDNodeHeadDefine *XHead;
	TriDNodeHeadDefine *XHead2;
	double m_TotalActivity;

	//����
	void INIReconCal();//��ʼ��
	void DeleteAllPointers();//����ؽ������е�����ָ�� 
	void OutputReconResults();//������
	bool ReadDetectionSetup(CString FullFileName);//��ȡ̽������

	bool VerifyDetection(CString Folder);//��֤̽�����Ч��	
	double AttenuationLength(CartesianCoordinate PA,double dYOfDrum,double HOfDrum,double ThicknessOfDrum);//�����ڷ���Ͱ�еĳ���
	double AttenuationLength2(CartesianCoordinate PA,double dYOfDrum,double HOfDrum,double D);//�����ڷ���Ͱ�еĳ���

	void AttenuationLength(CartesianCoordinate PA,double dYOfDrum,double HOfDrum,double ThicknessOfDrum,double LInSegments[]);//�����ڷ���Ͱ���β��е�˥������
	double CAttenuationPerMassOfPb(double Energy);//Ǧ������˥��ϵ����cm2/kg����������Χ0.1MeV��5MeV��Engergy��λΪMeV 
	double CAttenuationPerMass(double Energy);//���ʲ��ϵ�����˥��ϵ����cm2/kg����������Χ0.1MeV��5MeV��Engergy��λΪMeV 

	bool CreatDensityOfWD(int NuX,int NuY,int NuZ);//�Է���Ͱ�����ܶȾ���
	bool CreatActivityHead(int NuX,int NuY,int NuZ);//SGS�����ؽ�
	bool TransmissionReconstruction();//͸���ؽ�
	bool EmissionReconstruction();//�����ؽ�
	TriDNodeHeadDefine * CreatCoefArrayHead();//����ϵ������
	//ISGS�ؽ�
	double rInTime[20][20];
	double NurInTime[20][20][28];
	int nInTime;
	bool EmissionReconstructionISGS();//ISGS�����ؽ�
	bool CreatCoefInISGS();//����ISGS�õ�ָ��ռ�
	bool IdentifyRadiusInISGS();//�жϵ�Ч�뾶
	bool CreatCoefArrayHeadISGS();//����ISGS��ϵ������

	double OutputResults();//����ؽ����

	bool SetdataActivityHead(TriDNodeHeadDefine *headX);//�����洢�ؽ���Ƚ����ָ��ռ�
	//�����̼�¼�ļ�
	CString CreatingTime;
	CString DetectionTime;
	CString ReconstructionTime;
	int KeyOfResult;
	//�Ӳ�����¼���ȡʱ��
	CString ReadTime(CString CTemp);
	//ͳ��
	void Statistics();
	double m_H[20];
	double m_Density[20];
	double m_Activity[20][MaxNuDataInTriDNode];
	double m_MeanDensity;
	double m_MeanActivity[MaxNuDataInTriDNode];
	void ClassfyingFoldersInf(CString CurrentFolder);
	//��������������
	CString IdentifyNuclide(double Energy);

	//STGS͸���ؽ�,������ܶȣ��ټ���ƽ���ܶ�
	bool TransmissionReconstructionSTGS();
};

#endif // !defined(AFX_RSCALCULATION_H__A8C59F77_2E97_4472_AF98_1C077B94E5F8__INCLUDED_)
