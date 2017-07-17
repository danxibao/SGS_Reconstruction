// Calculation.cpp: implementation of the CCalculation class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Calculation.h"
#include <time.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif
#include <shlwapi.h>
#pragma comment(lib,"Shlwapi.lib")
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


CCalculation::CCalculation()
{

}

CCalculation::~CCalculation()
{

}

bool CCalculation::DECalculationControl()//�̵߳��ú���
{
	//���̽��Ч���ļ��Ƿ����
	bool keyCreation=true;
	IsDEFilesExist=false;
	CDEB.m_SaveEFolderDE=m_SaveEFolderDE;
	CDEB.m_FolderDE=m_DECalibrationFolder;
	if (!PathIsDirectory(m_SaveEFolderDE)) CreateDirectory(m_SaveEFolderDE,NULL);
	else 
	{
		CFileFind FileFind;
		if(FileFind.FindFile(m_SaveEFolderDE+_T("4SB_DE.dat"))) IsDEFilesExist=true;	
	}
	if(!IsDEFilesExist) //����ļ�������
	{	
		m_Inf+=_T("\r\n    ̽��Ч���ļ������ڣ���������Ч���ļ�...");
		OutputCalculationInf(_T("̽��Ч���ļ������ڣ���ʼ����Ч���ļ�"),false);
		CDEB.DECalSetup();
		//_NdY��CreatDEFiles.h���壬Ϊ̽����ƫ�ĸ�����ÿ��ƫ�Ĺ̶�Ϊ3.5cm
		for(int i=0;i<_NdY;i++)
		{
			if(!CDEB.DE_ProcessingEInCircle(i)) {keyCreation=false;break;}			
			if(!CDEB.DE_ProcessingEInCylinder(i)) {keyCreation=false;break;}
			if(!CDEB.DE_ProcessingEInCBody(i,false)) {keyCreation=false;break;}
			if(!CDEB.DE_ProcessingEInCBody(i,true)) {keyCreation=false;break;}
			m_NProgress=int(100.0*(i+1)/_NdY);
		}
			if(!CDEB.DE_ProcessingEInBody()) keyCreation=false;
		if(!keyCreation)
		{
			m_Inf+=_T("\r\n    ����Ч���ļ�����ʧ��");
			OutputCalculationInf(_T("����Ч���ļ�����ʧ��"),false);
			return false;
		}
		IsDEFilesExist=true;
		m_Inf+=_T("\r\n    ̽��Ч���ļ����ɽ���");
		OutputCalculationInf(_T("̽��Ч���ļ����ɽ���"),false);
	}
	m_Inf+=_T("\r\n    ̽��Ч���ļ�����");
	OutputCalculationInf(_T("̽��Ч���ļ������ɣ��ļ�����"),false);
	Sleep(1000);
	return true;
}

bool CCalculation::CalculationControl()//�̵߳��ú���
{
//	double Energy[1]={0.6617};
//	CDEB.DE_ProcessingEInCBody(0,8,true,Energy,1); 
//return true;

	CString CTemp;
	RSCB.m_SaveEFolderDE=m_SaveEFolderDE;
	RSCB.m_FolderDetection=m_DetectionFolder;

	ofstream output2;
	output2.open(RSCB.m_FolderDetection+_T("#ReconstructionReportList.txt"),ios::trunc);
	//Ctemp.Format(_T("%30s %16s "),_T("��Ŀ����"),_T("�ܶ�(g/cm^3)"),);
//	output2<<(LPCTSTR)Ctemp;
	output2.close();
	//�����ؽ��б����μ���
	while(m_NoOnRecon<m_NuFolders)
	{
		RSCB.m_ProjectName=m_ReconListFloders[m_NoOnRecon];
		m_ReconFolder=m_DetectionFolder+m_ReconListFloders[m_NoOnRecon]+_T("\\");//̽���ļ���
		RSCB.m_Folder=m_ReconFolder;

		CTemp.Format(_T("=> ��ʼ�ؽ���Ŀ%02d��%s"),m_NoOnRecon+1,m_ReconListFloders[m_NoOnRecon]);
		m_Inf+=_T("\r\n ")+CTemp;
		OutputCalculationInf(CTemp,false);

		m_TotalActivity=-1;//��ȳ�ֵ�����С���㣬��ʾ�ؽ�ʧ��
		if(ReconstructionOneProject())
		{
			m_Inf+=_T("\r\n    �ؽ����");
			OutputCalculationInf(_T("   �ؽ����"),false);
		}
		m_FinishedNoProject=m_NoOnRecon;//����ؽ�����Ŀ��
		m_FinishOneProject=true;//����ؽ�
		m_NoOnRecon++;
		m_NProgress=int(100.0*m_NoOnRecon/m_NuFolders);
		Sleep(2000);
	}

	m_Inf+=_T("\r\n ## ��� ##");;
	OutputCalculationInf(_T("## ��� ##"),false);

	Sleep(1000);
	return true;
}

bool CCalculation::ReconstructionOneProject()//��һ����Ŀ�����ؽ�����
{
	RSCB.INIReconCal();

	//��ȡ̽������
	if(!RSCB.VerifyDetection(m_ReconFolder))
		{
			m_Inf+=_T("\r\n    ����̽���ļ��������⣬δ������ؽ�");
			OutputCalculationInf(_T("����̽���ļ��������⣬δ������ؽ�"),false);
			return false;	
		}

	//͸���ؽ����㣬�����Ҫ����͸���ؽ�������У�����Ҫ����ظ�ֵ���ڽ����ܶ�ָ�봴��ʱ����ֵΪ�����ܶȣ�
		if(RSCB.DecSetupData.IsTransmission)//������н���͸������������͸������ؽ�
		{
			if(!RSCB.TransmissionReconstruction())
			{
				m_Inf+=_T("\r\n    ����͸���ؽ�ʧ�ܣ��Գ����ܶȽ��к����ؽ�");
				OutputCalculationInf(_T("����͸���ؽ�ʧ�ܣ��Գ����ܶȽ��к����ؽ�"),false);
			}
		}

		//�����ؽ�����
		RSCB.CreatDEFilesB.m_SaveEFolderDE=m_SaveEFolderDE;
		if(!RSCB.EmissionReconstruction())
		{
				m_Inf+=_T("\r\n    ���󣡷����ؽ�ʧ�ܣ�δ������ؽ�");
				OutputCalculationInf(_T("���󣡷����ؽ�ʧ�ܣ�δ������ؽ�"),false);
				return false;
		}
		//�����¼
		ofstream output;
		output.open(m_ReconFolder+_T("_DetInformation.txt"),ios::app);
		CTime time;
		time = CTime::GetCurrentTime();
		output<<time.Format(_T("\r\n=2= Finished Reconstruction:%Y-%m-%d %H:%M:%S;"))<<endl;
		output.close();

		//���������
		RSCB.OutputReconResults() ;
		m_TotalActivity=RSCB.OutputResults();		

		//ɾ���ڵ�		
		RSCB.DeleteAllPointers();

	return true;
}

void CCalculation::OutputCalculationInf(CString Inf,bool IsCover)//���������Ϣ���ļ�
{
	ofstream output;
	if(IsCover) output.open(m_DetectionFolder+_T("������Ϣ.txt"),ios::trunc);
	else output.open(m_DetectionFolder+_T("������Ϣ.txt"),ios::app);
	CTime time;
    time = CTime::GetCurrentTime();
	output<<time.Format(_T("%Y-%m-%d %H:%M:%S "))<<Inf<<"\n";
	output.close();
}

