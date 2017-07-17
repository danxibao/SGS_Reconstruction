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

bool CCalculation::DECalculationControl()//线程调用函数
{
	//检查探测效率文件是否存在
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
	if(!IsDEFilesExist) //如果文件不存在
	{	
		m_Inf+=_T("\r\n    探测效率文件不存在，正在生成效率文件...");
		OutputCalculationInf(_T("探测效率文件不存在，开始生成效率文件"),false);
		CDEB.DECalSetup();
		//_NdY在CreatDEFiles.h定义，为探测器偏心个数，每个偏心固定为3.5cm
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
			m_Inf+=_T("\r\n    错误，效率文件生成失败");
			OutputCalculationInf(_T("错误，效率文件生成失败"),false);
			return false;
		}
		IsDEFilesExist=true;
		m_Inf+=_T("\r\n    探测效率文件生成结束");
		OutputCalculationInf(_T("探测效率文件生成结束"),false);
	}
	m_Inf+=_T("\r\n    探测效率文件存在");
	OutputCalculationInf(_T("探测效率文件检测完成，文件存在"),false);
	Sleep(1000);
	return true;
}

bool CCalculation::CalculationControl()//线程调用函数
{
//	double Energy[1]={0.6617};
//	CDEB.DE_ProcessingEInCBody(0,8,true,Energy,1); 
//return true;

	CString CTemp;
	RSCB.m_SaveEFolderDE=m_SaveEFolderDE;
	RSCB.m_FolderDetection=m_DetectionFolder;

	ofstream output2;
	output2.open(RSCB.m_FolderDetection+_T("#ReconstructionReportList.txt"),ios::trunc);
	//Ctemp.Format(_T("%30s %16s "),_T("项目名称"),_T("密度(g/cm^3)"),);
//	output2<<(LPCTSTR)Ctemp;
	output2.close();
	//根据重建列表依次计算
	while(m_NoOnRecon<m_NuFolders)
	{
		RSCB.m_ProjectName=m_ReconListFloders[m_NoOnRecon];
		m_ReconFolder=m_DetectionFolder+m_ReconListFloders[m_NoOnRecon]+_T("\\");//探测文件夹
		RSCB.m_Folder=m_ReconFolder;

		CTemp.Format(_T("=> 开始重建项目%02d：%s"),m_NoOnRecon+1,m_ReconListFloders[m_NoOnRecon]);
		m_Inf+=_T("\r\n ")+CTemp;
		OutputCalculationInf(CTemp,false);

		m_TotalActivity=-1;//活度初值，如果小于零，表示重建失败
		if(ReconstructionOneProject())
		{
			m_Inf+=_T("\r\n    重建完成");
			OutputCalculationInf(_T("   重建完成"),false);
		}
		m_FinishedNoProject=m_NoOnRecon;//完成重建的项目号
		m_FinishOneProject=true;//完成重建
		m_NoOnRecon++;
		m_NProgress=int(100.0*m_NoOnRecon/m_NuFolders);
		Sleep(2000);
	}

	m_Inf+=_T("\r\n ## 完成 ##");;
	OutputCalculationInf(_T("## 完成 ##"),false);

	Sleep(1000);
	return true;
}

bool CCalculation::ReconstructionOneProject()//对一个项目进行重建计算
{
	RSCB.INIReconCal();

	//读取探测设置
	if(!RSCB.VerifyDetection(m_ReconFolder))
		{
			m_Inf+=_T("\r\n    错误！探测文件存在问题，未能完成重建");
			OutputCalculationInf(_T("错误！探测文件存在问题，未能完成重建"),false);
			return false;	
		}

	//透射重建计算，如果需要进行透射重建，则进行，不需要则沉重赋值（在进行密度指针创建时，初值为称重密度）
		if(RSCB.DecSetupData.IsTransmission)//如果进行进行透射测量，则进行透射测量重建
		{
			if(!RSCB.TransmissionReconstruction())
			{
				m_Inf+=_T("\r\n    错误！透射重建失败，以称重密度进行后续重建");
				OutputCalculationInf(_T("错误！透射重建失败，以称重密度进行后续重建"),false);
			}
		}

		//发射重建计算
		RSCB.CreatDEFilesB.m_SaveEFolderDE=m_SaveEFolderDE;
		if(!RSCB.EmissionReconstruction())
		{
				m_Inf+=_T("\r\n    错误！发生重建失败，未能完成重建");
				OutputCalculationInf(_T("错误！发生重建失败，未能完成重建"),false);
				return false;
		}
		//输出记录
		ofstream output;
		output.open(m_ReconFolder+_T("_DetInformation.txt"),ios::app);
		CTime time;
		time = CTime::GetCurrentTime();
		output<<time.Format(_T("\r\n=2= Finished Reconstruction:%Y-%m-%d %H:%M:%S;"))<<endl;
		output.close();

		//输出计算结果
		RSCB.OutputReconResults() ;
		m_TotalActivity=RSCB.OutputResults();		

		//删除节点		
		RSCB.DeleteAllPointers();

	return true;
}

void CCalculation::OutputCalculationInf(CString Inf,bool IsCover)//输出计算信息至文件
{
	ofstream output;
	if(IsCover) output.open(m_DetectionFolder+_T("计算信息.txt"),ios::trunc);
	else output.open(m_DetectionFolder+_T("计算信息.txt"),ios::app);
	CTime time;
    time = CTime::GetCurrentTime();
	output<<time.Format(_T("%Y-%m-%d %H:%M:%S "))<<Inf<<"\n";
	output.close();
}

