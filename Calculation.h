// Calculation.h: interface for the CCalculation class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CALCULATION_H__D8DBBC10_C318_4E0D_B240_B1DE7B1EA24D__INCLUDED_)
#define AFX_CALCULATION_H__D8DBBC10_C318_4E0D_B240_B1DE7B1EA24D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "TriData.h"
#include "CreatDEFiles.h"
#include "RSCalculation.h"

class CCalculation  
{
public:
	CCalculation();
	virtual ~CCalculation();
	CString	m_PresentFolder;//当前路径;
	CString	m_SaveEFolderDE;//效率文件存储文件夹;
	bool m_StopCalculation;

	bool m_IsOnCalculation;//是否计算;
	bool m_IsOnPause;//是否暂停;
	int m_NProgress;//进度条数据;
	CString m_Inf;//信息文字;
	CString m_InfSaved;//信息文字;
	RSCalculation RSCB;
	CCreatDEFiles CDEB;
	TriData TriDataB;
	void OutputCalculationInf(CString Inf,bool IsCover);//输出计算信息至文件
	//检查探测效率文件
	bool DECalculationControl();//计算探测效率的入口函数
	CString	m_DECalibrationFolder;//点源效率刻度文件夹
	bool IsDEFilesExist;

	//重建计算
	bool CalculationControl();//入口函数;
	CString	m_DetectionFolder;//探测文件夹
	CString	m_ReconFolder;//重建计算的文件夹，在m_DetectionFolder内
	CString m_ReconListFloders[101];
	int m_NuFolders;
	int m_NoOnRecon;
	bool ReconstructionOneProject();//对一个项目进行重建计算
	double m_TotalActivity;
	bool m_FinishOneProject;
	int m_FinishedNoProject;
};

#endif // !defined(AFX_CALCULATION_H__D8DBBC10_C318_4E0D_B240_B1DE7B1EA24D__INCLUDED_)
