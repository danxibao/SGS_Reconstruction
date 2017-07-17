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
	CString	m_PresentFolder;//��ǰ·��;
	CString	m_SaveEFolderDE;//Ч���ļ��洢�ļ���;
	bool m_StopCalculation;

	bool m_IsOnCalculation;//�Ƿ����;
	bool m_IsOnPause;//�Ƿ���ͣ;
	int m_NProgress;//����������;
	CString m_Inf;//��Ϣ����;
	CString m_InfSaved;//��Ϣ����;
	RSCalculation RSCB;
	CCreatDEFiles CDEB;
	TriData TriDataB;
	void OutputCalculationInf(CString Inf,bool IsCover);//���������Ϣ���ļ�
	//���̽��Ч���ļ�
	bool DECalculationControl();//����̽��Ч�ʵ���ں���
	CString	m_DECalibrationFolder;//��ԴЧ�ʿ̶��ļ���
	bool IsDEFilesExist;

	//�ؽ�����
	bool CalculationControl();//��ں���;
	CString	m_DetectionFolder;//̽���ļ���
	CString	m_ReconFolder;//�ؽ�������ļ��У���m_DetectionFolder��
	CString m_ReconListFloders[101];
	int m_NuFolders;
	int m_NoOnRecon;
	bool ReconstructionOneProject();//��һ����Ŀ�����ؽ�����
	double m_TotalActivity;
	bool m_FinishOneProject;
	int m_FinishedNoProject;
};

#endif // !defined(AFX_CALCULATION_H__D8DBBC10_C318_4E0D_B240_B1DE7B1EA24D__INCLUDED_)
