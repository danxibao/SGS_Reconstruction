// SGS_ReconstructionDlg.h : header file
//

#if !defined(AFX_SGS_RECONSTRUCTIONDLG_H__94294D1A_0A01_452B_BE25_4AD158A29068__INCLUDED_)
#define AFX_SGS_RECONSTRUCTIONDLG_H__94294D1A_0A01_452B_BE25_4AD158A29068__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

/////////////////////////////////////////////////////////////////////////////
// CSGS_ReconstructionDlg dialog
#include "Calculation.h"

struct Calculating_ThreadInf
{
public:
	CCalculation *CalculatingC;
	bool IsThreadFinishedSuccessful;
};
UINT Calculating_ThreadFunc(LPVOID lpParam);
UINT Calculating_ThreadFunc2(LPVOID lpParam);

class CSGS_ReconstructionDlg : public CDialog
{
// Construction
public:
	CSGS_ReconstructionDlg(CWnd* pParent = NULL);	// standard constructor
	void INIProcess();
	Calculating_ThreadInf Calculating_ThreadData;
	CCalculation CalculationB;
	CString	m_PresentFolder;//��ִ���ļ������ļ���

// Dialog Data
	//{{AFX_DATA(CSGS_ReconstructionDlg)
	enum { IDD = IDD_SGS_RECONSTRUCTION_DIALOG };
	CButton	m_StopB;
	CEdit	m_InfB;
	CListCtrl	m_List1B;
	CProgressCtrl	m_ProgressB;
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CSGS_ReconstructionDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	HICON m_hIcon;
	int m_LOfInf;
	bool m_IsFinisedIni;
	bool StartCalculation(void);

	void CheckDeFilesExist();//���̽���ļ��Ƿ���ڣ�����������Զ����ɣ������߳�2
	//��̽��ʱ��һЩ�ļ��н�������
	void SetCalculationFolder(void);
	//�����¼��Ϣ,��ֹͬʱ����ͬһ�����
	void OutputRecord(CString RecordText,bool IsTrunc);
	//void ReadReconList();//��ȡ�ؽ��б�
	bool ReadReconList();//��ȡҪ�ؽ����ļ���

	// Generated message map functions
	//{{AFX_MSG(CSGS_ReconstructionDlg)
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	afx_msg void OnTimer(UINT nIDEvent);
	afx_msg void OnButton1();
	afx_msg void OnButton2();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_SGS_RECONSTRUCTIONDLG_H__94294D1A_0A01_452B_BE25_4AD158A29068__INCLUDED_)
