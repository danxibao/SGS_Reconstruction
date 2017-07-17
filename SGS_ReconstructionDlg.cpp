// SGS_ReconstructionDlg.cpp : implementation file
//

#include "stdafx.h"
#include "SGS_Reconstruction.h"
#include "SGS_ReconstructionDlg.h"
#include <fstream.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	//{{AFX_DATA(CAboutDlg)
	enum { IDD = IDD_ABOUTBOX };
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CAboutDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	//{{AFX_MSG(CAboutDlg)
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
	//{{AFX_DATA_INIT(CAboutDlg)
	//}}AFX_DATA_INIT
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CAboutDlg)
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
	//{{AFX_MSG_MAP(CAboutDlg)
		// No message handlers
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CSGS_ReconstructionDlg dialog

CSGS_ReconstructionDlg::CSGS_ReconstructionDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CSGS_ReconstructionDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CSGS_ReconstructionDlg)
	//}}AFX_DATA_INIT
	// Note that LoadIcon does not require a subsequent DestroyIcon in Win32
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CSGS_ReconstructionDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CSGS_ReconstructionDlg)
	DDX_Control(pDX, IDC_BUTTON1, m_StopB);
	DDX_Control(pDX, IDC_EDIT1, m_InfB);
	DDX_Control(pDX, IDC_LIST1, m_List1B);
	DDX_Control(pDX, IDC_PROGRESS1, m_ProgressB);
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CSGS_ReconstructionDlg, CDialog)
	//{{AFX_MSG_MAP(CSGS_ReconstructionDlg)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_BUTTON1, OnButton1)
	ON_BN_CLICKED(IDC_BUTTON2, OnButton2)
	ON_WM_TIMER()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CSGS_ReconstructionDlg message handlers

BOOL CSGS_ReconstructionDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		CString strAboutMenu;
		strAboutMenu.LoadString(IDS_ABOUTBOX);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon
	
	// TODO: Add extra initialization here
	m_IsFinisedIni=false;
	INIProcess();//��ʼ����
	m_IsFinisedIni=true;
	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CSGS_ReconstructionDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else if((nID   &   0xfff0)   ==   SC_CLOSE)   //����;
	{
		OnButton2();//���ùرմ��ں���
	} 
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CSGS_ReconstructionDlg::OnPaint() 
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, (WPARAM) dc.GetSafeHdc(), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

// The system calls this to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CSGS_ReconstructionDlg::OnQueryDragIcon()
{
	return (HCURSOR) m_hIcon;
}


void CSGS_ReconstructionDlg::OnButton1() //��ͣ������
{
	if(CalculationB.m_IsOnCalculation)
	{
		if(CalculationB.m_StopCalculation)//���������ͣ״̬���ؼ���ʾ���������㡱
		{
			m_StopB.SetWindowText(_T("��������"));
			CalculationB.m_StopCalculation=false;
		}
		else//������ڼ���״̬���ؼ���ʾ����ͣ���㡱
		{
			m_StopB.SetWindowText(_T("��ͣ����"));
			CalculationB.m_StopCalculation=true;
		}
	}
}

void CSGS_ReconstructionDlg::OnButton2() //�˳�
{
	if(CalculationB.m_IsOnCalculation)//����������У�ѯ���Ƿ��˳�
	if(MessageBox(_T("ȷ��Ҫ�˳����㣿"),_T("�˳�����"), MB_ICONQUESTION|MB_OK|MB_OKCANCEL)==IDCANCEL)
	{
		return;//���ȡ�����ͷ���
	}
	time_t timeEnd;
	time(&timeEnd);//��¼ʱ��
	struct tm *localend;
	localend=localtime(&timeEnd);
	char tmpbufend[100];
    strftime(tmpbufend,100,"/ End at %Y-%m-%d %H:%M:%S",localend);
	OutputRecord(tmpbufend,false);

	CalculationB.m_IsOnCalculation=false;
	Sleep(1000);

	CDialog::OnCancel();//�˳�����
	
}
//�����¼��Ϣ,��ֹͬʱ����ͬһ�����
void CSGS_ReconstructionDlg::OutputRecord(CString RecordText,bool IsTrunc) 
{
	ofstream output;
	if(IsTrunc) output.open(m_PresentFolder+_T("LastRunningRecord"),ios::trunc);
	else output.open(m_PresentFolder+_T("LastRunningRecord"),ios::app);
	output<<RecordText;
	output.close();	
}
//��ʼ������ʼ����
void CSGS_ReconstructionDlg::INIProcess() 
{
	TCHAR path[MAX_PATH];
	::GetModuleFileName(NULL,path,MAX_PATH);
	CString p(path);
	CString subp;
	int nPos = p.ReverseFind('\\');
	m_PresentFolder=p.Left(nPos+1);
	SetCalculationFolder();//̽����ļ�·������
	CalculationB.IsDEFilesExist=true;
	CalculationB.m_FinishOneProject=false;
	CalculationB.m_StopCalculation=false;
	CalculationB.m_Inf.Empty();
	CalculationB.m_Inf+=_T(" �ؽ�·����")+CalculationB.m_DetectionFolder;
	//����·���ڵ���ͬ����Ƿ�������
	ifstream read(m_PresentFolder+_T("LastRunningRecord"));  //���ļ�;
	if(read)
	{
		char buf[200]={"\0"};	
		read.getline(buf,200);
		read.close();
		CString CTemp=buf;
		if(CTemp.GetLength()>2&&CTemp.Find(_T("End"))<0)
		{
			//AfxMessageBox(_T("������������У�ͬʱ���лḲ�Ǹ�·���ڵĽ���ļ�"));
			CDialog::OnCancel();//�˳�����		
		}
	}	
	//��¼��ʼ����ʱ��
	time_t time0;
	time(&time0);//��¼ʱ��
	struct tm *local;
	local=localtime(&time0);
	char tmpbuf[100];
    strftime(tmpbuf,100,"Start at %Y-%m-%d %H:%M:%S",local);
	OutputRecord(tmpbuf,true);


	//��ʱ��;
	SetTimer(1,1000,0);
	//�ؼ�����;
	m_ProgressB.SetRange(0,100);
	m_ProgressB.SetPos(0);
	m_LOfInf=0;

	//List
	::SendMessage(m_List1B.m_hWnd, LVM_SETEXTENDEDLISTVIEWSTYLE,LVS_EX_FULLROWSELECT,   LVS_EX_FULLROWSELECT);   
    m_List1B.SetExtendedStyle(
          LVS_EX_FLATSB    // ��ƽ������
        | LVS_EX_FULLROWSELECT    // ��������ѡ��
        | LVS_EX_HEADERDRAGDROP    // ���������ק
        | 0//LVS_EX_ONECLICKACTIVEATE    // ������ʾ
        | LVS_EX_GRIDLINES    // ����������
        );
	m_List1B.InsertColumn(0,"���",LVCFMT_LEFT,40);
	m_List1B.InsertColumn(1,"��Ŀ����",LVCFMT_LEFT,200);
	m_List1B.InsertColumn(2,"״̬",LVCFMT_LEFT,60);
	m_List1B.InsertColumn(3,"�ܻ�� (Bq)",LVCFMT_LEFT,100);
	//�ؽ�·���洢�ռ����
	for(int i=0;i<100;i++)
	CalculationB.m_ReconListFloders[i].Empty();
	CalculationB.m_NuFolders=-1;
	CalculationB.m_NoOnRecon=-1;
	CalculationB.OutputCalculationInf(_T("## ��ʼ ##"),true);
	//���̽��Ч���ļ��Ƿ����
	CalculationB.OutputCalculationInf(_T("���Ч���ļ�"),false);
	CalculationB.IsDEFilesExist=false;
	CalculationB.m_IsOnCalculation=false;
	CheckDeFilesExist();

}

//��̽��ʱ��һЩ�ļ��н�������
void CSGS_ReconstructionDlg::SetCalculationFolder(void)//;
{
	//�ؽ�������ڵ�·��
	CalculationB.m_PresentFolder=m_PresentFolder;
	//Ч�ʿ̶��ļ��洢��·����ͬ�ؽ��������·��ͬ�������ļ����⣩
	CString foldertemp=m_PresentFolder;
	int nPos = foldertemp.Find(_T("\\Reconstruction\\"));
	if(nPos>0) foldertemp=foldertemp.Left(nPos+1);
	CalculationB.m_DECalibrationFolder=foldertemp+("DECalibration\\DE\\");//Ч�ʿ̶ȵ�·��;
	//̽��Ч�ʴ洢��·����ͬ�ؽ������ͬһĿ¼��
	CalculationB.m_SaveEFolderDE=m_PresentFolder+_T("ESaved\\");
	//�洢̽�����ݵ��ļ�·�����ɵ�ǰĿ¼�µĴ洢�ļ�����
	CalculationB.m_DetectionFolder=foldertemp+("Detection\\");
	ifstream read(m_PresentFolder+_T("DetectionFolder"));  //���ļ�;
	char buf[200];
	CString CTemp;
	while(read)
	{
			for(int i=0;i<200;i++) buf[i]='\0';
			read.getline(buf,200);//���ж�ȡ�ļ�;
			CTemp.Empty();CTemp=buf;
			if(CTemp.GetLength()>2) 
			{
				CalculationB.m_DetectionFolder=CTemp;
				break;
			}
	}
}

//��ʼ�������
bool CSGS_ReconstructionDlg::StartCalculation(void)//;
{
	Calculating_ThreadData.CalculatingC=&CalculationB;//������
	Calculating_ThreadData.IsThreadFinishedSuccessful=true;
	AfxBeginThread(Calculating_ThreadFunc,&Calculating_ThreadData);//�����߳�
	return true;
}
//�߳�1�������ؽ�����
UINT Calculating_ThreadFunc(LPVOID lpParam)
{
	Calculating_ThreadInf * pInf=(Calculating_ThreadInf *)lpParam;
	pInf->CalculatingC->m_IsOnPause=false;//��ֵ
	pInf->CalculatingC->m_IsOnCalculation=true;//����״̬
	pInf->CalculatingC->m_NProgress=0;//��������ֵ
	pInf->IsThreadFinishedSuccessful=pInf->CalculatingC->CalculationControl();//������ڳ���
	pInf->CalculatingC->m_IsOnCalculation=false;//�������
	return 1;
}
//���̽��Ч���ļ��Ƿ���ڣ�����������Զ����ɣ������߳�2
void CSGS_ReconstructionDlg::CheckDeFilesExist() 
{
	Calculating_ThreadData.CalculatingC=&CalculationB;//������	
	Calculating_ThreadData.IsThreadFinishedSuccessful=true;
	CalculationB.m_Inf+=_T("\r\n >> ���̽��Ч���ļ�...");
	AfxBeginThread(Calculating_ThreadFunc2,&Calculating_ThreadData);//�����߳�
}
//�߳�2,���ڼ��̽��Ч���ļ��Ƿ���ڣ���������ڣ�����Ч���ļ�
UINT Calculating_ThreadFunc2(LPVOID lpParam)
{
	Calculating_ThreadInf * pInf=(Calculating_ThreadInf *)lpParam;
	pInf->CalculatingC->m_IsOnPause=false;//��ֵ
	pInf->CalculatingC->m_IsOnCalculation=true;//����״̬
	pInf->CalculatingC->m_NProgress=0;//��������ֵ
	pInf->CalculatingC->DECalculationControl();//������ڳ���
	pInf->CalculatingC->m_IsOnCalculation=false;//�������
	return 1;
}

//Timer��ÿ��һ�θ���
void CSGS_ReconstructionDlg::OnTimer(UINT nIDEvent) 
{
	if(m_IsFinisedIni)
	{//
		//�����Ϣ�ַ������ȷ����ı䣬���
		if(m_LOfInf!=CalculationB.m_Inf.GetLength())
		{
			m_InfB.SetWindowText(CalculationB.m_Inf);//ʵʱ��ʾ��Ϣ	
			m_LOfInf=CalculationB.m_Inf.GetLength();
			m_InfB.LineScroll(m_LOfInf,0);
		}
		if(CalculationB.m_IsOnCalculation)//����߳���������
		{
			m_ProgressB.SetPos(CalculationB.m_NProgress);//ʵʱ��ʾ����
		}
		else if(!CalculationB.IsDEFilesExist&&CalculationB.m_NuFolders==-1)//����߳�ֹͣ������Ч���ļ�������
		{
			CalculationB.m_Inf+=_T("\r\n ** Ч���ļ������ڣ����ܽ����ؽ����㣬������ٽ���");
			CalculationB.m_Inf+=_T("\r\n ## ���� ##");
			CalculationB.OutputCalculationInf(_T("Ч���ļ������ڣ����ܽ����ؽ�����"),false);
			CalculationB.OutputCalculationInf(_T("## ���� ##"),false);
			CalculationB.m_NuFolders=0;	
		}
		else//����߳�ֹͣ��Ч���ļ�����
		{
			if(CalculationB.m_NuFolders==-1) ReadReconList();//���m_NuFolders���ǳ�ʼֵ����δ��ȡ�ؽ��ļ��б���ȡ��Ŀ�б��ļ�
			if(CalculationB.m_NuFolders>0&&CalculationB.m_NuFolders<100&&CalculationB.m_NoOnRecon<CalculationB.m_NuFolders&&CalculationB.m_NoOnRecon!=-1)
				StartCalculation();//��ʼ�ؽ�		
		}
		if(CalculationB.m_FinishOneProject)//������һ����Ŀ���ؽ�
		{
			CString  CTemp;
			if(CalculationB.m_TotalActivity<0) m_List1B.SetItemText(CalculationB.m_FinishedNoProject,2, _T("����"));
			else
			{
				m_List1B.SetItemText(CalculationB.m_FinishedNoProject,2, _T("���"));
				CTemp.Format(_T("%.2e"),CalculationB.m_TotalActivity);
				m_List1B.SetItemText(CalculationB.m_FinishedNoProject,3, CTemp);
			}	
			CalculationB.m_FinishOneProject=false;
		}
	}//
	CDialog::OnTimer(nIDEvent);
}

//��ȡ���ؽ�����Ŀ·����
bool CSGS_ReconstructionDlg::ReadReconList() 
{
	CalculationB.m_Inf+=_T("\r\n >> ��ȡ̽����Ŀ��");
	CalculationB.OutputCalculationInf(_T("��ȡ̽����Ŀ��"),false);
	CalculationB.m_NuFolders=0;

	//��list�ļ�����
	CString OldFileName=CalculationB.m_DetectionFolder+_T("ReconList.txt");
	CString NewFileName=CalculationB.m_DetectionFolder+_T("ReconListOld.txt");
	FILE* fp = NULL;
    fp = fopen(NewFileName, "r");
    if(fp != NULL)//��ȡlist�ļ�������ڣ�ֱ���ؽ�
    {
		fclose(fp);
		DeleteFile(NewFileName);
	}
	MoveFile(OldFileName, NewFileName); //��list����Ϊoldlist
	//��ȡoldlist
	CString FullName=CalculationB.m_DetectionFolder+_T("ReconListOld.txt");

	ifstream read(FullName);  //���ļ�;
	if(read==NULL)
	{
		CalculationB.m_Inf+=_T("\r\n    ����δ�ҵ�̽����Ŀ�б��ļ�");
		CalculationB.m_Inf+=_T("\r\n ## ���� ##");
		CalculationB.OutputCalculationInf(_T("δ�ҵ�̽����Ŀ�б��ļ�"),false);
		CalculationB.OutputCalculationInf(_T("## ���� ##"),false);
		return false;	
	}
		int i;		
		char buf[200];
		CString CTemp;
		int NuFolders=0;
		while(read)
		{///
			for(i=0;i<200;i++) {buf[i]='\0';}
			read.getline(buf,200);//���ж�ȡ�ļ�;
			CTemp.Empty();CTemp=buf;
			if(CTemp.GetLength()>10&&(CTemp.Find(_T("SGS"))>=0||CTemp.Find(_T("TGS"))>=0))
				CalculationB.m_ReconListFloders[NuFolders++]=buf;
			if(NuFolders==100)
			{
				CalculationB.m_Inf+=_T("\r\n    ���棡��Ŀ���Ա����100�������ؽ�ǰ100����Ŀ");
				CalculationB.OutputCalculationInf(_T("���棡��Ŀ���Ա����100�������ؽ�ǰ100����Ŀ"),false);
				break;
			}
		}
		read.close();
	//��ȡҪ�ؽ����б�
	if(NuFolders>0)//�����Ŀ��������
	{
		CalculationB.m_NuFolders=NuFolders;
		CTemp.Format(_T("���� %d ����ؽ���Ŀ"),NuFolders);
		CalculationB.m_Inf+=_T("\r\n    ")+CTemp;
		CalculationB.OutputCalculationInf(CTemp,false);
		for(i=0;i<NuFolders;i++)
		{
			CTemp.Format(_T("%3d"),i+1); m_List1B.InsertItem(i,CTemp,LVCFMT_LEFT);
			m_List1B.SetItemText(i,1, CalculationB.m_ReconListFloders[i]);
			m_List1B.SetItemText(i,2, _T("�ȴ�"));
		}
	}
	else //�޴��ؽ���Ŀ
	{
		CalculationB.m_Inf+=_T("\r\n    �����޴��ؽ���Ŀ");
		CalculationB.m_Inf+=_T("\r\n ## ���� ##");
		CalculationB.OutputCalculationInf(_T("�޴��ؽ���Ŀ"),false);
		CalculationB.OutputCalculationInf(_T("## ���� ##"),false);
		return false;	
	}
	//��ʼ����
		CalculationB.m_Inf+=_T("\r\n >> ��ʼ�ؽ�����");
		CalculationB.OutputCalculationInf(_T("��ʼ�ؽ�����"),false);
		CalculationB.m_NoOnRecon=0;//�ӵ�0����Ŀ��ʼ�ؽ�
	return true;
}
