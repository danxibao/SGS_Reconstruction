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
	INIProcess();//起始程序
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
	else if((nID   &   0xfff0)   ==   SC_CLOSE)   //红叉叉;
	{
		OnButton2();//调用关闭窗口函数
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


void CSGS_ReconstructionDlg::OnButton1() //暂停、计算
{
	if(CalculationB.m_IsOnCalculation)
	{
		if(CalculationB.m_StopCalculation)//如果处于暂停状态，控件显示“继续计算”
		{
			m_StopB.SetWindowText(_T("继续计算"));
			CalculationB.m_StopCalculation=false;
		}
		else//如果处于计算状态，控件显示“暂停计算”
		{
			m_StopB.SetWindowText(_T("暂停计算"));
			CalculationB.m_StopCalculation=true;
		}
	}
}

void CSGS_ReconstructionDlg::OnButton2() //退出
{
	if(CalculationB.m_IsOnCalculation)//如果在运行中，询问是否退出
	if(MessageBox(_T("确定要退出计算？"),_T("退出警告"), MB_ICONQUESTION|MB_OK|MB_OKCANCEL)==IDCANCEL)
	{
		return;//如果取消，就返回
	}
	time_t timeEnd;
	time(&timeEnd);//记录时间
	struct tm *localend;
	localend=localtime(&timeEnd);
	char tmpbufend[100];
    strftime(tmpbufend,100,"/ End at %Y-%m-%d %H:%M:%S",localend);
	OutputRecord(tmpbufend,false);

	CalculationB.m_IsOnCalculation=false;
	Sleep(1000);

	CDialog::OnCancel();//退出窗口
	
}
//输出记录信息,防止同时开启同一个软件
void CSGS_ReconstructionDlg::OutputRecord(CString RecordText,bool IsTrunc) 
{
	ofstream output;
	if(IsTrunc) output.open(m_PresentFolder+_T("LastRunningRecord"),ios::trunc);
	else output.open(m_PresentFolder+_T("LastRunningRecord"),ios::app);
	output<<RecordText;
	output.close();	
}
//初始化并开始计算
void CSGS_ReconstructionDlg::INIProcess() 
{
	TCHAR path[MAX_PATH];
	::GetModuleFileName(NULL,path,MAX_PATH);
	CString p(path);
	CString subp;
	int nPos = p.ReverseFind('\\');
	m_PresentFolder=p.Left(nPos+1);
	SetCalculationFolder();//探测的文件路径设置
	CalculationB.IsDEFilesExist=true;
	CalculationB.m_FinishOneProject=false;
	CalculationB.m_StopCalculation=false;
	CalculationB.m_Inf.Empty();
	CalculationB.m_Inf+=_T(" 重建路径：")+CalculationB.m_DetectionFolder;
	//检查该路径内的相同软件是否在运行
	ifstream read(m_PresentFolder+_T("LastRunningRecord"));  //打开文件;
	if(read)
	{
		char buf[200]={"\0"};	
		read.getline(buf,200);
		read.close();
		CString CTemp=buf;
		if(CTemp.GetLength()>2&&CTemp.Find(_T("End"))<0)
		{
			//AfxMessageBox(_T("该软件正在运行，同时运行会覆盖该路径内的结果文件"));
			CDialog::OnCancel();//退出窗口		
		}
	}	
	//记录开始计算时间
	time_t time0;
	time(&time0);//记录时间
	struct tm *local;
	local=localtime(&time0);
	char tmpbuf[100];
    strftime(tmpbuf,100,"Start at %Y-%m-%d %H:%M:%S",local);
	OutputRecord(tmpbuf,true);


	//定时器;
	SetTimer(1,1000,0);
	//控件设置;
	m_ProgressB.SetRange(0,100);
	m_ProgressB.SetPos(0);
	m_LOfInf=0;

	//List
	::SendMessage(m_List1B.m_hWnd, LVM_SETEXTENDEDLISTVIEWSTYLE,LVS_EX_FULLROWSELECT,   LVS_EX_FULLROWSELECT);   
    m_List1B.SetExtendedStyle(
          LVS_EX_FLATSB    // 扁平风格滚动
        | LVS_EX_FULLROWSELECT    // 允许正航选中
        | LVS_EX_HEADERDRAGDROP    // 允许标题拖拽
        | 0//LVS_EX_ONECLICKACTIVEATE    // 高亮显示
        | LVS_EX_GRIDLINES    // 画出网格线
        );
	m_List1B.InsertColumn(0,"序号",LVCFMT_LEFT,40);
	m_List1B.InsertColumn(1,"项目名称",LVCFMT_LEFT,200);
	m_List1B.InsertColumn(2,"状态",LVCFMT_LEFT,60);
	m_List1B.InsertColumn(3,"总活度 (Bq)",LVCFMT_LEFT,100);
	//重建路径存储空间清空
	for(int i=0;i<100;i++)
	CalculationB.m_ReconListFloders[i].Empty();
	CalculationB.m_NuFolders=-1;
	CalculationB.m_NoOnRecon=-1;
	CalculationB.OutputCalculationInf(_T("## 开始 ##"),true);
	//检查探测效率文件是否存在
	CalculationB.OutputCalculationInf(_T("检查效率文件"),false);
	CalculationB.IsDEFilesExist=false;
	CalculationB.m_IsOnCalculation=false;
	CheckDeFilesExist();

}

//对探测时的一些文件夹进行设置
void CSGS_ReconstructionDlg::SetCalculationFolder(void)//;
{
	//重建软件所在的路径
	CalculationB.m_PresentFolder=m_PresentFolder;
	//效率刻度文件存储的路径，同重建软件所在路径同级（即文件夹外）
	CString foldertemp=m_PresentFolder;
	int nPos = foldertemp.Find(_T("\\Reconstruction\\"));
	if(nPos>0) foldertemp=foldertemp.Left(nPos+1);
	CalculationB.m_DECalibrationFolder=foldertemp+("DECalibration\\DE\\");//效率刻度的路径;
	//探测效率存储的路径，同重建软件在同一目录下
	CalculationB.m_SaveEFolderDE=m_PresentFolder+_T("ESaved\\");
	//存储探测数据的文件路径，由当前目录下的存储文件读入
	CalculationB.m_DetectionFolder=foldertemp+("Detection\\");
	ifstream read(m_PresentFolder+_T("DetectionFolder"));  //打开文件;
	char buf[200];
	CString CTemp;
	while(read)
	{
			for(int i=0;i<200;i++) buf[i]='\0';
			read.getline(buf,200);//挨行读取文件;
			CTemp.Empty();CTemp=buf;
			if(CTemp.GetLength()>2) 
			{
				CalculationB.m_DetectionFolder=CTemp;
				break;
			}
	}
}

//开始计算程序
bool CSGS_ReconstructionDlg::StartCalculation(void)//;
{
	Calculating_ThreadData.CalculatingC=&CalculationB;//计算类
	Calculating_ThreadData.IsThreadFinishedSuccessful=true;
	AfxBeginThread(Calculating_ThreadFunc,&Calculating_ThreadData);//调用线程
	return true;
}
//线程1，进行重建计算
UINT Calculating_ThreadFunc(LPVOID lpParam)
{
	Calculating_ThreadInf * pInf=(Calculating_ThreadInf *)lpParam;
	pInf->CalculatingC->m_IsOnPause=false;//初值
	pInf->CalculatingC->m_IsOnCalculation=true;//计算状态
	pInf->CalculatingC->m_NProgress=0;//进度条初值
	pInf->IsThreadFinishedSuccessful=pInf->CalculatingC->CalculationControl();//调用入口程序
	pInf->CalculatingC->m_IsOnCalculation=false;//计算结束
	return 1;
}
//检查探测效率文件是否存在，如果不存在自动生成，调用线程2
void CSGS_ReconstructionDlg::CheckDeFilesExist() 
{
	Calculating_ThreadData.CalculatingC=&CalculationB;//计算类	
	Calculating_ThreadData.IsThreadFinishedSuccessful=true;
	CalculationB.m_Inf+=_T("\r\n >> 检查探测效率文件...");
	AfxBeginThread(Calculating_ThreadFunc2,&Calculating_ThreadData);//调用线程
}
//线程2,用于检查探测效率文件是否存在，如果不存在，生成效率文件
UINT Calculating_ThreadFunc2(LPVOID lpParam)
{
	Calculating_ThreadInf * pInf=(Calculating_ThreadInf *)lpParam;
	pInf->CalculatingC->m_IsOnPause=false;//初值
	pInf->CalculatingC->m_IsOnCalculation=true;//计算状态
	pInf->CalculatingC->m_NProgress=0;//进度条初值
	pInf->CalculatingC->DECalculationControl();//调用入口程序
	pInf->CalculatingC->m_IsOnCalculation=false;//计算结束
	return 1;
}

//Timer，每秒一次更新
void CSGS_ReconstructionDlg::OnTimer(UINT nIDEvent) 
{
	if(m_IsFinisedIni)
	{//
		//如果信息字符串长度发生改变，输出
		if(m_LOfInf!=CalculationB.m_Inf.GetLength())
		{
			m_InfB.SetWindowText(CalculationB.m_Inf);//实时显示信息	
			m_LOfInf=CalculationB.m_Inf.GetLength();
			m_InfB.LineScroll(m_LOfInf,0);
		}
		if(CalculationB.m_IsOnCalculation)//如果线程在运行中
		{
			m_ProgressB.SetPos(CalculationB.m_NProgress);//实时显示进度
		}
		else if(!CalculationB.IsDEFilesExist&&CalculationB.m_NuFolders==-1)//如果线程停止，但是效率文件不存在
		{
			CalculationB.m_Inf+=_T("\r\n ** 效率文件不存在，不能进行重建计算，请检查后再进行");
			CalculationB.m_Inf+=_T("\r\n ## 结束 ##");
			CalculationB.OutputCalculationInf(_T("效率文件不存在，不能进行重建计算"),false);
			CalculationB.OutputCalculationInf(_T("## 结束 ##"),false);
			CalculationB.m_NuFolders=0;	
		}
		else//如果线程停止，效率文件存在
		{
			if(CalculationB.m_NuFolders==-1) ReadReconList();//如果m_NuFolders还是初始值，即未读取重建文件列表，读取项目列表文件
			if(CalculationB.m_NuFolders>0&&CalculationB.m_NuFolders<100&&CalculationB.m_NoOnRecon<CalculationB.m_NuFolders&&CalculationB.m_NoOnRecon!=-1)
				StartCalculation();//开始重建		
		}
		if(CalculationB.m_FinishOneProject)//如果完成一个项目的重建
		{
			CString  CTemp;
			if(CalculationB.m_TotalActivity<0) m_List1B.SetItemText(CalculationB.m_FinishedNoProject,2, _T("错误"));
			else
			{
				m_List1B.SetItemText(CalculationB.m_FinishedNoProject,2, _T("完成"));
				CTemp.Format(_T("%.2e"),CalculationB.m_TotalActivity);
				m_List1B.SetItemText(CalculationB.m_FinishedNoProject,3, CTemp);
			}	
			CalculationB.m_FinishOneProject=false;
		}
	}//
	CDialog::OnTimer(nIDEvent);
}

//读取待重建的项目路径名
bool CSGS_ReconstructionDlg::ReadReconList() 
{
	CalculationB.m_Inf+=_T("\r\n >> 读取探测项目组");
	CalculationB.OutputCalculationInf(_T("读取探测项目组"),false);
	CalculationB.m_NuFolders=0;

	//将list文件改名
	CString OldFileName=CalculationB.m_DetectionFolder+_T("ReconList.txt");
	CString NewFileName=CalculationB.m_DetectionFolder+_T("ReconListOld.txt");
	FILE* fp = NULL;
    fp = fopen(NewFileName, "r");
    if(fp != NULL)//读取list文件如果存在，直接重建
    {
		fclose(fp);
		DeleteFile(NewFileName);
	}
	MoveFile(OldFileName, NewFileName); //将list改名为oldlist
	//读取oldlist
	CString FullName=CalculationB.m_DetectionFolder+_T("ReconListOld.txt");

	ifstream read(FullName);  //打开文件;
	if(read==NULL)
	{
		CalculationB.m_Inf+=_T("\r\n    错误！未找到探测项目列表文件");
		CalculationB.m_Inf+=_T("\r\n ## 结束 ##");
		CalculationB.OutputCalculationInf(_T("未找到探测项目列表文件"),false);
		CalculationB.OutputCalculationInf(_T("## 结束 ##"),false);
		return false;	
	}
		int i;		
		char buf[200];
		CString CTemp;
		int NuFolders=0;
		while(read)
		{///
			for(i=0;i<200;i++) {buf[i]='\0';}
			read.getline(buf,200);//挨行读取文件;
			CTemp.Empty();CTemp=buf;
			if(CTemp.GetLength()>10&&(CTemp.Find(_T("SGS"))>=0||CTemp.Find(_T("TGS"))>=0))
				CalculationB.m_ReconListFloders[NuFolders++]=buf;
			if(NuFolders==100)
			{
				CalculationB.m_Inf+=_T("\r\n    警告！项目组成员超过100个，仅重建前100个项目");
				CalculationB.OutputCalculationInf(_T("警告！项目组成员超过100个，仅重建前100个项目"),false);
				break;
			}
		}
		read.close();
	//读取要重建的列表
	if(NuFolders>0)//如果项目数大于零
	{
		CalculationB.m_NuFolders=NuFolders;
		CTemp.Format(_T("读入 %d 组待重建项目"),NuFolders);
		CalculationB.m_Inf+=_T("\r\n    ")+CTemp;
		CalculationB.OutputCalculationInf(CTemp,false);
		for(i=0;i<NuFolders;i++)
		{
			CTemp.Format(_T("%3d"),i+1); m_List1B.InsertItem(i,CTemp,LVCFMT_LEFT);
			m_List1B.SetItemText(i,1, CalculationB.m_ReconListFloders[i]);
			m_List1B.SetItemText(i,2, _T("等待"));
		}
	}
	else //无待重建项目
	{
		CalculationB.m_Inf+=_T("\r\n    错误！无待重建项目");
		CalculationB.m_Inf+=_T("\r\n ## 结束 ##");
		CalculationB.OutputCalculationInf(_T("无待重建项目"),false);
		CalculationB.OutputCalculationInf(_T("## 结束 ##"),false);
		return false;	
	}
	//开始计算
		CalculationB.m_Inf+=_T("\r\n >> 开始重建计算");
		CalculationB.OutputCalculationInf(_T("开始重建计算"),false);
		CalculationB.m_NoOnRecon=0;//从第0个项目开始重建
	return true;
}
