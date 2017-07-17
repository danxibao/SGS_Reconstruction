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
	CString DetectionType;//探测类型
	int nDetectionType;//探测类型
	CString NoOfWD;//废物桶编号
	CString TypeOfWD;//废物桶类型
	CString TypeOfMatrix;//填充物类型
	double GWDD;//废物桶直径
	double GWDH;//废物桶高度
	double GHPGDX;//探测器测量初始位置x坐标，cm
	double GHPGDY;//探测器测量初始位置y坐标，cm
	double GHPGDZ;//探测器测量初始位置z坐标，cm
	int DSNumY;//平动次数
	int DSNumZ;//升降次数
	int DSNumA;//转动次数
	double DSLY;//平动步进距离，cm
	double DSLZ;//升降步进距离，cm
	double DSLA;//转动步进距离，cm
	double DSigleT;//探测时间，min
	double VORotation;//转动次数，min-1
	bool IsTransmission;//是否进行透射测量
	CString TransmissionSource;//填充物类型
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
	//变量
	CString m_Folder;
	CString m_FolderDetection;
	CString m_FolderDE;
	CString m_SaveEFolderDE;
	CString m_ProjectName;
	DecSetupDefine DecSetupData;

	//重建设置
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

	//函数
	void INIReconCal();//初始化
	void DeleteAllPointers();//清除重建过程中的所有指针 
	void OutputReconResults();//输出结果
	bool ReadDetectionSetup(CString FullFileName);//读取探测设置

	bool VerifyDetection(CString Folder);//验证探测的有效性	
	double AttenuationLength(CartesianCoordinate PA,double dYOfDrum,double HOfDrum,double ThicknessOfDrum);//射线在废物桶中的长度
	double AttenuationLength2(CartesianCoordinate PA,double dYOfDrum,double HOfDrum,double D);//射线在废物桶中的长度

	void AttenuationLength(CartesianCoordinate PA,double dYOfDrum,double HOfDrum,double ThicknessOfDrum,double LInSegments[]);//射线在废物桶各段层中的衰减长度
	double CAttenuationPerMassOfPb(double Energy);//铅的质量衰减系数（cm2/kg），能量范围0.1MeV至5MeV，Engergy单位为MeV 
	double CAttenuationPerMass(double Energy);//轻质材料的质量衰减系数（cm2/kg），能量范围0.1MeV至5MeV，Engergy单位为MeV 

	bool CreatDensityOfWD(int NuX,int NuY,int NuZ);//对废物桶创建密度矩阵
	bool CreatActivityHead(int NuX,int NuY,int NuZ);//SGS发射重建
	bool TransmissionReconstruction();//透射重建
	bool EmissionReconstruction();//发射重建
	TriDNodeHeadDefine * CreatCoefArrayHead();//创建系数矩阵
	//ISGS重建
	double rInTime[20][20];
	double NurInTime[20][20][28];
	int nInTime;
	bool EmissionReconstructionISGS();//ISGS发射重建
	bool CreatCoefInISGS();//创建ISGS用的指针空间
	bool IdentifyRadiusInISGS();//判断等效半径
	bool CreatCoefArrayHeadISGS();//创建ISGS的系数矩阵

	double OutputResults();//输出重建结果

	bool SetdataActivityHead(TriDNodeHeadDefine *headX);//创建存储重建活度结果的指针空间
	//读过程记录文件
	CString CreatingTime;
	CString DetectionTime;
	CString ReconstructionTime;
	int KeyOfResult;
	//从测量记录里读取时间
	CString ReadTime(CString CTemp);
	//统计
	void Statistics();
	double m_H[20];
	double m_Density[20];
	double m_Activity[20][MaxNuDataInTriDNode];
	double m_MeanDensity;
	double m_MeanActivity[MaxNuDataInTriDNode];
	void ClassfyingFoldersInf(CString CurrentFolder);
	//根据能量辨别核素
	CString IdentifyNuclide(double Energy);

	//STGS透射重建,求各环密度，再计算平均密度
	bool TransmissionReconstructionSTGS();
};

#endif // !defined(AFX_RSCALCULATION_H__A8C59F77_2E97_4472_AF98_1C077B94E5F8__INCLUDED_)
