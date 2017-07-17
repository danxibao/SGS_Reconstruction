// CreatDEFiles.h: interface for the CCreatDEFiles class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CREATDEFILES_H__6E408369_9D9C_4B72_87DC_97724C8C27DB__INCLUDED_)
#define AFX_CREATDEFILES_H__6E408369_9D9C_4B72_87DC_97724C8C27DB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "TriData.h"
#include <fstream.h>
#include <iostream.h>
#include <iomanip.h>
#define MaxNuDivisionInR 28
#define _NdY 10
class CCreatDEFiles  
{
public:
	CCreatDEFiles();
	virtual ~CCreatDEFiles();
	TriData TriDataB;
	CString m_Folder;
	CString m_FolderDE;
	CString m_SaveEFolderDE;
	void DECalSetup();//读取设置文件，调用函数时给出文件名FullFileName
	double AttenuationLength(CartesianCoordinate PA,double dYOfDrum,double HOfDrum,double ThicknessOfDrum);//射线在废物桶中的衰减长度 
	double CAttenuationPerMassOfPb(double Energy);//铅的质量衰减系数（cm2/g），能量范围0.1MeV至5MeV，Engergy单位为MeV 
	double CAttenuationPerMass(double Energy);//轻质材料的质量衰减系数（cm2/g），能量范围0.1MeV至5MeV，Engergy单位为MeV
	bool DE_ProcessingEInCircle(int NoPos);//一个探测器位置时环形线源的探测效率 
	bool DE_ProcessingEInCylinder(int NoPos);//一个探测器位置时环形线源的探测效率
	bool DE_ProcessingEInCBody(int NoPos,bool IsEqualVolume);//一个探测器位置时一个段层的探测效率
	TriDNodeHeadDefine *DE_ProcessingEInCBody(int NoPos,int NuDivisionInR,bool IsEqualVolume);//一个探测器位置时段层内环形体的探测效率 
	TriDNodeHeadDefine *DE_ProcessingEInCBody(int NoPos,int NuDivisionInR,bool IsEqualVolume,double Energy[],int NuEnergy);//一个探测器位置时段层内环形体的探测效率,给固定的能量
	bool DE_ProcessingEInBody();//一个探测器位置时段层内环形体的探测效率
	TriDNodeHeadDefine *ReadRSE_ISGS(double dy);//ISGS内读取dy=0和偏心位置的环面效率

	//探测器的一些参数
	double XOfCollimator;//准直器X坐标，以废物桶圆心为坐标轴中心，准直器在X正轴区域，单位cm;
	double ROfWasteDrum;//废物桶半径，单位cm;
	double LOfCollimator;//准直器长度，单位cm;
	double ROfCollimator;//准直器半径，单位cm;
	double XOfCrystal;//晶体X坐标，以废物桶圆心为坐标轴中心，准直器在X正轴区域，单位cm;
	double dHOfWasteDrum;//废物桶半径，单位cm;
};

#endif // !defined(AFX_CREATDEFILES_H__6E408369_9D9C_4B72_87DC_97724C8C27DB__INCLUDED_)
