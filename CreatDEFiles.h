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
	void DECalSetup();//��ȡ�����ļ������ú���ʱ�����ļ���FullFileName
	double AttenuationLength(CartesianCoordinate PA,double dYOfDrum,double HOfDrum,double ThicknessOfDrum);//�����ڷ���Ͱ�е�˥������ 
	double CAttenuationPerMassOfPb(double Energy);//Ǧ������˥��ϵ����cm2/g����������Χ0.1MeV��5MeV��Engergy��λΪMeV 
	double CAttenuationPerMass(double Energy);//���ʲ��ϵ�����˥��ϵ����cm2/g����������Χ0.1MeV��5MeV��Engergy��λΪMeV
	bool DE_ProcessingEInCircle(int NoPos);//һ��̽����λ��ʱ������Դ��̽��Ч�� 
	bool DE_ProcessingEInCylinder(int NoPos);//һ��̽����λ��ʱ������Դ��̽��Ч��
	bool DE_ProcessingEInCBody(int NoPos,bool IsEqualVolume);//һ��̽����λ��ʱһ���β��̽��Ч��
	TriDNodeHeadDefine *DE_ProcessingEInCBody(int NoPos,int NuDivisionInR,bool IsEqualVolume);//һ��̽����λ��ʱ�β��ڻ������̽��Ч�� 
	TriDNodeHeadDefine *DE_ProcessingEInCBody(int NoPos,int NuDivisionInR,bool IsEqualVolume,double Energy[],int NuEnergy);//һ��̽����λ��ʱ�β��ڻ������̽��Ч��,���̶�������
	bool DE_ProcessingEInBody();//һ��̽����λ��ʱ�β��ڻ������̽��Ч��
	TriDNodeHeadDefine *ReadRSE_ISGS(double dy);//ISGS�ڶ�ȡdy=0��ƫ��λ�õĻ���Ч��

	//̽������һЩ����
	double XOfCollimator;//׼ֱ��X���꣬�Է���ͰԲ��Ϊ���������ģ�׼ֱ����X�������򣬵�λcm;
	double ROfWasteDrum;//����Ͱ�뾶����λcm;
	double LOfCollimator;//׼ֱ�����ȣ���λcm;
	double ROfCollimator;//׼ֱ���뾶����λcm;
	double XOfCrystal;//����X���꣬�Է���ͰԲ��Ϊ���������ģ�׼ֱ����X�������򣬵�λcm;
	double dHOfWasteDrum;//����Ͱ�뾶����λcm;
};

#endif // !defined(AFX_CREATDEFILES_H__6E408369_9D9C_4B72_87DC_97724C8C27DB__INCLUDED_)
