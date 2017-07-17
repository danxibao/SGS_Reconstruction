
#pragma once

#define MaxNuDataInTriDNode 10//�ڵ�洢���ݵĸ���;
#define LineChar_MaxLength 3000//�ַ�����ʱ��������;
#define PI 3.14159265358979323846

class CartesianCoordinate
{
public:
	double x;
	double y;
	double z;
};

class TriDNodeDefine
{
public:
	int NoX;
	int NoY;
	int NoZ;
	double x;
	double y;
	double z;
	double data[MaxNuDataInTriDNode];
	double s;
	double xx;
	double yy;
	double zz;
	double SumS;
	CString Inf;
	TriDNodeDefine *NextX;
	TriDNodeDefine *LastX;
	TriDNodeDefine *NextY;
	TriDNodeDefine *LastY;
	TriDNodeDefine *NextZ;
	TriDNodeDefine *LastZ;
	TriDNodeDefine *Link;
	int NuLink;
};

class TriDNodeHeadDefine
{
protected:

public:
	int DimensionOfNodes;//ά����һάֵΪ1��2άΪ2��3άΪ3;
	int NuNodesInX;//X����ڵ����;
	int NuNodesInY;//Y����ڵ����;
	int NuNodesInZ;//Z����ڵ����;
	int NuNodesTotal;//�ڵ�����;
	CString Title;//���ݿռ������;
	CString Inf;//�����ݿռ��ע��;
	CString VARIABLES[MaxNuDataInTriDNode];//double�����ڸ���������������;
	CString TitleXYZ[3];
	double VARIABLESData[MaxNuDataInTriDNode];//ÿ���ڵ������˵��;
	double MaxData[MaxNuDataInTriDNode];
	double MinData[MaxNuDataInTriDNode];
	TriDNodeDefine *FirstNode;
	TriDNodeDefine *EndNode;
	TriDNodeDefine *UsingNode;
	TriDNodeDefine *UsedNode;
	int NuDataUsed;//data����������ռ��ʹ�ø�����ӦС��MaxNuDataInTriDNode;
	int NuNodesUsed;//ʹ�ýڵ������;
	int NuDecimal;//С����λ��;
	bool IsScientific;//�Ƿ��ѧ������;
	bool IsOutputInf;//�Ƿ��ѧ������;
};

class TriData  
{
public:
	int KeyDirection[3];
	int NuIJK;//��ȡdata�ļ�ʱ,IJK���ڵĸ�����������Ϊ3;
	int IsExistI;//��ȡdata�ļ�ʱ���Ƿ����I=,������Ϊ-1������ΪIJK���д����ţ��㿪ʼ;
	int IsExistJ;//��ȡdata�ļ�ʱ���Ƿ����I=,������Ϊ-1������ΪIJK���д����ţ��㿪ʼ;
	int IsExistK;//��ȡdata�ļ�ʱ���Ƿ����I=,������Ϊ-1������ΪIJK���д����ţ��㿪ʼ;
	TriDNodeHeadDefine *TriDataUsing;
	TriDNodeHeadDefine *TriDataUsingTemp;
	double m_RelaxFactor;
	//��ά���ݿռ�ʹ�ó�ʼ��;
	void TriDataInitialization(void);
	//�����ڵ�,x,y,z���������ϵ���������ʼ���ݶ�����;
	TriDNodeHeadDefine * CreatNodes(int NuX, int NuY, int NuZ);
	//����ͷ�ڵ�;
	TriDNodeHeadDefine * CreatNodeHead();
	//����head�ڵ�;
	void CopyHead(TriDNodeHeadDefine *nodehead1,TriDNodeHeadDefine *nodehead2);
	//��head1�е���Ϣ���Ƹ�head2�����ǽڵ�������;
	void CopyHeadInf(TriDNodeHeadDefine *nodehead1,TriDNodeHeadDefine *nodehead2);
	//������ά���ݽڵ㣬ͷ��Ϣ������һ��;
	TriDNodeHeadDefine *TriData::CopyNodes(TriDNodeHeadDefine *nodehead);
	//���ɲ���ʼ���ڵ�,�ڵ�š��������㣬����ָ����NULL;
	TriDNodeDefine * CreatNode();
	//�����ڵ���ͷ���,x,y,z���������ϵ���������ʼ���ݶ�����;
	TriDNodeDefine * CreatNodesWithoutHead(int NuX, int NuY, int NuZ);
	//�ƶ��ڵ㣬�ƶ�һ��,��ǰ�ڵ㣬�ƶ����򣬷����ƶ���Ľڵ�ָ�룬����1��2��3�ֱ�Ϊx��y��z���򣬸��ű�ʾ������;
	TriDNodeDefine * MoveNode(TriDNodeDefine *node,int direction);
	//�����ڵ�֮������ӣ�node1ָ��node2,direction=1 x����2 y����3 z����;
	void ConnectTwoNodes(TriDNodeDefine *node1,TriDNodeDefine *node2,int direction);
	//���нڵ�֮������ӣ�node1���ڵ���ָ��node2���ڵ���,direction=1 x����2 y����3 z����;
	void ConnectTwoLineNodes(TriDNodeDefine *node1,TriDNodeDefine *node2,int directionMove,int directionLink);
	//����ڵ�֮������ӣ�node1������ָ��node2������,direction=1 x����2 y����3 z����;
	void ConnectTwoFaceNodes(TriDNodeDefine *node1,TriDNodeDefine *node2,int directionLink);
	//���ͷָ��ɾ���ڵ�;
	void DeleteNodes(TriDNodeHeadDefine *nodehead);
	//���ڵ������ȫ������;
	void SetNodesData(TriDNodeHeadDefine *nodehead,double data);
	//����ڵ����ݵ��ļ�����ά��ʽ��������IJK����ʹΪ1��,��ʽ�ɱ�tecplot���룬����Ϊ��ͷָ�룬ȫ��������·���������ݴ洢�����ȷ����硰XYZ������ZYX���ȣ�;
	void outputdata_3d(TriDNodeHeadDefine *nodehead,CString FullName,CString Direction);
	void outputdata_3dJZB(TriDNodeHeadDefine *nodehead,CString FullName);//����������£����㿨�����������(���ܶ���);
	//����ڵ����ݵ��ļ�����ά��ʽ��,��ʽ�ɱ�tecplot���룬����Ϊ��ͷָ�룬ȫ��������·������ƽ�棨�硰XY������ZY���ȣ����ڼ�ƽ��;
	void outputdata_2d(TriDNodeHeadDefine *nodehead,CString FullName,CString Plane,int NoPlane);
	//���ڵ������ȫ������ĳ��;
	void MultiplyNodesData(TriDNodeHeadDefine *nodehead,double data);
	//����ڵ����ݵ��ļ���һά��ʽ��,��ʽ�ɱ�tecplot���룬����Ϊ��ͷָ�룬ȫ��������·������xyz��������(���㿪ʼ)������һ��С��0ʱ������ȫ���;
	void outputdata_1d(TriDNodeHeadDefine *nodehead,CString FullName,int NoX,int NoY,int NoZ);
	//����ά������ƽ��ķ�ʽ������ļ�����ά��ʽ��,��ʽ�ɱ�tecplot���룬����Ϊ��ͷָ�룬ȫ��������·������ƽ�棨�硰XY������ZY���ȣ����ڼ�ƽ��,�ڼ�������;
	void outputdata_2d(TriDNodeHeadDefine *nodehead,CString FullName,CString Plane,int NoPlane,int NoData);
	void outputdata_2d2(TriDNodeHeadDefine *nodehead,CString FullName,CString Plane,int NoPlane);
	//��������Ϊһά�������������������ַ����硰full��;
	void outputdata_1d(TriDNodeHeadDefine *nodehead,CString FullName,CString DirectionTemp,int NoY,int NoZ);
	void outputdata_1d(TriDNodeHeadDefine *nodehead,CString FullName,int NoX,CString DirectionTemp,int NoZ);
	void outputdata_1d(TriDNodeHeadDefine *nodehead,CString FullName,int NoX,int NoY,CString DirectionTemp);
	void outputdata_Array(CString FullName,double x[],double y[],int N);
	void outputdata_Array(CString FullName,double y[],int N);

	void CreatData(int NuX, int NuY, int NuZ);
	void DeleteData();
	void ReadDataFromFile(CString FileFullName);
	//���ļ�FullName��ȡ���ݣ����ļ���ʽͬtecplot��ȡ��ʽ;
	TriDNodeHeadDefine * ReadData(CString FullName);
	//���ַ�����linetext�ֽ�double����������datatemp������ʶ�����ݵĸ���;
	//��VARIABLES�洢���ַ�������ֵ����;
	void ATOF_VARIABLES(TriDNodeHeadDefine *nodehead);
	int IdentifyData(char linetext[LineChar_MaxLength],double datatemp[]);
	//�����ļ���һ�д�������ȡ�ļ���;
	CString Process1stLine(char linetext[LineChar_MaxLength]);
	//�����ļ��ڶ��д�����ȡ������;
	int Process2ndLine(char linetext[LineChar_MaxLength],CString VARIABLESTemp[],CString VARIABLESDataTemp[]);
	//�����ļ������д�����ȡIJK��Ϣ;
	bool Process3rdLine(char linetext[LineChar_MaxLength]);
	//��ȡ����ʱ�������¶�ȡ�ڵ�����һ�ڵ��б�洢��xyz����;
	int ProcessDirection(TriDNodeDefine *nodenew,TriDNodeDefine *nodeold);
	//��node1��data�����ݼ����꿽����node2;
	bool CopyDataNode1toNode2(TriDNodeDefine *node1,TriDNodeDefine *node2);
	//Ѱ�ҽڵ㣬IJKΪ�ڵ��ţ����㿪ʼ;
	TriDNodeDefine * FindNode(TriDNodeHeadDefine *nodehead,int I, int J, int K);
	//��������Խڵ㶨λ��XYZΪ���������ϵ����꣬�ҵ��õ���Χ��8�����㣬�������굥������;
	void LocationNode(TriDNodeHeadDefine *nodehead,double X, double Y, double Z,TriDNodeDefine *node[8]);
	//�����ڵ�֮�����Բ�ֵ��һ���ڵ�,��������һ���ڵ㣬�ڶ����ڵ㣬��ֵ�����꣬����1,2,3����xyz����dataʹ�����ݵĸ���(����֪������д����999);
	TriDNodeDefine * InsertOneNodeBetweenTwoNodes(TriDNodeDefine *VertexNode1,TriDNodeDefine *VertexNode2,double xAxis,int Direction,int NuData);
	//����x,y,z���꣬����ά���ݽṹ�ڶ����ڵ�8���ڵ�������Բ�ֵ����ȡ���ݣ�������8���ڵ�ָ�룬xyz���꣬��ֵ������ݴ���ռ䣬���ݸ���;
	void InterpolationFrom8Nodes(TriDNodeDefine *VertexNode[8],double x, double y, double z,double data[],int NuData);
	//����x,y,z���꣬����ά���ݿռ���Ѱ�����ڽڵ㣬���������Բ�ֵ����ȡ���ݣ�������data������;
	bool InterpolationFromNodes(TriDNodeHeadDefine *nodehead,double x, double y, double z,double data[]);
	//��һ���������У��ҵ�x���ڵ�λ�ã������i��i+1֮�䣬����i�����Ϲ涨����-2������߷���-1,���ұ߷���NYT-1;
	int LocationInLineArray(double xAxis[],int NYT,double x);
	//���Բ�ֵ,��������x����������y��������������ֵx���飬��ֵ�洢�Ŀռ�y����ֵ����;
	double LineInterpolation(double xAxis[], double yAxis[], int NYT, double x);
	//���Բ�ֵ,��������x����������y��������������ֵx�����ز�ֵ����;
	bool LineInterpolation(double xAxis[], double yAxis[], int NYT, double x[],double y[],int Nyt);
	//���Բ�ֵ,����x1,����y1������x2,����y2����������,��ֵx���飬��ֵ�洢�Ŀռ�y���飬��ֵ����;
	void LineInterpolation(double xAxis1,double xAxis2, double yAxis1, double yAxis2,double x[],double y[],int Nyt);
	//���Բ�ֵ,����x1,����y1������x2,����y2����ֵx�����ز�ֵ����;
	double LineInterpolation(double xAxis1,double xAxis2, double yAxis1, double yAxis2,double x);
	//������ֵ,��������x����������y��������������ֵx���飬��ֵ�洢�Ŀռ�y����ֵ����;
	bool SplineInterpolation(double xAxis[], double yAxis[], int NYT, double x[],double y[],int Nyt);
	//������ֵ,��������x����������y��������������ֵx���飬��ֵ�洢�Ŀռ�y����ֵ�������Ƿ�yֵָ���任;
	bool  SplineInterpolation(double xAxis[], double yAxis[], int NYT, double x[],double y[],int Nyt,bool IsExponent);
	//�������ղ�ֵ,��������x����������y��������������ֵx���飬��ֵ�洢�Ŀռ�y����ֵ����;
	bool LagrangeInterpolation(double xAxis[], double yAxis[], int NYT, double x[],double y[],int Nyt);
	//��һ���ռ����һ����������ݲ�ֵ����
	TriDNodeHeadDefine * ExpendDatabufInY(TriDNodeHeadDefine *head,int Beishu);

	//��B=AX;
	TriDNodeHeadDefine *BCalculation(TriDNodeHeadDefine *headA,CString DirectionA,TriDNodeHeadDefine *headX,CString DirectionX);
	//��A������B�������ڻ�;
	double ABVectorsProduct(TriDNodeDefine *nodeA,int directionA,TriDNodeDefine *nodeB,int directionB,int NoData);
	//����A��ƽ����ƽ��;
	double VectorsMeanSquare(TriDNodeDefine *nodeA,int directionA,int NoData);
	//�󷽳̲в�,RMS(B-AX),ϵ������A��ͷ�ڵ㣬�����ұ�Bͷ�ڵ㣬δ֪��Xͷ�ڵ㣬����в�ĵ�NoData������;
	bool ResidualOfABX(TriDNodeHeadDefine *headA,int directionAI,int directionAJ,TriDNodeHeadDefine *headB,int directionB,TriDNodeHeadDefine *headX,int directionX,double Resid[]);
	//��������,Sirt,ϵ������A����������ƽ�棨�硰XY��ָXΪ����YΪ���򣩣��������B������δ֪�����ݿռ�X������,�洢�в������;
	bool MLEM_Iteration(TriDNodeHeadDefine *headA,CString DirectionA,TriDNodeHeadDefine *headB,CString DirectionB,TriDNodeHeadDefine *headX,CString DirectionX,double Resid[]);
	//��������,ART,ϵ������A����������ƽ�棨�硰XY��ָXΪ����YΪ���򣩣��������B������δ֪�����ݿռ�X������,�洢�в������;
	bool ART_Iteration(TriDNodeHeadDefine *headA,CString DirectionA,TriDNodeHeadDefine *headB,CString DirectionB,TriDNodeHeadDefine *headX,CString DirectionX,double Resid[]);
	bool Gauss_Iteration(TriDNodeHeadDefine *headA,CString DirectionA,TriDNodeHeadDefine *headB,CString DirectionB,TriDNodeHeadDefine *headX,CString DirectionX,double  Resid[]);
	//X����ֵ,X=AT*B,ϵ������A����������ƽ�棨�硰XY��ָXΪ����YΪ���򣩣��������B������δ֪�����ݿռ�X������;
	bool InitializationX(TriDNodeHeadDefine *headA,CString DirectionA,TriDNodeHeadDefine *headB,CString DirectionB,TriDNodeHeadDefine *headX,CString DirectionX);
	//��A�ռ��ڵ����ݳ���B�ռ������ݵı�ֵ��A/B�������ر�ֵ�����ݿռ�;
	TriDNodeHeadDefine *RatioOfAToB(TriDNodeHeadDefine *headA,TriDNodeHeadDefine *headB);
	bool IterationSGS(TriDNodeHeadDefine *headA,CString DirectionA,TriDNodeHeadDefine *headB,CString DirectionB,TriDNodeHeadDefine *headX,CString DirectionX,double  Resid[]);

	//����С���˵�ART
	TriDNodeHeadDefine *headAZXEC;
	TriDNodeHeadDefine *headBZXEC;
	bool zuixiaoercheng(TriDNodeHeadDefine *headA,CString DirectionA,TriDNodeHeadDefine *headB,CString DirectionB);

	bool LST_Changed(TriDNodeHeadDefine *headA,CString DirectionA,TriDNodeHeadDefine *headB,CString DirectionB);
};


