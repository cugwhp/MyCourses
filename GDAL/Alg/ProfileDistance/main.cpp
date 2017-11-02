/****************************************************************************
* ��������������DEMդ�������ϣ������������߾���							*
* ʵ�ֲ��裺1)���������ת����ͼ������ϵ������ֵ������					*
*			 2)��������ľ��루����Ϊ��λ������һ����Ԫ��С�����߶�			*
*			 3)�ڲ��߶���ÿ��ĸ߳�ֵZ���˴���Ϊ���ڽ��ڲ壩				*
*			 4)��ÿ�����ͼ������ת��Ϊͼ�����꣬���������������ת��		*
*			 5)�ۻ������߶���������ŷʽ���룬��Ϊ�����������߳���		*
* Author: Wang Hongping														*
* Date/Time: 2017-11-12 13:17PM												*
****************************************************************************/

// ͷ�ļ�
#include <iostream>
#include <gdal_priv.h>	//gdal�⣬����Ӱ���д������任
#include <cmath>
#include <string>
#include <vector>
#include "cpl_conv.h"
#include "ogr_srs_api.h"	//����任
#include "ogr_spatialref.h"	//�ռ�ο�

using namespace std;	// namespace

// �����ά������ṹ��
struct DOT_2D
{
	double	dx;	//x����
	double	dy;	//y����
};

// ������ά������ṹ��
struct DOT_3D
{
	double	dx;	// x
	double	dy;	// y
	double	dz;	// z
};

//==========================================================================//
// ����DEM���������������												//
// �����б�GDALDatasetH hDS - DEMդ�����ݼ�								//
// 			 const DOT_2D& dot1 - ��ʼ��1									//
//			 const DOT_2D& dot2 - �յ�2										//
// ����ֵ��double - ������룬���󷵻�ֵ-1.0 < 0.0							//
//==========================================================================//
double CalcProfileDistance(GDALDatasetH hDS, const DOT_2D& dot1, const DOT_2D& dot2);


//==========================================================================//
// ��������������������ͼ�񣬵��ýӿ�										//
//==========================================================================//
int main(int argc, char** argv)
{
	if (argc < 2)
	{
		cerr << "No input Raster Datsset Path." << endl;
		return -1;
	}

	GDALAllRegister();	//ʹ��GDAL

	GDALDatasetH	hDS;	//դ������

	hDS = GDALOpen(argv[1], GA_ReadOnly);		//���ļ�
	if (hDS == NULL)
	{
		cerr << "Open Raster File Failure. " << argv[1] << endl;
		return 1;
	}

	// ��ʼ����ʼ������
	DOT_2D	dot1, dot2;

	dot1.dx = 115.1;
	dot1.dy = 30.1;

	dot2.dx = dot1.dx + 0.52;
	dot2.dy = dot1.dy + 0.31;

	// ���ú��������������߳���
	double	dDist = CalcProfileDistance(hDS, dot1, dot2);
	cout << "Profile Distance is : " << dDist << endl;

	// �ر�դ������
	if (hDS)
	{
		GDALClose(hDS);
		hDS = NULL;
	}

	return 0;
}

//==========================================================================//
// ����DEM���������������												//
// �����б�GDALDatasetH hDS - DEMդ�����ݼ�								//
// 			 const DOT_2D& dot1 - ��ʼ��1									//
//			 const DOT_2D& dot2 - �յ�2										//
// ����ֵ��double - ������룬���󷵻�ֵ-1.0 < 0.0							//
//==========================================================================//
double CalcProfileDistance(GDALDatasetH hDS, const DOT_2D& dot1, const DOT_2D& dot2)
{
	double			dProfileDist = 0.0f;	//Profile Distance

	if (hDS == NULL)
	{
		return -1.0f;	//
	}

	//--------------------------------------------------------------//
	//					 1. ����ʼ��ľ�γ�� -> ͼ������ֵ			//
	//--------------------------------------------------------------//
	double	dszTransform[6];	// ����ת������
	double	dszInvTransform[6];	// ��任����

	GDALGetGeoTransform(hDS, dszTransform);	// ��ȡ����任ϵ��
	GDALInvGeoTransform(dszTransform, dszInvTransform);	//��任ϵ��

	// ��ʼ��γ������ -> ͼ������ֵ����
	DOT_2D	dBegRasDot, dEndRasDot;	
	GDALApplyGeoTransform(dszInvTransform, dot1.dx, dot1.dy, &dBegRasDot.dx, &dBegRasDot.dy);
	GDALApplyGeoTransform(dszInvTransform, dot2.dx, dot2.dy, &dEndRasDot.dx, &dEndRasDot.dy);

	//--------------------------------------------------------------//
	// 2. ��ͼ������ϵ�£�����ʼ�㰴1����ԪΪ��࣬�ڲ��м�������	//
	//==============================================================//
	double		dx = dEndRasDot.dx - dBegRasDot.dx;		//��ʱ����x,y,z
	double		dy = dEndRasDot.dy - dBegRasDot.dy;
	double		dz;
	double		dPixelDist = sqrt(dx*dx+dy*dy);			//�߶ξ��루��λ���أ�
	int			nDotNum = (long) ceil(dPixelDist) + 1;	//�߶εĽڵ�����
	vector<DOT_3D>		vecProfileNodes;	//�����߽ڵ��б�
	DOT_3D		dot3d;

	// ������ʼ��
	dot3d.dx = dBegRasDot.dx;
	dot3d.dy = dBegRasDot.dy;
	vecProfileNodes.push_back(dot3d);

	// ��ÿ���ڵ����Ϊ1�����أ����ܵ�����
	for (int i=1; i<nDotNum-1; ++i)
	{
		// �������꣬��ӵ�����
		dot3d.dx = dBegRasDot.dx + i*(dEndRasDot.dx - dBegRasDot.dx)/dPixelDist;
		dot3d.dy = dBegRasDot.dy + i*(dEndRasDot.dy - dBegRasDot.dy)/dPixelDist;
		vecProfileNodes.push_back(dot3d);
	}

	// ����β�ڵ�
	dot3d.dx = dEndRasDot.dx;
	dot3d.dy = dEndRasDot.dy;
	vecProfileNodes.push_back(dot3d);

	//-------------------------------------------------------------------------//
	// 3. �ز����ȵ�ÿ����ĸ߳�
	//=========================================================================//
	GDALRasterBandH	hDemBand = GDALGetRasterBand(hDS, 1);	//��ȡ�̲߳���
	if (hDemBand == NULL)
		return -2.0f;	//

	for (int i=0; i<nDotNum; ++i)
	{
		// ���ü򵥵����ڽ��ڲ�߳�Z
		GDALRasterIO(hDemBand, GF_Read, int(vecProfileNodes[i].dx), int(vecProfileNodes[i].dy), 1, 1, 
			(void*)&vecProfileNodes[i].dz, 1, 1, GDT_Float64, 0, 0);

		// ͼ������ֵ -> ��������
		dx = vecProfileNodes[i].dx;
		dy = vecProfileNodes[i].dy;
		GDALApplyGeoTransform(dszTransform, dx, dy, &vecProfileNodes[i].dx, &vecProfileNodes[i].dy);
	}

	//----------------------------------------------------------------------//
	//					 4. ͶӰ�任������ÿ����ĸ�˹���� 					//
	OGRSpatialReference ogrRasSrs;		//ͼ��ռ�ο�
	OGRSpatialReference	ogrTMSrs;		//TMͶӰ
	OGRCoordinateTransformation*	pogrCT = NULL;	//����ת��

	// ��Dataset��ȡͶӰ��Ϣ�������ռ�ο�
	const char * pszProjection = GDALGetProjectionRef(hDS);
	ogrRasSrs.importFromWkt((char**)&pszProjection);
//	cout << ogrRasSrs.IsGeographic() << endl;

	// ����TMͶӰ
	ogrTMSrs.SetProjCS("TM WGS84");
	ogrTMSrs.CopyGeogCSFrom(&ogrRasSrs);
	ogrTMSrs.SetTM(0, 114, 1.0, 500000, 0);	//����λ�õ������뾭��
//	cout << ogrTMSrs.IsProjected() << endl;

	// ��������ת������
	pogrCT = OGRCreateCoordinateTransformation(&ogrRasSrs, &ogrTMSrs);
	if (pogrCT == NULL)
		return -3.0f;	//����ת������ʧ��

	for (int i=0; i<nDotNum; ++i)	//ת��ÿһ����
	{
		pogrCT->Transform(1, &vecProfileNodes[i].dx, &vecProfileNodes[i].dy, &vecProfileNodes[i].dz);
	}

	// ��������ת������
	OGRCoordinateTransformation::DestroyCT(pogrCT);
	pogrCT = NULL;

	// 5. �ֶ�ͳ�ƾ��룬�����ܵľ���
	for (int i=1; i<nDotNum; ++i)
	{
		dx = vecProfileNodes[i].dx - vecProfileNodes[i-1].dx;
		dy = vecProfileNodes[i].dy - vecProfileNodes[i-1].dy;
		dz = vecProfileNodes[i].dz - vecProfileNodes[i-1].dz;

		dProfileDist += sqrt(dx*dx + dy*dy + dz*dz);	//�����
	}

	return dProfileDist;	//���������߳���
}