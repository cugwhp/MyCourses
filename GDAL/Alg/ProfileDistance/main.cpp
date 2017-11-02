/****************************************************************************
* 功能描述：计算DEM栅格数据上，两点间的剖面线距离							*
* 实现步骤：1)将点的坐标转换到图像坐标系（行列值）处理；					*
*			 2)计算两点的距离（像素为单位），按一个像元大小加密线段			*
*			 3)内插线段上每点的高程值Z（此处简化为最邻近内插）				*
*			 4)将每个点的图像坐标转化为图形坐标，并对坐标进行坐标转换		*
*			 5)累积计算线段上两点间的欧式距离，即为两点间的剖面线长度		*
* Author: Wang Hongping														*
* Date/Time: 2017-11-12 13:17PM												*
****************************************************************************/

// 头文件
#include <iostream>
#include <gdal_priv.h>	//gdal库，用于影像读写和坐标变换
#include <cmath>
#include <string>
#include <vector>
#include "cpl_conv.h"
#include "ogr_srs_api.h"	//坐标变换
#include "ogr_spatialref.h"	//空间参考

using namespace std;	// namespace

// 定义二维点坐标结构体
struct DOT_2D
{
	double	dx;	//x坐标
	double	dy;	//y坐标
};

// 定义三维点坐标结构体
struct DOT_3D
{
	double	dx;	// x
	double	dy;	// y
	double	dz;	// z
};

//==========================================================================//
// 计算DEM上两点间的剖面距离												//
// 参数列表：GDALDatasetH hDS - DEM栅格数据集								//
// 			 const DOT_2D& dot1 - 起始点1									//
//			 const DOT_2D& dot2 - 终点2										//
// 返回值：double - 剖面距离，错误返回值-1.0 < 0.0							//
//==========================================================================//
double CalcProfileDistance(GDALDatasetH hDS, const DOT_2D& dot1, const DOT_2D& dot2);


//==========================================================================//
// 主函数，测试用例，打开图像，调用接口										//
//==========================================================================//
int main(int argc, char** argv)
{
	if (argc < 2)
	{
		cerr << "No input Raster Datsset Path." << endl;
		return -1;
	}

	GDALAllRegister();	//使用GDAL

	GDALDatasetH	hDS;	//栅格数据

	hDS = GDALOpen(argv[1], GA_ReadOnly);		//打开文件
	if (hDS == NULL)
	{
		cerr << "Open Raster File Failure. " << argv[1] << endl;
		return 1;
	}

	// 初始化起始点坐标
	DOT_2D	dot1, dot2;

	dot1.dx = 115.1;
	dot1.dy = 30.1;

	dot2.dx = dot1.dx + 0.52;
	dot2.dy = dot1.dy + 0.31;

	// 调用函数，计算剖面线长度
	double	dDist = CalcProfileDistance(hDS, dot1, dot2);
	cout << "Profile Distance is : " << dDist << endl;

	// 关闭栅格数据
	if (hDS)
	{
		GDALClose(hDS);
		hDS = NULL;
	}

	return 0;
}

//==========================================================================//
// 计算DEM上两点间的剖面距离												//
// 参数列表：GDALDatasetH hDS - DEM栅格数据集								//
// 			 const DOT_2D& dot1 - 起始点1									//
//			 const DOT_2D& dot2 - 终点2										//
// 返回值：double - 剖面距离，错误返回值-1.0 < 0.0							//
//==========================================================================//
double CalcProfileDistance(GDALDatasetH hDS, const DOT_2D& dot1, const DOT_2D& dot2)
{
	double			dProfileDist = 0.0f;	//Profile Distance

	if (hDS == NULL)
	{
		return -1.0f;	//
	}

	//--------------------------------------------------------------//
	//					 1. 将起始点的经纬度 -> 图像行列值			//
	//--------------------------------------------------------------//
	double	dszTransform[6];	// 坐标转换参数
	double	dszInvTransform[6];	// 逆变换参数

	GDALGetGeoTransform(hDS, dszTransform);	// 获取仿射变换系数
	GDALInvGeoTransform(dszTransform, dszInvTransform);	//逆变换系数

	// 起始经纬度坐标 -> 图像行列值坐标
	DOT_2D	dBegRasDot, dEndRasDot;	
	GDALApplyGeoTransform(dszInvTransform, dot1.dx, dot1.dy, &dBegRasDot.dx, &dBegRasDot.dy);
	GDALApplyGeoTransform(dszInvTransform, dot2.dx, dot2.dy, &dEndRasDot.dx, &dEndRasDot.dy);

	//--------------------------------------------------------------//
	// 2. 在图像坐标系下，将起始点按1个像元为间距，内插中间点的坐标	//
	//==============================================================//
	double		dx = dEndRasDot.dx - dBegRasDot.dx;		//临时变量x,y,z
	double		dy = dEndRasDot.dy - dBegRasDot.dy;
	double		dz;
	double		dPixelDist = sqrt(dx*dx+dy*dy);			//线段距离（单位像素）
	int			nDotNum = (long) ceil(dPixelDist) + 1;	//线段的节点数量
	vector<DOT_3D>		vecProfileNodes;	//剖面线节点列表
	DOT_3D		dot3d;

	// 插入起始点
	dot3d.dx = dBegRasDot.dx;
	dot3d.dy = dBegRasDot.dy;
	vecProfileNodes.push_back(dot3d);

	// 按每个节点距离为1个像素，加密点坐标
	for (int i=1; i<nDotNum-1; ++i)
	{
		// 加密坐标，添加到队列
		dot3d.dx = dBegRasDot.dx + i*(dEndRasDot.dx - dBegRasDot.dx)/dPixelDist;
		dot3d.dy = dBegRasDot.dy + i*(dEndRasDot.dy - dBegRasDot.dy)/dPixelDist;
		vecProfileNodes.push_back(dot3d);
	}

	// 插入尾节点
	dot3d.dx = dEndRasDot.dx;
	dot3d.dy = dEndRasDot.dy;
	vecProfileNodes.push_back(dot3d);

	//-------------------------------------------------------------------------//
	// 3. 重采样等到每个点的高程
	//=========================================================================//
	GDALRasterBandH	hDemBand = GDALGetRasterBand(hDS, 1);	//获取高程波段
	if (hDemBand == NULL)
		return -2.0f;	//

	for (int i=0; i<nDotNum; ++i)
	{
		// 采用简单的最邻近内插高程Z
		GDALRasterIO(hDemBand, GF_Read, int(vecProfileNodes[i].dx), int(vecProfileNodes[i].dy), 1, 1, 
			(void*)&vecProfileNodes[i].dz, 1, 1, GDT_Float64, 0, 0);

		// 图像行列值 -> 地理坐标
		dx = vecProfileNodes[i].dx;
		dy = vecProfileNodes[i].dy;
		GDALApplyGeoTransform(dszTransform, dx, dy, &vecProfileNodes[i].dx, &vecProfileNodes[i].dy);
	}

	//----------------------------------------------------------------------//
	//					 4. 投影变换，计算每个点的高斯坐标 					//
	OGRSpatialReference ogrRasSrs;		//图像空间参考
	OGRSpatialReference	ogrTMSrs;		//TM投影
	OGRCoordinateTransformation*	pogrCT = NULL;	//坐标转换

	// 从Dataset读取投影信息，构建空间参考
	const char * pszProjection = GDALGetProjectionRef(hDS);
	ogrRasSrs.importFromWkt((char**)&pszProjection);
//	cout << ogrRasSrs.IsGeographic() << endl;

	// 构建TM投影
	ogrTMSrs.SetProjCS("TM WGS84");
	ogrTMSrs.CopyGeogCSFrom(&ogrRasSrs);
	ogrTMSrs.SetTM(0, 114, 1.0, 500000, 0);	//根据位置调整中央经线
//	cout << ogrTMSrs.IsProjected() << endl;

	// 构建坐标转换对象
	pogrCT = OGRCreateCoordinateTransformation(&ogrRasSrs, &ogrTMSrs);
	if (pogrCT == NULL)
		return -3.0f;	//坐标转换工作失败

	for (int i=0; i<nDotNum; ++i)	//转换每一个点
	{
		pogrCT->Transform(1, &vecProfileNodes[i].dx, &vecProfileNodes[i].dy, &vecProfileNodes[i].dz);
	}

	// 销毁坐标转换对象
	OGRCoordinateTransformation::DestroyCT(pogrCT);
	pogrCT = NULL;

	// 5. 分段统计距离，计算总的距离
	for (int i=1; i<nDotNum; ++i)
	{
		dx = vecProfileNodes[i].dx - vecProfileNodes[i-1].dx;
		dy = vecProfileNodes[i].dy - vecProfileNodes[i-1].dy;
		dz = vecProfileNodes[i].dz - vecProfileNodes[i-1].dz;

		dProfileDist += sqrt(dx*dx + dy*dy + dz*dz);	//距离和
	}

	return dProfileDist;	//返回剖面线长度
}