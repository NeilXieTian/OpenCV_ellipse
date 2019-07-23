/*
This code is intended for academic use only.
You are free to use and modify the code, at your own risk.

If you use this code, or find it useful, please refer to the paper:


The comments in the code refer to the abovementioned paper.
If you need further details about the code or the algorithm, please contact me at:

lianbosong@foxmail.com

last update: 
*/

#include "CNEllipseDetector.h"


bool myselect1=true;
bool myselect2=false;
bool myselect3=false;
float tCNC=0.3f;
float tTCNl=0.05f;//过滤太直的曲线
////////////////
//5.1Xt新增Ransac
//生成[0,1]之间符合均匀分布的数
double uniformRandom(void)
{
	return (double)rand() / (double)RAND_MAX;
}

//5.12新增  自写的最小二乘法拟合曲线
/*************************************************************************
最小二乘法拟合直线，y = a*x + b; n组数据; r-相关系数[-1,1],fabs(r)->1,说明x,y之间线性关系好，fabs(r)->0，x,y之间无线性关系，拟合无意义
a = (n*C - B*D) / (n*A - B*B)
b = (A*D - B*C) / (n*A - B*B)
r = E / F
其中：
A = sum(Xi * Xi)
B = sum(Xi)
C = sum(Xi * Yi)
D = sum(Yi)
E = sum((Xi - Xmean)*(Yi - Ymean))
F = sqrt(sum((Xi - Xmean)*(Xi - Xmean))) * sqrt(sum((Yi - Ymean)*(Yi - Ymean)))
**************************************************************************/
void LineFitLeastSquares(vector<Point2f> ptsF, vector<float> &vResult)
{
	float A = 0.0;
	float B = 0.0;
	float C = 0.0;
	float D = 0.0;
	float E = 0.0;
	float F = 0.0;

	int data_n = ptsF.size();

	for (int i = 0; i < data_n; i++)
	{
		ptsF[i].x;
		ptsF[i].y;
		A += ptsF[i].x * ptsF[i].x;
		B += ptsF[i].x;
		C += ptsF[i].x * ptsF[i].y;
		D += ptsF[i].y;
	}

	// 计算斜率a和截距b
	float a, b, temp = 0;
	if (temp = (data_n*A - B * B))// 判断分母不为0
	{
		a = (data_n*C - B * D) / temp;
		b = (A*D - B * C) / temp;
	}
	else
	{
		a = 1;
		b = 0;
	}

	// 计算相关系数r
	float Xmean, Ymean;
	Xmean = B / data_n;
	Ymean = D / data_n;

	float tempSumXX = 0.0, tempSumYY = 0.0;
	for (int i = 0; i < data_n; i++)
	{
		tempSumXX += (ptsF[i].x - Xmean) * (ptsF[i].x - Xmean);
		tempSumYY += (ptsF[i].y - Ymean) * (ptsF[i].y - Ymean);
		E += (ptsF[i].x - Xmean) * (ptsF[i].y - Ymean);
	}
	F = sqrt(tempSumXX) * sqrt(tempSumYY);

	float r;
	r = E / F;

	vResult.push_back(a);
	vResult.push_back(b);
	vResult.push_back(r*r);
}

//根据点集拟合直线ax+by+c=0，res为残差
void calcLinePara(vector<Point2d> pts, double &a, double &b, double &c, double &res)
{
	res = 0;
	Vec4f line;
	vector<Point2f> ptsF;
	for (unsigned int i = 0; i < pts.size(); i++)
		ptsF.push_back(pts[i]);

	InputArray ia = pts;
	Mat points = ia.getMat();
	int npoints2 = points.checkVector(2, -1, false);
	cout << "npoints2 = " << npoints2 << endl;

	vector<float> vResult;
	LineFitLeastSquares(ptsF, vResult);
	a = vResult[0];
	b = -1;
	c = vResult[1];

	/*
	fitLine(ptsF, line, CV_DIST_L2, 0, 1e-2, 1e-2);
	a = line[1];
	b = -line[0];
	c = line[0] * line[3] - line[1] * line[2];
	*/
	for (unsigned int i = 0; i < pts.size(); i++)
	{
		double resid_ = fabs(pts[i].x * a + pts[i].y * b + c);
		res += resid_;
	}
	res /= pts.size();
}

//得到直线拟合样本，即在直线采样点集上随机选2个点
bool getSample(vector<int> set, vector<int> &sset)
{
	int i[2];
	if (set.size() > 2)
	{
		do
		{
			for (int n = 0; n < 2; n++)
				i[n] = int(uniformRandom() * (set.size() - 1));
		} while (!(i[1] != i[0]));
		for (int n = 0; n < 2; n++)
		{
			sset.push_back(i[n]);
		}
	}
	else
	{
		return false;
	}
	return true;
}

//直线样本中两随机点位置不能太近
bool verifyComposition(const vector<Point2d> pts)
{
	cv::Point2d pt1 = pts[0];
	cv::Point2d pt2 = pts[1];
	if (abs(pt1.x - pt2.x) < 5 && abs(pt1.y - pt2.y) < 5)
		return false;

	return true;
}

//RANSAC直线拟合
void fitLineRANSAC(vector<Point2f> ptSet, double &a, double &b, double &c, vector<bool> &inlierFlag)
{
	double residual_error = 2.99; //内点阈值

	bool stop_loop = false;
	int maximum = 0;  //最大内点数

	//最终内点标识及其残差
	inlierFlag = vector<bool>(ptSet.size(), false);
	vector<double> resids_(ptSet.size(), 3);
	int sample_count = 0;
	int N = 500;

	double res = 0;

	// RANSAC
	srand((unsigned int)time(NULL)); //设置随机数种子
	vector<int> ptsID;
	for (unsigned int i = 0; i < ptSet.size(); i++)
		ptsID.push_back(i);
	while (N > sample_count && !stop_loop)
	{
		vector<bool> inlierstemp;
		vector<double> residualstemp;
		vector<int> ptss;
		int inlier_count = 0;
		if (!getSample(ptsID, ptss))
		{
			stop_loop = true;
			continue;
		}

		vector<Point2d> pt_sam;
		pt_sam.push_back(ptSet[ptss[0]]);
		pt_sam.push_back(ptSet[ptss[1]]);

		if (!verifyComposition(pt_sam))
		{
			++sample_count;
			continue;
		}

		// 计算直线方程
		calcLinePara(pt_sam, a, b, c, res);
		//内点检验
		for (unsigned int i = 0; i < ptSet.size(); i++)
		{
			Point2d pt = ptSet[i];
			double resid_ = fabs(pt.x * a + pt.y * b + c);
			residualstemp.push_back(resid_);
			inlierstemp.push_back(false);
			if (resid_ < residual_error)
			{
				++inlier_count;
				inlierstemp[i] = true;
			}
		}
		// 找到最佳拟合直线
		if (inlier_count >= maximum)
		{
			maximum = inlier_count;
			resids_ = residualstemp;
			inlierFlag = inlierstemp;
		}
		// 更新RANSAC迭代次数，以及内点概率
		if (inlier_count == 0)
		{
			N = 500;
		}
		else
		{
			double epsilon = 1.0 - double(inlier_count) / (double)ptSet.size(); //野值点比例
			double p = 0.99; //所有样本中存在1个好样本的概率
			double s = 2.0;
			N = int(log(1.0 - p) / log(1.0 - pow((1.0 - epsilon), s)));
		}
		++sample_count;
	}

	//利用所有内点重新拟合直线
	vector<Point2d> pset;
	for (unsigned int i = 0; i < ptSet.size(); i++)
	{
		if (inlierFlag[i])
			pset.push_back(ptSet[i]);
	}

	calcLinePara(pset, a, b, c, res);
}

////////////////
CNEllipseDetector::CNEllipseDetector(void) : _times(6, 0.0), _timesHelper(6, 0.0)
{
	// Default Parameters Settings
	_szPreProcessingGaussKernelSize = Size(5, 5);
	_dPreProcessingGaussSigma = 1.0;
	_fThPosition = 1.0f;
	_fMaxCenterDistance = 100.0f * 0.05f;
	_fMaxCenterDistance2 = _fMaxCenterDistance * _fMaxCenterDistance;
	_iMinEdgeLength = 16;
	_fMinOrientedRectSide = 3.0f;
	_fDistanceToEllipseContour = 0.1f;
	_fMinScore =0.7f;
	_fMinReliability =0.5;
	_uNs = 16;

	srand(unsigned(time(NULL)));
}


CNEllipseDetector::~CNEllipseDetector(void)
{
}

void CNEllipseDetector::SetParameters(Size	szPreProcessingGaussKernelSize,
										 double	dPreProcessingGaussSigma,
										 float 	fThPosition,
										 float	fMaxCenterDistance,
										 int		iMinEdgeLength,
										 float	fMinOrientedRectSide,
										 float	fDistanceToEllipseContour,
										 float	fMinScore,
										 float	fMinReliability,
										 int     iNs
										 )
{
	_szPreProcessingGaussKernelSize = szPreProcessingGaussKernelSize;
	_dPreProcessingGaussSigma = dPreProcessingGaussSigma;
	_fThPosition = fThPosition;
	_fMaxCenterDistance = fMaxCenterDistance;
	_iMinEdgeLength = iMinEdgeLength;
	_fMinOrientedRectSide = fMinOrientedRectSide;
	_fDistanceToEllipseContour = fDistanceToEllipseContour;
	_fMinScore = fMinScore;
	_fMinReliability = fMinReliability;
	_uNs = iNs;

	_fMaxCenterDistance2 = _fMaxCenterDistance * _fMaxCenterDistance;

}

uint inline CNEllipseDetector::GenerateKey(uchar pair, ushort u, ushort v)
{
	return (pair << 30) + (u << 15) + v;
};

int CNEllipseDetector::FindMaxK(const int* v) const
{
	int max_val = 0;
	int max_idx = 0;
	for (int i = 0; i<ACC_R_SIZE; ++i)
	{
		(v[i] > max_val) ? max_val = v[i], max_idx = i : NULL;
	}

	return max_idx + 90;
};

int CNEllipseDetector::FindMaxN(const int* v) const
{
	int max_val = 0;
	int max_idx = 0;
	for (int i = 0; i<ACC_N_SIZE; ++i)
	{
		(v[i] > max_val) ? max_val = v[i], max_idx = i : NULL;
	}

	return max_idx;
};

int CNEllipseDetector::FindMaxA(const int* v) const
{
	int max_val = 0;
	int max_idx = 0;
	for (int i = 0; i<ACC_A_SIZE; ++i)
	{
		(v[i] > max_val) ? max_val = v[i], max_idx = i : NULL;
	}

	return max_idx;
};

float CNEllipseDetector::GetMedianSlope(vector<Point2f>& med, Point2f& M, vector<float>& slopes)
{//input med slopes ;output:M return  
	// med		: vector of points  
	// M		: centroid of the points in med  
	// slopes	: vector of the slopes  

	unsigned iNofPoints = med.size();
	//CV_Assert(iNofPoints >= 2);

	unsigned halfSize = iNofPoints >> 1;
	unsigned quarterSize = halfSize >> 1;

	vector<float> xx, yy;
	slopes.reserve(halfSize);
	xx.reserve(iNofPoints);
	yy.reserve(iNofPoints);

	for (unsigned i = 0; i < halfSize; ++i)
	{
		Point2f& p1 = med[i];
		Point2f& p2 = med[halfSize + i];

		xx.push_back(p1.x);
		xx.push_back(p2.x);
		yy.push_back(p1.y);
		yy.push_back(p2.y);

		float den = (p2.x - p1.x);
		float num = (p2.y - p1.y);

		if (den == 0) den = 0.00001f;

		slopes.push_back(num / den);
	}

	nth_element(slopes.begin(), slopes.begin() + quarterSize, slopes.end());
	nth_element(xx.begin(), xx.begin() + halfSize, xx.end());
	nth_element(yy.begin(), yy.begin() + halfSize, yy.end());
	M.x = xx[halfSize];
	M.y = yy[halfSize];

	return slopes[quarterSize];
};


void CNEllipseDetector::GetFastCenter(vector<Point>& e1, vector<Point>& e2, EllipseData& data)
{
	countsOfGetFastCenter++;
	data.isValid = true;

	unsigned size_1 = unsigned(e1.size());
	unsigned size_2 = unsigned(e2.size());

	unsigned hsize_1 = size_1 >> 1;
	unsigned hsize_2 = size_2 >> 1;

	Point& med1 = e1[hsize_1];
	Point& med2 = e2[hsize_2];

	Point2f M12, M34;
	float q2, q4;

	{// First to second Reference slope
		float dx_ref = float(e1[0].x - med2.x);
		float dy_ref = float(e1[0].y - med2.y);

		if (dy_ref == 0) dy_ref = 0.00001f;

		float m_ref = dy_ref / dx_ref;
		data.ra = m_ref;

		// Find points with same slope as reference
		vector<Point2f> med;
		med.reserve(hsize_2);

		unsigned minPoints = (_uNs < hsize_2) ? _uNs : hsize_2;//平行弦个数

		vector<uint> indexes(minPoints);
		if (_uNs < hsize_2)
		{//半个弧点数多于设定的最少平行弦数
			unsigned iSzBin = hsize_2 / unsigned(_uNs);
			unsigned iIdx = hsize_2 + (iSzBin / 2);

			for (unsigned i = 0; i<_uNs; ++i)
			{
				indexes[i] = iIdx;
				iIdx += iSzBin;
			}
		}
		else
		{
			iota(indexes.begin(), indexes.end(), hsize_2);//转成unsigned
		}
		for (uint ii = 0; ii<minPoints; ++ii)
		{//弧2中每个点找到对应平行弦点
			uint i = indexes[ii];

			float x1 = float(e2[i].x);
			float y1 = float(e2[i].y);

			uint begin = 0;
			uint end = size_1 - 1;
			//一个点与弧1上的开始后最后点在平行弦两侧
			float xb = float(e1[begin].x);
			float yb = float(e1[begin].y);
			float res_begin = ((xb - x1) * dy_ref) - ((yb - y1) * dx_ref);
			int sign_begin = sgn(res_begin);
			if (sign_begin == 0)
			{
				//found
				med.push_back(Point2f((xb + x1)* 0.5f, (yb + y1)* 0.5f));
				continue;
			}

			float xe = float(e1[end].x);
			float ye = float(e1[end].y);
			float res_end = ((xe - x1) * dy_ref) - ((ye - y1) * dx_ref);
			int sign_end = sgn(res_end);
			if (sign_end == 0)
			{
				//found
				med.push_back(Point2f((xe + x1)* 0.5f, (ye + y1)* 0.5f));
				continue;
			}
			//如果在一侧
			if ((sign_begin + sign_end) != 0)
			{//区域内没有平行的弦
				continue;
			}

			//二分法查找平行弦
			uint j = (begin + end) >> 1;
			while (end - begin > 2)
			{
				float x2 = float(e1[j].x);
				float y2 = float(e1[j].y);
				float res = ((x2 - x1) * dy_ref) - ((y2 - y1) * dx_ref);
				int sign_res = sgn(res);

				if (sign_res == 0)
				{
					//found
					med.push_back(Point2f((x2 + x1)* 0.5f, (y2 + y1)* 0.5f));
					break;
				}

				if (sign_res + sign_begin == 0)
				{
					sign_end = sign_res;
					end = j;
				}
				else
				{
					sign_begin = sign_res;
					begin = j;
				}
				j = (begin + end) >> 1;
			}
			//二分法查找平行弦 end 有错
			med.push_back(Point2f((e1[j].x + x1)* 0.5f, (e1[j].y + y1)* 0.5f));
		}

		if (med.size() < 2)
		{
			data.isValid = false;
			return;
		}

		q2 = GetMedianSlope(med, M12, data.Sa);//得到Sa ta=q2 Ma
		//cout << "************************************************" << endl;
		//cout << "1-2原方法得到的中点连线斜率为" << q2 << endl;
		double A, B, C;
		vector<bool> inliers;
		fitLineRANSAC(med, A, B, C, inliers);
		//cout << "RANSAC得到的A和B为：" << A << "	" << B << endl;
		//cout << "斜率是" << -A / B << endl;
		//cout << "************************************************" << endl;
	}

	{// Second to first
		// Reference slope
		float dx_ref = float(med1.x - e2[0].x);
		float dy_ref = float(med1.y - e2[0].y);

		if (dy_ref == 0) dy_ref = 0.00001f;

		float m_ref = dy_ref / dx_ref;
		data.rb = m_ref;

		// Find points with same slope as reference
		vector<Point2f> med;
		med.reserve(hsize_1);

		uint minPoints = (_uNs < hsize_1) ? _uNs : hsize_1;

		vector<uint> indexes(minPoints);
		if (_uNs < hsize_1)
		{
			unsigned iSzBin = hsize_1 / unsigned(_uNs);
			unsigned iIdx = hsize_1 + (iSzBin / 2);

			for (unsigned i = 0; i<_uNs; ++i)
			{
				indexes[i] = iIdx;
				iIdx += iSzBin;
			}
		}
		else
		{
			iota(indexes.begin(), indexes.end(), hsize_1);
		}


		for (uint ii = 0; ii<minPoints; ++ii)
		{
			uint i = indexes[ii];

			float x1 = float(e1[i].x);
			float y1 = float(e1[i].y);

			uint begin = 0;
			uint end = size_2 - 1;

			float xb = float(e2[begin].x);
			float yb = float(e2[begin].y);
			float res_begin = ((xb - x1) * dy_ref) - ((yb - y1) * dx_ref);
			int sign_begin = sgn(res_begin);
			if (sign_begin == 0)
			{
				//found
				med.push_back(Point2f((xb + x1)* 0.5f, (yb + y1)* 0.5f));
				continue;
			}

			float xe = float(e2[end].x);
			float ye = float(e2[end].y);
			float res_end = ((xe - x1) * dy_ref) - ((ye - y1) * dx_ref);
			int sign_end = sgn(res_end);
			if (sign_end == 0)
			{
				//found
				med.push_back(Point2f((xe + x1)* 0.5f, (ye + y1)* 0.5f));
				continue;
			}

			if ((sign_begin + sign_end) != 0)
			{
				continue;
			}

			uint j = (begin + end) >> 1;

			while (end - begin > 2)
			{
				float x2 = float(e2[j].x);
				float y2 = float(e2[j].y);
				float res = ((x2 - x1) * dy_ref) - ((y2 - y1) * dx_ref);
				int sign_res = sgn(res);

				if (sign_res == 0)
				{
					//found
					med.push_back(Point2f((x2 + x1)* 0.5f, (y2 + y1)* 0.5f));
					break;
				}

				if (sign_res + sign_begin == 0)
				{
					sign_end = sign_res;
					end = j;
				}
				else
				{
					sign_begin = sign_res;
					begin = j;
				}
				j = (begin + end) >> 1;
			}

			med.push_back(Point2f((e2[j].x + x1)* 0.5f, (e2[j].y + y1)* 0.5f));
		}

		if (med.size() < 2)
		{
			data.isValid = false;
			return;
		}
		q4 = GetMedianSlope(med, M34, data.Sb);
		//cout << "************************************************" << endl;
		//cout << "1-2原方法得到的中点连线斜率为" << q4 << endl;
		double A, B, C;
		vector<bool> inliers;
		fitLineRANSAC(med, A, B, C, inliers);
		//cout << "RANSAC得到的A和B为：" << A << "	" << B << endl;
		//cout << "斜率是" << -A / B << endl;
		//cout << "************************************************" << endl;
	}

	if (q2 == q4)
	{
		data.isValid = false;
		return;
	}

	float invDen = 1 / (q2 - q4);
	data.Cab.x = (M34.y - q4*M34.x - M12.y + q2*M12.x) * invDen;
	data.Cab.y = (q2*M34.y - q4*M12.y + q2*q4*(M12.x - M34.x)) * invDen;
	data.ta = q2;
	data.tb = q4;
	data.Ma = M12;
	data.Mb = M34;
};

#define	DISCARD_TCN1
void CNEllipseDetector::DetectEdges13(Mat1b& DP, VVP& points_1, VVP& points_3)
{
	// Vector of connected edge points
	VVP contours;
	int countedges=0;
	// Labeling 8-connected edge points, discarding edge too small
	Labeling(DP, contours, _iMinEdgeLength);//把在一个弧的点放到一个无序集合中
	int iContoursSize = int(contours.size());

	// For each edge
	for (int i = 0; i < iContoursSize; ++i)
	{
		VP& edgeSegment = contours[i];
#ifndef DISCARD_CONSTRAINT_OBOX

		// Selection strategy - Step 1 - See Sect [3.1.2] of the paper
		// Constraint on axes aspect ratio
		RotatedRect oriented = minAreaRect(edgeSegment);
		float o_min = min(oriented.size.width, oriented.size.height);

		if (o_min < _fMinOrientedRectSide)
		{
			countedges++;
			continue;
		}
#endif
		// Order edge points of the same arc
		sort(edgeSegment.begin(), edgeSegment.end(), SortTopLeft2BottomRight);//有序连接点
		int iEdgeSegmentSize = unsigned(edgeSegment.size());

		// Get extrema of the arc
		Point& left = edgeSegment[0];
		Point& right = edgeSegment[iEdgeSegmentSize - 1];
#ifndef DISCARD_TCN
#ifndef DISCARD_TCN2
		int flag=0;
		for(int j=0;j<iEdgeSegmentSize;j++){
			Point& mid = edgeSegment[j];
			float data[] = {left.x, left.y, 1,mid.x, mid.y, 1,right.x, right.y, 1};
			Mat threePoints(Size(3, 3), CV_32FC1, data);
			double ans = determinant(threePoints);

			float dx=1.0f*(left.x-right.x);
			float dy=1.0f*(left.y-right.y);
			float edgelength2=dx*dx+dy*dy;
			//double TCNl=ans/edgelength2;
			double TCNl=ans/(2*sqrt(edgelength2));
			if (abs(TCNl)>tTCNl){
				flag=1;
				break;
			}
		}
		if(0==flag){
			countedges++;
			continue;
		}
#endif
#ifndef DISCARD_TCN1
		Point& mid = edgeSegment[iEdgeSegmentSize/2];
		float data[] = {left.x, left.y, 1,mid.x, mid.y, 1,right.x, right.y, 1};
		Mat threePoints(Size(3, 3), CV_32FC1, data);
		double ans = determinant(threePoints);

		float dx=1.0f*(left.x-right.x);
		float dy=1.0f*(left.y-right.y);
		float edgelength2=dx*dx+dy*dy;
		double TCNl=ans/edgelength2;
		//double TCNl=ans/(2*sqrt(edgelength2));
		if (abs(TCNl)<tTCNl){
			countedges++;
			continue;
		}
#endif
#endif
		// Find convexity - See Sect [3.1.3] of the paper
		int iCountTop = 0;
		int xx = left.x;
		for (int k = 1; k < iEdgeSegmentSize; ++k)
		{
			if (edgeSegment[k].x == xx) continue;

			iCountTop += (edgeSegment[k].y - left.y);
			xx = edgeSegment[k].x;
		}

		int width = abs(right.x - left.x) + 1;
		int height = abs(right.y - left.y) + 1;
		int iCountBottom = (width * height) - iEdgeSegmentSize - iCountTop;

		if (iCountBottom > iCountTop)
		{	//1
			points_1.push_back(edgeSegment);
		}
		else if (iCountBottom < iCountTop)
		{	//3
			points_3.push_back(edgeSegment);
		}
	}
};


void CNEllipseDetector::DetectEdges24(Mat1b& DN, VVP& points_2, VVP& points_4 )
{
	// Vector of connected edge points
	VVP contours;
	int countedges=0;
	/// Labeling 8-connected edge points, discarding edge too small
	Labeling(DN, contours, _iMinEdgeLength);

	int iContoursSize = unsigned(contours.size());

	// For each edge
	for (int i = 0; i < iContoursSize; ++i)
	{
		VP& edgeSegment = contours[i];

#ifndef DISCARD_CONSTRAINT_OBOX

		// Selection strategy - Step 1 - See Sect [3.1.2] of the paper
		// Constraint on axes aspect ratio
		RotatedRect oriented = minAreaRect(edgeSegment);
		float o_min = min(oriented.size.width, oriented.size.height);

		if (o_min < _fMinOrientedRectSide)
		{
			countedges++;
			continue;
		}
#endif
		// Order edge points of the same arc
		sort(edgeSegment.begin(), edgeSegment.end(), SortBottomLeft2TopRight);
		int iEdgeSegmentSize = unsigned(edgeSegment.size());

		// Get extrema of the arc
		Point& left = edgeSegment[0];
		Point& right = edgeSegment[iEdgeSegmentSize - 1];
#ifndef DISCARD_TCN
#ifndef DISCARD_TCN2
		int flag=0;
		for(int j=0;j<iEdgeSegmentSize;j++){
			Point& mid = edgeSegment[j];
			float data[] = {left.x, left.y, 1,mid.x, mid.y, 1,right.x, right.y, 1};
			Mat threePoints(Size(3, 3), CV_32FC1, data);
			double ans = determinant(threePoints);

			float dx=1.0f*(left.x-right.x);
			float dy=1.0f*(left.y-right.y);
			float edgelength2=dx*dx+dy*dy;
			//double TCNl=ans/edgelength2;
			double TCNl=ans/(2*sqrt(edgelength2));
			if (abs(TCNl)>tTCNl){
				flag=1;
				break;
			}
		}
		if(0==flag){
			countedges++;
			continue;
		}else{
		}
#endif
#ifndef DISCARD_TCN1
		Point& mid = edgeSegment[iEdgeSegmentSize/2];
		float data[] = {left.x, left.y, 1,mid.x, mid.y, 1,right.x, right.y, 1};
		Mat threePoints(Size(3, 3), CV_32FC1, data);
		double ans = determinant(threePoints);

		float dx=1.0f*(left.x-right.x);
		float dy=1.0f*(left.y-right.y);
		float edgelength2=dx*dx+dy*dy;
		double TCNl=ans/edgelength2;
		//double TCNl=ans/(2*sqrt(edgelength2));
		if (abs(TCNl)<tTCNl){
			countedges++;
			continue;
		}
#endif
#endif
		// Find convexity - See Sect [3.1.3] of the paper
		//判断凹凸性
		int iCountBottom = 0;
		int xx = left.x;
		for (int k = 1; k < iEdgeSegmentSize; ++k)
		{
			if (edgeSegment[k].x == xx) continue;

			iCountBottom += (left.y - edgeSegment[k].y);
			xx = edgeSegment[k].x;
		}

		int width = abs(right.x - left.x) + 1;
		int height = abs(right.y - left.y) + 1;
		int iCountTop = (width *height) - iEdgeSegmentSize - iCountBottom;

		if (iCountBottom > iCountTop)
		{
			//2
			points_2.push_back(edgeSegment);
		}
		else if (iCountBottom < iCountTop)
		{
			//4
			points_4.push_back(edgeSegment);
		}
	}
};

void CNEllipseDetector::DetectEdges13_zrk(Mat1b& DP, VVP& points_1, VVP& points_3)
{
	// Vector of connected edge points
	VVP contours;
	int countedges = 0;
	// Labeling 8-connected edge points, discarding edge too small
	Labeling_zrk(DP, contours, _iMinEdgeLength);//把在一个弧的点放到一个无序集合中
	int iContoursSize = int(contours.size());

	// For each edge
	for (int i = 0; i < iContoursSize; ++i)
	{
		VP& edgeSegment = contours[i];
#ifndef DISCARD_CONSTRAINT_OBOX

		// Selection strategy - Step 1 - See Sect [3.1.2] of the paper
		// Constraint on axes aspect ratio
		RotatedRect oriented = minAreaRect(edgeSegment);
		float o_min = min(oriented.size.width, oriented.size.height);

		if (o_min < _fMinOrientedRectSide)
		{
			countedges++;
			continue;
		}
#endif
		// Order edge points of the same arc
		sort(edgeSegment.begin(), edgeSegment.end(), SortTopLeft2BottomRight);//有序连接点
		int iEdgeSegmentSize = unsigned(edgeSegment.size());

		// Get extrema of the arc
		Point& left = edgeSegment[0];
		Point& right = edgeSegment[iEdgeSegmentSize - 1];
#ifndef DISCARD_TCN
#ifndef DISCARD_TCN2
		int flag = 0;
		for (int j = 0; j < iEdgeSegmentSize; j++) {
			Point& mid = edgeSegment[j];
			float data[] = { left.x, left.y, 1,mid.x, mid.y, 1,right.x, right.y, 1 };
			Mat threePoints(Size(3, 3), CV_32FC1, data);
			double ans = determinant(threePoints);

			float dx = 1.0f*(left.x - right.x);
			float dy = 1.0f*(left.y - right.y);
			float edgelength2 = dx * dx + dy * dy;
			//double TCNl=ans/edgelength2;
			double TCNl = ans / (2 * sqrt(edgelength2));
			if (abs(TCNl) > tTCNl) {
				flag = 1;
				break;
			}
		}
		if (0 == flag) {
			countedges++;
			continue;
		}
#endif
#ifndef DISCARD_TCN1
		Point& mid = edgeSegment[iEdgeSegmentSize / 2];
		float data[] = { left.x, left.y, 1,mid.x, mid.y, 1,right.x, right.y, 1 };
		Mat threePoints(Size(3, 3), CV_32FC1, data);
		double ans = determinant(threePoints);

		float dx = 1.0f*(left.x - right.x);
		float dy = 1.0f*(left.y - right.y);
		float edgelength2 = dx * dx + dy * dy;
		double TCNl = ans / edgelength2;
		//double TCNl=ans/(2*sqrt(edgelength2));
		if (abs(TCNl) < tTCNl) {
			countedges++;
			continue;
		}
#endif
#endif
		// Find convexity - See Sect [3.1.3] of the paper
		int iCountTop = 0;
		int xx = left.x;
		for (int k = 1; k < iEdgeSegmentSize; ++k)
		{
			if (edgeSegment[k].x == xx) continue;

			iCountTop += (edgeSegment[k].y - left.y);
			xx = edgeSegment[k].x;
		}

		int width = abs(right.x - left.x) + 1;
		int height = abs(right.y - left.y) + 1;
		int iCountBottom = (width * height) - iEdgeSegmentSize - iCountTop;

		if (iCountBottom > iCountTop)
		{	//1
			points_1.push_back(edgeSegment);
		}
		else if (iCountBottom < iCountTop)
		{	//3
			points_3.push_back(edgeSegment);
		}
	}
};

void CNEllipseDetector::DetectEdges24_zrk(Mat1b& DN, VVP& points_2, VVP& points_4)
{
	// Vector of connected edge points
	VVP contours;
	int countedges = 0;
	///Labeling 8-connected edge points, discarding edge too small
	Labeling_zrk(DN, contours, _iMinEdgeLength);

	int iContoursSize = unsigned(contours.size());

	// For each edge
	for (int i = 0; i < iContoursSize; ++i)
	{
		VP& edgeSegment = contours[i];

#ifndef DISCARD_CONSTRAINT_OBOX

		// Selection strategy - Step 1 - See Sect [3.1.2] of the paper
		// Constraint on axes aspect ratio
		RotatedRect oriented = minAreaRect(edgeSegment);
		float o_min = min(oriented.size.width, oriented.size.height);

		if (o_min < _fMinOrientedRectSide)
		{
			countedges++;
			continue;
		}
#endif
		// Order edge points of the same arc
		sort(edgeSegment.begin(), edgeSegment.end(), SortBottomLeft2TopRight);
		int iEdgeSegmentSize = unsigned(edgeSegment.size());

		// Get extrema of the arc
		Point& left = edgeSegment[0];
		Point& right = edgeSegment[iEdgeSegmentSize - 1];
#ifndef DISCARD_TCN
#ifndef DISCARD_TCN2
		int flag = 0;
		for (int j = 0; j < iEdgeSegmentSize; j++) {
			Point& mid = edgeSegment[j];
			float data[] = { left.x, left.y, 1,mid.x, mid.y, 1,right.x, right.y, 1 };
			Mat threePoints(Size(3, 3), CV_32FC1, data);
			double ans = determinant(threePoints);

			float dx = 1.0f*(left.x - right.x);
			float dy = 1.0f*(left.y - right.y);
			float edgelength2 = dx * dx + dy * dy;
			//double TCNl=ans/edgelength2;
			double TCNl = ans / (2 * sqrt(edgelength2));
			if (abs(TCNl) > tTCNl) {
				flag = 1;
				break;
			}
		}
		if (0 == flag) {
			countedges++;
			continue;
		}
		else {
		}
#endif
#ifndef DISCARD_TCN1
		Point& mid = edgeSegment[iEdgeSegmentSize / 2];
		float data[] = { left.x, left.y, 1,mid.x, mid.y, 1,right.x, right.y, 1 };
		Mat threePoints(Size(3, 3), CV_32FC1, data);
		double ans = determinant(threePoints);

		float dx = 1.0f*(left.x - right.x);
		float dy = 1.0f*(left.y - right.y);
		float edgelength2 = dx * dx + dy * dy;
		double TCNl = ans / edgelength2;
		//double TCNl=ans/(2*sqrt(edgelength2));
		if (abs(TCNl) < tTCNl) {
			countedges++;
			continue;
		}
#endif
#endif
		// Find convexity - See Sect [3.1.3] of the paper
		//判断凹凸性
		int iCountBottom = 0;
		int xx = left.x;
		for (int k = 1; k < iEdgeSegmentSize; ++k)
		{
			if (edgeSegment[k].x == xx) continue;

			iCountBottom += (left.y - edgeSegment[k].y);
			xx = edgeSegment[k].x;
		}

		int width = abs(right.x - left.x) + 1;
		int height = abs(right.y - left.y) + 1;
		int iCountTop = (width *height) - iEdgeSegmentSize - iCountBottom;

		if (iCountBottom > iCountTop)
		{
			//2
			points_2.push_back(edgeSegment);
		}
		else if (iCountBottom < iCountTop)
		{
			//4
			points_4.push_back(edgeSegment);
		}
	}
};

// Most important function for detecting ellipses. See Sect[3.2.3] of the paper
void CNEllipseDetector::FindEllipses(	Point2f& center,
										VP& edge_i,VP& edge_j,VP& edge_k,
										EllipseData& data_ij,EllipseData& data_ik,
										vector<Ellipse>& ellipses)
{//根据三段弧计算椭圆参数，并计算得到的椭圆参数的准确定。
	countsOfFindEllipse++;
	// Find ellipse parameters

	// 0-initialize accumulators
	memset(accN, 0, sizeof(int)*ACC_N_SIZE);
	memset(accR, 0, sizeof(int)*ACC_R_SIZE);
	memset(accA, 0, sizeof(int)*ACC_A_SIZE);

	Tac(3); //estimation

	// Get size of the 4 vectors of slopes (2 pairs of arcs)
	int sz_ij1 = int(data_ij.Sa.size());
	int sz_ij2 = int(data_ij.Sb.size());
	int sz_ik1 = int(data_ik.Sa.size());
	int sz_ik2 = int(data_ik.Sb.size());

	// Get the size of the 3 arcs
	size_t sz_ei = edge_i.size();
	size_t sz_ej = edge_j.size();
	size_t sz_ek = edge_k.size();

	// Center of the estimated ellipse
	float a0 = center.x;
	float b0 = center.y;


	// Estimation of remaining parameters
	// Uses 4 combinations of parameters. See Table 1 and Sect [3.2.3] of the paper.
	//ij1 and ik
	{
		float q1 = data_ij.ra;
		float q3 = data_ik.ra;
		float q5 = data_ik.rb;

		for (int ij1 = 0; ij1 < sz_ij1; ++ij1)
		{
			float q2 = data_ij.Sa[ij1];//有必要遍历吗？

			float q1xq2 = q1*q2;
			// ij1 and ik1
			for (int ik1 = 0; ik1 < sz_ik1; ++ik1)
			{
				float q4 = data_ik.Sa[ik1];//有必要遍历吗？

				float q3xq4 = q3*q4;

				// See Eq. [13-18] in the paper

				float a = (q1xq2 - q3xq4);//gama
				float b = (q3xq4 + 1)*(q1 + q2) - (q1xq2 + 1)*(q3 + q4);//beta
				float Kp = (-b + sqrt(b*b + 4 * a*a)) / (2 * a);//K+
				float zplus = ((q1 - Kp)*(q2 - Kp)) / ((1 + q1*Kp)*(1 + q2*Kp));
				// 判断有效性  zplus 和K不是线性关系
				if (zplus >= 0.0f) continue;

				float Np = sqrt(-zplus);//N+
				float rho = atan(Kp);//rho tmp
				int rhoDeg;
				if (Np > 1.f)
				{
					Np = 1.f / Np;
					rhoDeg = cvRound((rho * 180 / CV_PI) + 180) % 180; // [0,180)					
				}
				else
				{
					rhoDeg = cvRound((rho * 180 / CV_PI) + 90) % 180; // [0,180)//rho 角度表示，归一化
				}

				int iNp = cvRound(Np * 100); // [0, 100]取整

				if (0 <= iNp	&& iNp < ACC_N_SIZE &&
					0 <= rhoDeg	&& rhoDeg < ACC_R_SIZE
					)
				{//为什么都遍历，只遍历一个可以吗  因为 zplus 和K不是线性关系
					++accN[iNp];	// Increment N accumulator
					++accR[rhoDeg];	// Increment R accumulator
				}
			}
			// ij1 and ik2
			for (int ik2 = 0; ik2 < sz_ik2; ++ik2)
			{
				float q4 = data_ik.Sb[ik2];

				float q5xq4 = q5*q4;

				// See Eq. [13-18] in the paper

				float a = (q1xq2 - q5xq4);
				float b = (q5xq4 + 1)*(q1 + q2) - (q1xq2 + 1)*(q5 + q4);
				float Kp = (-b + sqrt(b*b + 4 * a*a)) / (2 * a);
				float zplus = ((q1 - Kp)*(q2 - Kp)) / ((1 + q1*Kp)*(1 + q2*Kp));

				if (zplus >= 0.0f)
				{
					continue;
				}

				float Np = sqrt(-zplus);
				float rho = atan(Kp);
				int rhoDeg;
				if (Np > 1.f)
				{
					Np = 1.f / Np;
					rhoDeg = cvRound((rho * 180 / CV_PI) + 180) % 180; // [0,180)					
				}
				else
				{
					rhoDeg = cvRound((rho * 180 / CV_PI) + 90) % 180; // [0,180)
				}

				int iNp = cvRound(Np * 100); // [0, 100]

				if (0 <= iNp	&& iNp < ACC_N_SIZE &&
					0 <= rhoDeg	&& rhoDeg < ACC_R_SIZE
					)
				{
					++accN[iNp];		// Increment N accumulator
					++accR[rhoDeg];		// Increment R accumulator
				}
			}

		}
	}

	//ij2 and ik
	{
		float q1 = data_ij.rb;
		float q3 = data_ik.rb;
		float q5 = data_ik.ra;

		for (int ij2 = 0; ij2 < sz_ij2; ++ij2)
		{
			float q2 = data_ij.Sb[ij2];

			float q1xq2 = q1*q2;
			//ij2 and ik2
			for (int ik2 = 0; ik2 < sz_ik2; ++ik2)
			{
				float q4 = data_ik.Sb[ik2];

				float q3xq4 = q3*q4;

				// See Eq. [13-18] in the paper

				float a = (q1xq2 - q3xq4);
				float b = (q3xq4 + 1)*(q1 + q2) - (q1xq2 + 1)*(q3 + q4);
				float Kp = (-b + sqrt(b*b + 4 * a*a)) / (2 * a);
				float zplus = ((q1 - Kp)*(q2 - Kp)) / ((1 + q1*Kp)*(1 + q2*Kp));

				if (zplus >= 0.0f)
				{
					continue;
				}

				float Np = sqrt(-zplus);
				float rho = atan(Kp);
				int rhoDeg;
				if (Np > 1.f)
				{
					Np = 1.f / Np;
					rhoDeg = cvRound((rho * 180 / CV_PI) + 180) % 180; // [0,180)
				}
				else
				{
					rhoDeg = cvRound((rho * 180 / CV_PI) + 90) % 180; // [0,180)
				}

				int iNp = cvRound(Np * 100); // [0, 100]

				if (0 <= iNp	&& iNp < ACC_N_SIZE &&
					0 <= rhoDeg	&& rhoDeg < ACC_R_SIZE
					)
				{
					++accN[iNp];		// Increment N accumulator
					++accR[rhoDeg];		// Increment R accumulator
				}
			}

			//ij2 and ik1
			for (int ik1 = 0; ik1 < sz_ik1; ++ik1)
			{
				float q4 = data_ik.Sa[ik1];

				float q5xq4 = q5*q4;

				// See Eq. [13-18] in the paper

				float a = (q1xq2 - q5xq4);
				float b = (q5xq4 + 1)*(q1 + q2) - (q1xq2 + 1)*(q5 + q4);
				float Kp = (-b + sqrt(b*b + 4 * a*a)) / (2 * a);
				float zplus = ((q1 - Kp)*(q2 - Kp)) / ((1 + q1*Kp)*(1 + q2*Kp));

				if (zplus >= 0.0f)
				{
					continue;
				}

				float Np = sqrt(-zplus);
				float rho = atan(Kp);
				int rhoDeg;
				if (Np > 1.f)
				{
					Np = 1.f / Np;
					rhoDeg = cvRound((rho * 180 / CV_PI) + 180) % 180; // [0,180)
				}
				else
				{
					rhoDeg = cvRound((rho * 180 / CV_PI) + 90) % 180; // [0,180)
				}

				int iNp = cvRound(Np * 100); // [0, 100]

				if (0 <= iNp	&& iNp < ACC_N_SIZE &&
					0 <= rhoDeg	&& rhoDeg < ACC_R_SIZE
					)
				{
					++accN[iNp];		// Increment N accumulator
					++accR[rhoDeg];		// Increment R accumulator
				}
			}

		}
	}

	// Find peak in N and K accumulator
	int iN = FindMaxN(accN);
	int iK = FindMaxK(accR);

	// Recover real values
	float fK = float(iK);
	float Np = float(iN) * 0.01f;
	float rho = fK * float(CV_PI) / 180.f;	//deg 2 rad
	float Kp = tan(rho);

	// Estimate A. See Eq. [19 - 22] in Sect [3.2.3] of the paper  
	// 三段弧上的点都用来估算A
	//可以优化   
	for (ushort l = 0; l < sz_ei; ++l)
	{
		Point& pp = edge_i[l];
		float sk = 1.f / sqrt(Kp*Kp + 1.f);//cos rho
		float x0 = ((pp.x - a0) * sk) + (((pp.y - b0)*Kp) * sk);//可优化
		float y0 = -(((pp.x - a0) * Kp) * sk) + ((pp.y - b0) * sk);//可优化
		float Ax = sqrt((x0*x0*Np*Np + y0*y0) / ((Np*Np)*(1.f + Kp*Kp)));
		int A = cvRound(abs(Ax / cos(rho)));//可以优化
		if ((0 <= A) && (A < ACC_A_SIZE))
		{
			++accA[A];
		}
	}

	for (ushort l = 0; l < sz_ej; ++l)
	{
		Point& pp = edge_j[l];
		float sk = 1.f / sqrt(Kp*Kp + 1.f);
		float x0 = ((pp.x - a0) * sk) + (((pp.y - b0)*Kp) * sk);
		float y0 = -(((pp.x - a0) * Kp) * sk) + ((pp.y - b0) * sk);
		float Ax = sqrt((x0*x0*Np*Np + y0*y0) / ((Np*Np)*(1.f + Kp*Kp)));
		int A = cvRound(abs(Ax / cos(rho)));
		if ((0 <= A) && (A < ACC_A_SIZE))
		{
			++accA[A];
		}
	}

	for (ushort l = 0; l < sz_ek; ++l)
	{
		Point& pp = edge_k[l];
		float sk = 1.f / sqrt(Kp*Kp + 1.f);
		float x0 = ((pp.x - a0) * sk) + (((pp.y - b0)*Kp) * sk);
		float y0 = -(((pp.x - a0) * Kp) * sk) + ((pp.y - b0) * sk);
		float Ax = sqrt((x0*x0*Np*Np + y0*y0) / ((Np*Np)*(1.f + Kp*Kp)));
		int A = cvRound(abs(Ax / cos(rho)));
		if ((0 <= A) && (A < ACC_A_SIZE))
		{
			++accA[A];
		}
	}

	// Find peak in A accumulator
	int A = FindMaxA(accA);
	float fA = float(A);

	// Find B value. See Eq [23] in the paper
	float fB = abs(fA * Np);

	// Got all ellipse parameters!
	Ellipse ell(a0, b0, fA, fB, fmod(rho + float(CV_PI)*2.f, float(CV_PI)));

	Toc(3); //estimation
	Tac(4); //validation

	// Get the score. See Sect [3.3.1] in the paper

	// Find the number of edge pixel lying on the ellipse
	float _cos = cos(-ell._rad);
	float _sin = sin(-ell._rad);

	float invA2 = 1.f / (ell._a * ell._a);
	float invB2 = 1.f / (ell._b * ell._b);

	float invNofPoints = 1.f / float(sz_ei + sz_ej + sz_ek);
	int counter_on_perimeter = 0;

	for (ushort l = 0; l < sz_ei; ++l)
	{
		float tx = float(edge_i[l].x) - ell._xc;
		float ty = float(edge_i[l].y) - ell._yc;
		float rx = (tx*_cos - ty*_sin);
		float ry = (tx*_sin + ty*_cos);

		float h = (rx*rx)*invA2 + (ry*ry)*invB2;
		if (abs(h - 1.f) < _fDistanceToEllipseContour)
		{
			++counter_on_perimeter;
		}
	}

	for (ushort l = 0; l < sz_ej; ++l)
	{
		float tx = float(edge_j[l].x) - ell._xc;
		float ty = float(edge_j[l].y) - ell._yc;
		float rx = (tx*_cos - ty*_sin);
		float ry = (tx*_sin + ty*_cos);

		float h = (rx*rx)*invA2 + (ry*ry)*invB2;
		if (abs(h - 1.f) < _fDistanceToEllipseContour)
		{
			++counter_on_perimeter;
		}
	}

	for (ushort l = 0; l < sz_ek; ++l)
	{
		float tx = float(edge_k[l].x) - ell._xc;
		float ty = float(edge_k[l].y) - ell._yc;
		float rx = (tx*_cos - ty*_sin);
		float ry = (tx*_sin + ty*_cos);

		float h = (rx*rx)*invA2 + (ry*ry)*invB2;
		if (abs(h - 1.f) < _fDistanceToEllipseContour)
		{
			++counter_on_perimeter;
		}
	}

	//no points found on the ellipse
	if (counter_on_perimeter <= 0)
	{
		Toc(4); //validation
		return;
	}


	// Compute score
	float score = float(counter_on_perimeter) * invNofPoints; 
	cout << "score " << score << endl;
	if (score < _fMinScore)
	{
		Toc(4); //validation
		return;
	}

	// Compute reliability	
	// this metric is not described in the paper, mostly due to space limitations.
	// The main idea is that for a given ellipse (TD) even if the score is high, the arcs 
	// can cover only a small amount of the contour of the estimated ellipse. 
	// A low reliability indicate that the arcs form an elliptic shape by chance, but do not underlie
	// an actual ellipse. The value is normalized between 0 and 1. 
	// The default value is 0.4.

	// It is somehow similar to the "Angular Circumreference Ratio" saliency criteria 
	// as in the paper: 
	// D. K. Prasad, M. K. Leung, S.-Y. Cho, Edge curvature and convexity
	// based ellipse detection method, Pattern Recognition 45 (2012) 3204-3221.

	float di, dj, dk;
	{
		Point2f p1(float(edge_i[0].x), float(edge_i[0].y));
		Point2f p2(float(edge_i[sz_ei - 1].x), float(edge_i[sz_ei - 1].y));
		p1.x -= ell._xc;
		p1.y -= ell._yc;
		p2.x -= ell._xc;
		p2.y -= ell._yc;
		Point2f r1((p1.x*_cos - p1.y*_sin), (p1.x*_sin + p1.y*_cos));
		Point2f r2((p2.x*_cos - p2.y*_sin), (p2.x*_sin + p2.y*_cos));
		di = abs(r2.x - r1.x) + abs(r2.y - r1.y);
	}
	{
		Point2f p1(float(edge_j[0].x), float(edge_j[0].y));
		Point2f p2(float(edge_j[sz_ej - 1].x), float(edge_j[sz_ej - 1].y));
		p1.x -= ell._xc;
		p1.y -= ell._yc;
		p2.x -= ell._xc;
		p2.y -= ell._yc;
		Point2f r1((p1.x*_cos - p1.y*_sin), (p1.x*_sin + p1.y*_cos));
		Point2f r2((p2.x*_cos - p2.y*_sin), (p2.x*_sin + p2.y*_cos));
		dj = abs(r2.x - r1.x) + abs(r2.y - r1.y);
	}
	{
		Point2f p1(float(edge_k[0].x), float(edge_k[0].y));
		Point2f p2(float(edge_k[sz_ek - 1].x), float(edge_k[sz_ek - 1].y));
		p1.x -= ell._xc;
		p1.y -= ell._yc;
		p2.x -= ell._xc;
		p2.y -= ell._yc;
		Point2f r1((p1.x*_cos - p1.y*_sin), (p1.x*_sin + p1.y*_cos));
		Point2f r2((p2.x*_cos - p2.y*_sin), (p2.x*_sin + p2.y*_cos));
		dk = abs(r2.x - r1.x) + abs(r2.y - r1.y);
	}

	// This allows to get rid of thick edges
	float rel = min(1.f, ((di + dj + dk) / (3 * (ell._a + ell._b))));
	cout << "rel " << rel << endl;
	if (rel < _fMinReliability)
	{
		Toc(4); //validation
		return;
	}

	// Assign the new score!
	ell._score = (score + rel) * 0.5f;//need to change

	// The tentative detection has been confirmed. Save it!
	ellipses.push_back(ell);

	Toc(4); // Validation
};

bool CNEllipseDetector::FindEllipses(int qi, int qj, int qk, int i, int k, int j, Point2f& center,
	VP& edge_i, VP& edge_j, VP& edge_k,
	EllipseData& data_ij, EllipseData& data_ik,
	vector<Ellipse>& ellipses
)
{//根据三段弧计算椭圆参数，并计算得到的椭圆参数的准确定。
	countsOfFindEllipse++;
	// Find ellipse parameters

	// 0-initialize accumulators
	memset(accN, 0, sizeof(int)*ACC_N_SIZE);
	memset(accR, 0, sizeof(int)*ACC_R_SIZE);
	memset(accA, 0, sizeof(int)*ACC_A_SIZE);

	Tac(3); //estimation

	// Get size of the 4 vectors of slopes (2 pairs of arcs)
	int sz_ij1 = int(data_ij.Sa.size());
	int sz_ij2 = int(data_ij.Sb.size());
	int sz_ik1 = int(data_ik.Sa.size());
	int sz_ik2 = int(data_ik.Sb.size());

	// Get the size of the 3 arcs
	size_t sz_ei = edge_i.size();
	size_t sz_ej = edge_j.size();
	size_t sz_ek = edge_k.size();

	// Center of the estimated ellipse
	float a0 = center.x;
	float b0 = center.y;


	// Estimation of remaining parameters
	// Uses 4 combinations of parameters. See Table 1 and Sect [3.2.3] of the paper.
	//ij1 and ik
	{
		float q1 = data_ij.ra;
		float q3 = data_ik.ra;
		float q5 = data_ik.rb;

		for (int ij1 = 0; ij1 < sz_ij1; ++ij1)
		{
			float q2 = data_ij.Sa[ij1];//有必要遍历吗？

			float q1xq2 = q1 * q2;
			// ij1 and ik1
			for (int ik1 = 0; ik1 < sz_ik1; ++ik1)
			{
				float q4 = data_ik.Sa[ik1];//有必要遍历吗？

				float q3xq4 = q3 * q4;

				// See Eq. [13-18] in the paper

				float a = (q1xq2 - q3xq4);//gama
				float b = (q3xq4 + 1)*(q1 + q2) - (q1xq2 + 1)*(q3 + q4);//beta
				float Kp = (-b + sqrt(b*b + 4 * a*a)) / (2 * a);//K+
				float zplus = ((q1 - Kp)*(q2 - Kp)) / ((1 + q1 * Kp)*(1 + q2 * Kp));
				// 判断有效性  zplus 和K不是线性关系
				if (zplus >= 0.0f) continue;

				float Np = sqrt(-zplus);//N+
				float rho = atan(Kp);//rho tmp
				int rhoDeg;
				if (Np > 1.f)
				{
					Np = 1.f / Np;
					rhoDeg = cvRound((rho * 180 / CV_PI) + 180) % 180; // [0,180)					
				}
				else
				{
					rhoDeg = cvRound((rho * 180 / CV_PI) + 90) % 180; // [0,180)//rho 角度表示，归一化
				}

				int iNp = cvRound(Np * 100); // [0, 100]取整

				if (0 <= iNp && iNp < ACC_N_SIZE &&
					0 <= rhoDeg && rhoDeg < ACC_R_SIZE
					)
				{//为什么都遍历，只遍历一个可以吗  因为 zplus 和K不是线性关系
					++accN[iNp];	// Increment N accumulator
					++accR[rhoDeg];	// Increment R accumulator
				}
			}
			// ij1 and ik2
			for (int ik2 = 0; ik2 < sz_ik2; ++ik2)
			{
				float q4 = data_ik.Sb[ik2];

				float q5xq4 = q5 * q4;

				// See Eq. [13-18] in the paper

				float a = (q1xq2 - q5xq4);
				float b = (q5xq4 + 1)*(q1 + q2) - (q1xq2 + 1)*(q5 + q4);
				float Kp = (-b + sqrt(b*b + 4 * a*a)) / (2 * a);
				float zplus = ((q1 - Kp)*(q2 - Kp)) / ((1 + q1 * Kp)*(1 + q2 * Kp));

				if (zplus >= 0.0f)
				{
					continue;
				}

				float Np = sqrt(-zplus);
				float rho = atan(Kp);
				int rhoDeg;
				if (Np > 1.f)
				{
					Np = 1.f / Np;
					rhoDeg = cvRound((rho * 180 / CV_PI) + 180) % 180; // [0,180)					
				}
				else
				{
					rhoDeg = cvRound((rho * 180 / CV_PI) + 90) % 180; // [0,180)
				}

				int iNp = cvRound(Np * 100); // [0, 100]

				if (0 <= iNp && iNp < ACC_N_SIZE &&
					0 <= rhoDeg && rhoDeg < ACC_R_SIZE
					)
				{
					++accN[iNp];		// Increment N accumulator
					++accR[rhoDeg];		// Increment R accumulator
				}
			}

		}
	}

	//ij2 and ik
	{
		float q1 = data_ij.rb;
		float q3 = data_ik.rb;
		float q5 = data_ik.ra;

		for (int ij2 = 0; ij2 < sz_ij2; ++ij2)
		{
			float q2 = data_ij.Sb[ij2];

			float q1xq2 = q1 * q2;
			//ij2 and ik2
			for (int ik2 = 0; ik2 < sz_ik2; ++ik2)
			{
				float q4 = data_ik.Sb[ik2];

				float q3xq4 = q3 * q4;

				// See Eq. [13-18] in the paper

				float a = (q1xq2 - q3xq4);
				float b = (q3xq4 + 1)*(q1 + q2) - (q1xq2 + 1)*(q3 + q4);
				float Kp = (-b + sqrt(b*b + 4 * a*a)) / (2 * a);
				float zplus = ((q1 - Kp)*(q2 - Kp)) / ((1 + q1 * Kp)*(1 + q2 * Kp));

				if (zplus >= 0.0f)
				{
					continue;
				}

				float Np = sqrt(-zplus);
				float rho = atan(Kp);
				int rhoDeg;
				if (Np > 1.f)
				{
					Np = 1.f / Np;
					rhoDeg = cvRound((rho * 180 / CV_PI) + 180) % 180; // [0,180)
				}
				else
				{
					rhoDeg = cvRound((rho * 180 / CV_PI) + 90) % 180; // [0,180)
				}

				int iNp = cvRound(Np * 100); // [0, 100]

				if (0 <= iNp && iNp < ACC_N_SIZE &&
					0 <= rhoDeg && rhoDeg < ACC_R_SIZE
					)
				{
					++accN[iNp];		// Increment N accumulator
					++accR[rhoDeg];		// Increment R accumulator
				}
			}

			//ij2 and ik1
			for (int ik1 = 0; ik1 < sz_ik1; ++ik1)
			{
				float q4 = data_ik.Sa[ik1];

				float q5xq4 = q5 * q4;

				// See Eq. [13-18] in the paper

				float a = (q1xq2 - q5xq4);
				float b = (q5xq4 + 1)*(q1 + q2) - (q1xq2 + 1)*(q5 + q4);
				float Kp = (-b + sqrt(b*b + 4 * a*a)) / (2 * a);
				float zplus = ((q1 - Kp)*(q2 - Kp)) / ((1 + q1 * Kp)*(1 + q2 * Kp));

				if (zplus >= 0.0f)
				{
					continue;
				}

				float Np = sqrt(-zplus);
				float rho = atan(Kp);
				int rhoDeg;
				if (Np > 1.f)
				{
					Np = 1.f / Np;
					rhoDeg = cvRound((rho * 180 / CV_PI) + 180) % 180; // [0,180)
				}
				else
				{
					rhoDeg = cvRound((rho * 180 / CV_PI) + 90) % 180; // [0,180)
				}

				int iNp = cvRound(Np * 100); // [0, 100]

				if (0 <= iNp && iNp < ACC_N_SIZE &&
					0 <= rhoDeg && rhoDeg < ACC_R_SIZE
					)
				{
					++accN[iNp];		// Increment N accumulator
					++accR[rhoDeg];		// Increment R accumulator
				}
			}

		}
	}

	// Find peak in N and K accumulator
	int iN = FindMaxN(accN);
	int iK = FindMaxK(accR);

	// Recover real values
	float fK = float(iK);
	float Np = float(iN) * 0.01f;
	float rho = fK * float(CV_PI) / 180.f;	//deg 2 rad
	float Kp = tan(rho);

	// Estimate A. See Eq. [19 - 22] in Sect [3.2.3] of the paper  
	// 三段弧上的点都用来估算A
	//可以优化   
	for (ushort l = 0; l < sz_ei; ++l)
	{
		Point& pp = edge_i[l];
		float sk = 1.f / sqrt(Kp*Kp + 1.f);//cos rho
		float x0 = ((pp.x - a0) * sk) + (((pp.y - b0)*Kp) * sk);//可优化
		float y0 = -(((pp.x - a0) * Kp) * sk) + ((pp.y - b0) * sk);//可优化
		float Ax = sqrt((x0*x0*Np*Np + y0 * y0) / ((Np*Np)*(1.f + Kp * Kp)));
		int A = cvRound(abs(Ax / cos(rho)));//可以优化
		if ((0 <= A) && (A < ACC_A_SIZE))
		{
			++accA[A];
		}
	}

	for (ushort l = 0; l < sz_ej; ++l)
	{
		Point& pp = edge_j[l];
		float sk = 1.f / sqrt(Kp*Kp + 1.f);
		float x0 = ((pp.x - a0) * sk) + (((pp.y - b0)*Kp) * sk);
		float y0 = -(((pp.x - a0) * Kp) * sk) + ((pp.y - b0) * sk);
		float Ax = sqrt((x0*x0*Np*Np + y0 * y0) / ((Np*Np)*(1.f + Kp * Kp)));
		int A = cvRound(abs(Ax / cos(rho)));
		if ((0 <= A) && (A < ACC_A_SIZE))
		{
			++accA[A];
		}
	}

	for (ushort l = 0; l < sz_ek; ++l)
	{
		Point& pp = edge_k[l];
		float sk = 1.f / sqrt(Kp*Kp + 1.f);
		float x0 = ((pp.x - a0) * sk) + (((pp.y - b0)*Kp) * sk);
		float y0 = -(((pp.x - a0) * Kp) * sk) + ((pp.y - b0) * sk);
		float Ax = sqrt((x0*x0*Np*Np + y0 * y0) / ((Np*Np)*(1.f + Kp * Kp)));
		int A = cvRound(abs(Ax / cos(rho)));
		if ((0 <= A) && (A < ACC_A_SIZE))
		{
			++accA[A];
		}
	}

	// Find peak in A accumulator
	int A = FindMaxA(accA);
	float fA = float(A);

	// Find B value. See Eq [23] in the paper
	float fB = abs(fA * Np);

	// Got all ellipse parameters!
	Ellipse ell(a0, b0, fA, fB, fmod(rho + float(CV_PI)*2.f, float(CV_PI)));

	Toc(3); //estimation
	Tac(4); //validation

	// Get the score. See Sect [3.3.1] in the paper

	// Find the number of edge pixel lying on the ellipse
	float _cos = cos(-ell._rad);
	float _sin = sin(-ell._rad);

	float invA2 = 1.f / (ell._a * ell._a);
	float invB2 = 1.f / (ell._b * ell._b);

	float invNofPoints = 1.f / float(sz_ei + sz_ej + sz_ek);
	int counter_on_perimeter = 0;

	for (ushort l = 0; l < sz_ei; ++l)
	{
		float tx = float(edge_i[l].x) - ell._xc;
		float ty = float(edge_i[l].y) - ell._yc;
		float rx = (tx*_cos - ty * _sin);
		float ry = (tx*_sin + ty * _cos);

		float h = (rx*rx)*invA2 + (ry*ry)*invB2;
		if (abs(h - 1.f) < _fDistanceToEllipseContour)
		{
			++counter_on_perimeter;
		}
	}

	for (ushort l = 0; l < sz_ej; ++l)
	{
		float tx = float(edge_j[l].x) - ell._xc;
		float ty = float(edge_j[l].y) - ell._yc;
		float rx = (tx*_cos - ty * _sin);
		float ry = (tx*_sin + ty * _cos);

		float h = (rx*rx)*invA2 + (ry*ry)*invB2;
		if (abs(h - 1.f) < _fDistanceToEllipseContour)
		{
			++counter_on_perimeter;
		}
	}

	for (ushort l = 0; l < sz_ek; ++l)
	{
		float tx = float(edge_k[l].x) - ell._xc;
		float ty = float(edge_k[l].y) - ell._yc;
		float rx = (tx*_cos - ty * _sin);
		float ry = (tx*_sin + ty * _cos);

		float h = (rx*rx)*invA2 + (ry*ry)*invB2;
		if (abs(h - 1.f) < _fDistanceToEllipseContour)
		{
			++counter_on_perimeter;
		}
	}

	//no points found on the ellipse
	if (counter_on_perimeter <= 0)
	{
		Toc(4); //validation
		return false;
	}


	// Compute score
	float score = float(counter_on_perimeter) * invNofPoints;
	//cout << "in validation score is "<< score << endl;
	if (score < _fMinScore)
	{
		Toc(4); //validation
		return false;
	}

	// Compute reliability	
	// this metric is not described in the paper, mostly due to space limitations.
	// The main idea is that for a given ellipse (TD) even if the score is high, the arcs 
	// can cover only a small amount of the contour of the estimated ellipse. 
	// A low reliability indicate that the arcs form an elliptic shape by chance, but do not underlie
	// an actual ellipse. The value is normalized between 0 and 1. 
	// The default value is 0.4.

	// It is somehow similar to the "Angular Circumreference Ratio" saliency criteria 
	// as in the paper: 
	// D. K. Prasad, M. K. Leung, S.-Y. Cho, Edge curvature and convexity
	// based ellipse detection method, Pattern Recognition 45 (2012) 3204-3221.

	float di, dj, dk;
	{
		Point2f p1(float(edge_i[0].x), float(edge_i[0].y));
		Point2f p2(float(edge_i[sz_ei - 1].x), float(edge_i[sz_ei - 1].y));
		p1.x -= ell._xc;
		p1.y -= ell._yc;
		p2.x -= ell._xc;
		p2.y -= ell._yc;
		Point2f r1((p1.x*_cos - p1.y*_sin), (p1.x*_sin + p1.y*_cos));
		Point2f r2((p2.x*_cos - p2.y*_sin), (p2.x*_sin + p2.y*_cos));
		di = abs(r2.x - r1.x) + abs(r2.y - r1.y);
	}
	{
		Point2f p1(float(edge_j[0].x), float(edge_j[0].y));
		Point2f p2(float(edge_j[sz_ej - 1].x), float(edge_j[sz_ej - 1].y));
		p1.x -= ell._xc;
		p1.y -= ell._yc;
		p2.x -= ell._xc;
		p2.y -= ell._yc;
		Point2f r1((p1.x*_cos - p1.y*_sin), (p1.x*_sin + p1.y*_cos));
		Point2f r2((p2.x*_cos - p2.y*_sin), (p2.x*_sin + p2.y*_cos));
		dj = abs(r2.x - r1.x) + abs(r2.y - r1.y);
	}
	{
		Point2f p1(float(edge_k[0].x), float(edge_k[0].y));
		Point2f p2(float(edge_k[sz_ek - 1].x), float(edge_k[sz_ek - 1].y));
		p1.x -= ell._xc;
		p1.y -= ell._yc;
		p2.x -= ell._xc;
		p2.y -= ell._yc;
		Point2f r1((p1.x*_cos - p1.y*_sin), (p1.x*_sin + p1.y*_cos));
		Point2f r2((p2.x*_cos - p2.y*_sin), (p2.x*_sin + p2.y*_cos));
		dk = abs(r2.x - r1.x) + abs(r2.y - r1.y);
	}

	// This allows to get rid of thick edges
	float rel = min(1.f, ((di + dj + dk) / (3 * (ell._a + ell._b))));
	//cout << "in validation rel is" << rel << endl;
	if (rel < _fMinReliability)
	{
		Toc(4); //validation
		return false;
	}

	// Assign the new score!
	ell._score = (score + rel) * 0.5f;//need to change

	// The tentative detection has been confirmed. Save it!
	//zrk add rect into the ell
	ellipses.push_back(ell);

	ArcsInfo a(qi, i);
	ArcsInfo b(qj, j);
	ArcsInfo c(qk, k);
	ell._arcs[0] = a;
	ell._arcs[1] = b;
	ell._arcs[2] = c;

	Toc(4); // Validation
	return true;
};

bool CNEllipseDetector::FindEllipses(Point2f& center,
	VP& edge_i, VP& edge_j, VP& edge_k,
	EllipseData& data_ij, EllipseData& data_ik,
	vector<Ellipse>& ellipses,
	Rectangle rect)
{//根据三段弧计算椭圆参数，并计算得到的椭圆参数的准确定。

	countsOfFindEllipse++;
	// Find ellipse parameters

	// 0-initialize accumulators
	memset(accN, 0, sizeof(int)*ACC_N_SIZE);
	memset(accR, 0, sizeof(int)*ACC_R_SIZE);
	memset(accA, 0, sizeof(int)*ACC_A_SIZE);

	Tac(3); //estimation

	// Get size of the 4 vectors of slopes (2 pairs of arcs)
	int sz_ij1 = int(data_ij.Sa.size());
	int sz_ij2 = int(data_ij.Sb.size());
	int sz_ik1 = int(data_ik.Sa.size());
	int sz_ik2 = int(data_ik.Sb.size());

	// Get the size of the 3 arcs
	size_t sz_ei = edge_i.size();
	size_t sz_ej = edge_j.size();
	size_t sz_ek = edge_k.size();

	// Center of the estimated ellipse
	float a0 = center.x;
	float b0 = center.y;


	// Estimation of remaining parameters
	// Uses 4 combinations of parameters. See Table 1 and Sect [3.2.3] of the paper.
	//ij1 and ik
	{
		float q1 = data_ij.ra;
		float q3 = data_ik.ra;
		float q5 = data_ik.rb;

		for (int ij1 = 0; ij1 < sz_ij1; ++ij1)
		{
			float q2 = data_ij.Sa[ij1];//有必要遍历吗？

			float q1xq2 = q1 * q2;
			// ij1 and ik1
			for (int ik1 = 0; ik1 < sz_ik1; ++ik1)
			{
				float q4 = data_ik.Sa[ik1];//有必要遍历吗？

				float q3xq4 = q3 * q4;

				// See Eq. [13-18] in the paper

				float a = (q1xq2 - q3xq4);//gama
				float b = (q3xq4 + 1)*(q1 + q2) - (q1xq2 + 1)*(q3 + q4);//beta
				float Kp = (-b + sqrt(b*b + 4 * a*a)) / (2 * a);//K+
				float zplus = ((q1 - Kp)*(q2 - Kp)) / ((1 + q1 * Kp)*(1 + q2 * Kp));
				// 判断有效性  zplus 和K不是线性关系
				if (zplus >= 0.0f) continue;

				float Np = sqrt(-zplus);//N+
				float rho = atan(Kp);//rho tmp
				int rhoDeg;
				if (Np > 1.f)
				{
					Np = 1.f / Np;
					rhoDeg = cvRound((rho * 180 / CV_PI) + 180) % 180; // [0,180)					
				}
				else
				{
					rhoDeg = cvRound((rho * 180 / CV_PI) + 90) % 180; // [0,180)//rho 角度表示，归一化
				}

				int iNp = cvRound(Np * 100); // [0, 100]取整

				if (0 <= iNp && iNp < ACC_N_SIZE &&
					0 <= rhoDeg && rhoDeg < ACC_R_SIZE
					)
				{//为什么都遍历，只遍历一个可以吗  因为 zplus 和K不是线性关系
					++accN[iNp];	// Increment N accumulator
					++accR[rhoDeg];	// Increment R accumulator
				}
			}
			// ij1 and ik2
			for (int ik2 = 0; ik2 < sz_ik2; ++ik2)
			{
				float q4 = data_ik.Sb[ik2];

				float q5xq4 = q5 * q4;

				// See Eq. [13-18] in the paper

				float a = (q1xq2 - q5xq4);
				float b = (q5xq4 + 1)*(q1 + q2) - (q1xq2 + 1)*(q5 + q4);
				float Kp = (-b + sqrt(b*b + 4 * a*a)) / (2 * a);
				float zplus = ((q1 - Kp)*(q2 - Kp)) / ((1 + q1 * Kp)*(1 + q2 * Kp));

				if (zplus >= 0.0f)
				{
					continue;
				}

				float Np = sqrt(-zplus);
				float rho = atan(Kp);
				int rhoDeg;
				if (Np > 1.f)
				{
					Np = 1.f / Np;
					rhoDeg = cvRound((rho * 180 / CV_PI) + 180) % 180; // [0,180)					
				}
				else
				{
					rhoDeg = cvRound((rho * 180 / CV_PI) + 90) % 180; // [0,180)
				}

				int iNp = cvRound(Np * 100); // [0, 100]

				if (0 <= iNp && iNp < ACC_N_SIZE &&
					0 <= rhoDeg && rhoDeg < ACC_R_SIZE
					)
				{
					++accN[iNp];		// Increment N accumulator
					++accR[rhoDeg];		// Increment R accumulator
				}
			}

		}
	}

	//ij2 and ik
	{
		float q1 = data_ij.rb;
		float q3 = data_ik.rb;
		float q5 = data_ik.ra;

		for (int ij2 = 0; ij2 < sz_ij2; ++ij2)
		{
			float q2 = data_ij.Sb[ij2];

			float q1xq2 = q1 * q2;
			//ij2 and ik2
			for (int ik2 = 0; ik2 < sz_ik2; ++ik2)
			{
				float q4 = data_ik.Sb[ik2];

				float q3xq4 = q3 * q4;

				// See Eq. [13-18] in the paper

				float a = (q1xq2 - q3xq4);
				float b = (q3xq4 + 1)*(q1 + q2) - (q1xq2 + 1)*(q3 + q4);
				float Kp = (-b + sqrt(b*b + 4 * a*a)) / (2 * a);
				float zplus = ((q1 - Kp)*(q2 - Kp)) / ((1 + q1 * Kp)*(1 + q2 * Kp));

				if (zplus >= 0.0f)
				{
					continue;
				}

				float Np = sqrt(-zplus);
				float rho = atan(Kp);
				int rhoDeg;
				if (Np > 1.f)
				{
					Np = 1.f / Np;
					rhoDeg = cvRound((rho * 180 / CV_PI) + 180) % 180; // [0,180)
				}
				else
				{
					rhoDeg = cvRound((rho * 180 / CV_PI) + 90) % 180; // [0,180)
				}

				int iNp = cvRound(Np * 100); // [0, 100]

				if (0 <= iNp && iNp < ACC_N_SIZE &&
					0 <= rhoDeg && rhoDeg < ACC_R_SIZE
					)
				{
					++accN[iNp];		// Increment N accumulator
					++accR[rhoDeg];		// Increment R accumulator
				}
			}

			//ij2 and ik1
			for (int ik1 = 0; ik1 < sz_ik1; ++ik1)
			{
				float q4 = data_ik.Sa[ik1];

				float q5xq4 = q5 * q4;

				// See Eq. [13-18] in the paper

				float a = (q1xq2 - q5xq4);
				float b = (q5xq4 + 1)*(q1 + q2) - (q1xq2 + 1)*(q5 + q4);
				float Kp = (-b + sqrt(b*b + 4 * a*a)) / (2 * a);
				float zplus = ((q1 - Kp)*(q2 - Kp)) / ((1 + q1 * Kp)*(1 + q2 * Kp));

				if (zplus >= 0.0f)
				{
					continue;
				}

				float Np = sqrt(-zplus);
				float rho = atan(Kp);
				int rhoDeg;
				if (Np > 1.f)
				{
					Np = 1.f / Np;
					rhoDeg = cvRound((rho * 180 / CV_PI) + 180) % 180; // [0,180)
				}
				else
				{
					rhoDeg = cvRound((rho * 180 / CV_PI) + 90) % 180; // [0,180)
				}

				int iNp = cvRound(Np * 100); // [0, 100]

				if (0 <= iNp && iNp < ACC_N_SIZE &&
					0 <= rhoDeg && rhoDeg < ACC_R_SIZE
					)
				{
					++accN[iNp];		// Increment N accumulator
					++accR[rhoDeg];		// Increment R accumulator
				}
			}

		}
	}

	// Find peak in N and K accumulator
	int iN = FindMaxN(accN);
	int iK = FindMaxK(accR);

	// Recover real values
	float fK = float(iK);
	float Np = float(iN) * 0.01f;
	float rho = fK * float(CV_PI) / 180.f;	//deg 2 rad
	float Kp = tan(rho);

	// Estimate A. See Eq. [19 - 22] in Sect [3.2.3] of the paper  
	// 三段弧上的点都用来估算A
	//可以优化   
	for (ushort l = 0; l < sz_ei; ++l)
	{
		Point& pp = edge_i[l];
		float sk = 1.f / sqrt(Kp*Kp + 1.f);//cos rho
		float x0 = ((pp.x - a0) * sk) + (((pp.y - b0)*Kp) * sk);//可优化
		float y0 = -(((pp.x - a0) * Kp) * sk) + ((pp.y - b0) * sk);//可优化
		float Ax = sqrt((x0*x0*Np*Np + y0 * y0) / ((Np*Np)*(1.f + Kp * Kp)));
		int A = cvRound(abs(Ax / cos(rho)));//可以优化
		if ((0 <= A) && (A < ACC_A_SIZE))
		{
			++accA[A];
		}
	}

	for (ushort l = 0; l < sz_ej; ++l)
	{
		Point& pp = edge_j[l];
		float sk = 1.f / sqrt(Kp*Kp + 1.f);
		float x0 = ((pp.x - a0) * sk) + (((pp.y - b0)*Kp) * sk);
		float y0 = -(((pp.x - a0) * Kp) * sk) + ((pp.y - b0) * sk);
		float Ax = sqrt((x0*x0*Np*Np + y0 * y0) / ((Np*Np)*(1.f + Kp * Kp)));
		int A = cvRound(abs(Ax / cos(rho)));
		if ((0 <= A) && (A < ACC_A_SIZE))
		{
			++accA[A];
		}
	}

	for (ushort l = 0; l < sz_ek; ++l)
	{
		Point& pp = edge_k[l];
		float sk = 1.f / sqrt(Kp*Kp + 1.f);
		float x0 = ((pp.x - a0) * sk) + (((pp.y - b0)*Kp) * sk);
		float y0 = -(((pp.x - a0) * Kp) * sk) + ((pp.y - b0) * sk);
		float Ax = sqrt((x0*x0*Np*Np + y0 * y0) / ((Np*Np)*(1.f + Kp * Kp)));
		int A = cvRound(abs(Ax / cos(rho)));
		if ((0 <= A) && (A < ACC_A_SIZE))
		{
			++accA[A];
		}
	}

	// Find peak in A accumulator
	int A = FindMaxA(accA);
	float fA = float(A);

	// Find B value. See Eq [23] in the paper
	float fB = abs(fA * Np);

	// Got all ellipse parameters!
	Ellipse ell(a0, b0, fA, fB, fmod(rho + float(CV_PI)*2.f, float(CV_PI)));

	Toc(3); //estimation
	Tac(4); //validation

	// Get the score. See Sect [3.3.1] in the paper

	// Find the number of edge pixel lying on the ellipse
	float _cos = cos(-ell._rad);
	float _sin = sin(-ell._rad);

	float invA2 = 1.f / (ell._a * ell._a);
	float invB2 = 1.f / (ell._b * ell._b);

	float invNofPoints = 1.f / float(sz_ei + sz_ej + sz_ek);
	int counter_on_perimeter = 0;

	for (ushort l = 0; l < sz_ei; ++l)
	{
		float tx = float(edge_i[l].x) - ell._xc;
		float ty = float(edge_i[l].y) - ell._yc;
		float rx = (tx*_cos - ty * _sin);
		float ry = (tx*_sin + ty * _cos);

		float h = (rx*rx)*invA2 + (ry*ry)*invB2;
		if (abs(h - 1.f) < _fDistanceToEllipseContour)
		{
			++counter_on_perimeter;
		}
	}

	for (ushort l = 0; l < sz_ej; ++l)
	{
		float tx = float(edge_j[l].x) - ell._xc;
		float ty = float(edge_j[l].y) - ell._yc;
		float rx = (tx*_cos - ty * _sin);
		float ry = (tx*_sin + ty * _cos);

		float h = (rx*rx)*invA2 + (ry*ry)*invB2;
		if (abs(h - 1.f) < _fDistanceToEllipseContour)
		{
			++counter_on_perimeter;
		}
	}

	for (ushort l = 0; l < sz_ek; ++l)
	{
		float tx = float(edge_k[l].x) - ell._xc;
		float ty = float(edge_k[l].y) - ell._yc;
		float rx = (tx*_cos - ty * _sin);
		float ry = (tx*_sin + ty * _cos);

		float h = (rx*rx)*invA2 + (ry*ry)*invB2;
		if (abs(h - 1.f) < _fDistanceToEllipseContour)
		{
			++counter_on_perimeter;
		}
	}

	//no points found on the ellipse
	if (counter_on_perimeter <= 0)
	{
		Toc(4); //validation
		return false;
	}


	// Compute score
	float score = float(counter_on_perimeter) * invNofPoints;
	cout << "in validation score is " << score << endl;
	if (score < _fMinScore)
	{
		Toc(4); //validation
		return false;
	}

	// Compute reliability	
	// this metric is not described in the paper, mostly due to space limitations.
	// The main idea is that for a given ellipse (TD) even if the score is high, the arcs 
	// can cover only a small amount of the contour of the estimated ellipse. 
	// A low reliability indicate that the arcs form an elliptic shape by chance, but do not underlie
	// an actual ellipse. The value is normalized between 0 and 1. 
	// The default value is 0.4.

	// It is somehow similar to the "Angular Circumreference Ratio" saliency criteria 
	// as in the paper: 
	// D. K. Prasad, M. K. Leung, S.-Y. Cho, Edge curvature and convexity
	// based ellipse detection method, Pattern Recognition 45 (2012) 3204-3221.

	float di, dj, dk;
	{
		Point2f p1(float(edge_i[0].x), float(edge_i[0].y));
		Point2f p2(float(edge_i[sz_ei - 1].x), float(edge_i[sz_ei - 1].y));
		p1.x -= ell._xc;
		p1.y -= ell._yc;
		p2.x -= ell._xc;
		p2.y -= ell._yc;
		Point2f r1((p1.x*_cos - p1.y*_sin), (p1.x*_sin + p1.y*_cos));
		Point2f r2((p2.x*_cos - p2.y*_sin), (p2.x*_sin + p2.y*_cos));
		di = abs(r2.x - r1.x) + abs(r2.y - r1.y);
	}
	{
		Point2f p1(float(edge_j[0].x), float(edge_j[0].y));
		Point2f p2(float(edge_j[sz_ej - 1].x), float(edge_j[sz_ej - 1].y));
		p1.x -= ell._xc;
		p1.y -= ell._yc;
		p2.x -= ell._xc;
		p2.y -= ell._yc;
		Point2f r1((p1.x*_cos - p1.y*_sin), (p1.x*_sin + p1.y*_cos));
		Point2f r2((p2.x*_cos - p2.y*_sin), (p2.x*_sin + p2.y*_cos));
		dj = abs(r2.x - r1.x) + abs(r2.y - r1.y);
	}
	{
		Point2f p1(float(edge_k[0].x), float(edge_k[0].y));
		Point2f p2(float(edge_k[sz_ek - 1].x), float(edge_k[sz_ek - 1].y));
		p1.x -= ell._xc;
		p1.y -= ell._yc;
		p2.x -= ell._xc;
		p2.y -= ell._yc;
		Point2f r1((p1.x*_cos - p1.y*_sin), (p1.x*_sin + p1.y*_cos));
		Point2f r2((p2.x*_cos - p2.y*_sin), (p2.x*_sin + p2.y*_cos));
		dk = abs(r2.x - r1.x) + abs(r2.y - r1.y);
	}

	// This allows to get rid of thick edges
	float rel = min(1.f, ((di + dj + dk) / (3 * (ell._a + ell._b))));
	cout << "in validation rel is " << rel << endl;
	if (rel < _fMinReliability)
	{
		Toc(4); //validation
		return false;
	}
	//cout << "_fMinReliability" << _fMinReliability << endl;
	// Assign the new score!
	ell._score = (score + rel) * 0.5f;//need to change
	//cout << "assign the new score is " << ell._score << endl;
	// The tentative detection has been confirmed. Save it!
	//zrk add rect into the ell
	ell._rect = rect;
	//构成椭圆，输出一个输出一个

	ellipses.push_back(ell);

	Toc(4); // Validation
	return true;
};

// Get the coordinates of the center, given the intersection of the estimated lines. See Fig. [8] in Sect [3.2.3] in the paper.
Point2f CNEllipseDetector::GetCenterCoordinates(EllipseData& data_ij, EllipseData& data_ik)
{
	float xx[7];
	float yy[7];

	xx[0] = data_ij.Cab.x;
	xx[1] = data_ik.Cab.x;
	yy[0] = data_ij.Cab.y;
	yy[1] = data_ik.Cab.y;

	{
		//1-1
		float q2 = data_ij.ta;
		float q4 = data_ik.ta;
		Point2f& M12 = data_ij.Ma;
		Point2f& M34 = data_ik.Ma;

		float invDen = 1 / (q2 - q4);
		xx[2] = (M34.y - q4*M34.x - M12.y + q2*M12.x) * invDen;
		yy[2] = (q2*M34.y - q4*M12.y + q2*q4*(M12.x - M34.x)) * invDen;
	}

	{
		//1-2
		float q2 = data_ij.ta;
		float q4 = data_ik.tb;
		Point2f& M12 = data_ij.Ma;
		Point2f& M34 = data_ik.Mb;

		float invDen = 1 / (q2 - q4);
		xx[3] = (M34.y - q4*M34.x - M12.y + q2*M12.x) * invDen;
		yy[3] = (q2*M34.y - q4*M12.y + q2*q4*(M12.x - M34.x)) * invDen;
	}

	{
		//2-2
		float q2 = data_ij.tb;
		float q4 = data_ik.tb;
		Point2f& M12 = data_ij.Mb;
		Point2f& M34 = data_ik.Mb;

		float invDen = 1 / (q2 - q4);
		xx[4] = (M34.y - q4*M34.x - M12.y + q2*M12.x) * invDen;
		yy[4] = (q2*M34.y - q4*M12.y + q2*q4*(M12.x - M34.x)) * invDen;
	}

	{
		//2-1
		float q2 = data_ij.tb;
		float q4 = data_ik.ta;
		Point2f& M12 = data_ij.Mb;
		Point2f& M34 = data_ik.Ma;

		float invDen = 1 / (q2 - q4);
		xx[5] = (M34.y - q4*M34.x - M12.y + q2*M12.x) * invDen;
		yy[5] = (q2*M34.y - q4*M12.y + q2*q4*(M12.x - M34.x)) * invDen;
	}

	xx[6] = (xx[0] + xx[1]) * 0.5f;
	yy[6] = (yy[0] + yy[1]) * 0.5f;


	// Median
	nth_element(xx, xx + 3, xx + 7);
	nth_element(yy, yy + 3, yy + 7);
	float xc = xx[3];
	float yc = yy[3];

	return Point2f(xc, yc);
};

//123456 124 80 48 246
//#define T124 pil,pim,pif,pjl,pjm,pjf

#define T124 pjf,pjm,pjl,pif,pim,pil //原来的
#define T231 pil,pim,pif,pjf,pjm,pjl
#define T342 pif,pim,pil,pjf,pjm,pjl
#define T413 pif,pim,pil,pjl,pjm,pjf



// Verify triplets of arcs with convexity: i=1, j=2, k=4

void CNEllipseDetector::Triplets124(VVP& pi,
									   VVP& pj,
									   VVP& pk,
									   unordered_map<uint, EllipseData>& data,
									   vector<Ellipse>& ellipses
									   )
{
	// get arcs length
	ushort sz_i = ushort(pi.size());
	ushort sz_j = ushort(pj.size());
	ushort sz_k = ushort(pk.size());

	// For each edge i
	for (ushort i = 0; i < sz_i; ++i)
	{
		VP& edge_i = pi[i];
		ushort sz_ei = ushort(edge_i.size());

		Point& pif = edge_i[0];
		Point& pim = edge_i[sz_ei/2];
		Point& pil = edge_i[sz_ei - 1];

		// 1,2 -> reverse 1, swap
		VP rev_i(edge_i.size());
		reverse_copy(edge_i.begin(), edge_i.end(), rev_i.begin());

		// For each edge j
		for (ushort j = 0; j < sz_j; ++j)
		{
			VP& edge_j = pj[j];
			ushort sz_ej = ushort(edge_j.size());

			Point& pjf = edge_j[0];
			Point& pjm = edge_j[sz_ej/2];
			Point& pjl = edge_j[sz_ej - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
			// CONSTRAINTS on position
			if (pjl.x > pif.x + _fThPosition) //is right	
			{
				//discard
				continue;
			}
#endif
#ifndef DISCARD_CONSTRAINT_CNC
			// 添加限定条件  2se se1//pil,pim,pif,pjf,pjm,pjl pjf,pjm,pjl,pif,pim,pil
			if(myselect1&&fabs(value4SixPoints(T124)-1)>tCNC)
			{
				continue;
			}
#endif
			uint key_ij = GenerateKey(PAIR_12, i, j);

			//for each edge k
			for (ushort k = 0; k < sz_k; ++k)
			{
				VP& edge_k = pk[k];
				ushort sz_ek = ushort(edge_k.size());

				Point& pkf = edge_k[0];
				Point& pkm = edge_k[sz_ek/2];
				Point& pkl = edge_k[sz_ek - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
				//CONSTRAINTS on position
				if (pkl.y < pil.y - _fThPosition)
				{
					//discard
					continue;
				}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
				// 添加限定条件2
				if(myselect2&&fabs(value4SixPoints(pif,pim,pil,pkf,pkm,pkl)-1)>tCNC)
				{
					continue;
				}
				// 添加限定条件3
				if(myselect3&&fabs(value4SixPoints(pjf,pjm,pjl,pkf,pkm,pkl)-1)>tCNC)
				{
					continue;
				}
#endif

				uint key_ik = GenerateKey(PAIR_14, i, k);

				// Find centers

				EllipseData data_ij, data_ik;

				// If the data for the pair i-j have not been computed yet
				if (data.count(key_ij) == 0)
				{
					//1,2 -> reverse 1, swap

					// Compute data!
					GetFastCenter(edge_j, rev_i, data_ij);
					// Insert computed data in the hash table
					data.insert(pair<uint, EllipseData>(key_ij, data_ij));
				}
				else
				{
					// Otherwise, just lookup the data in the hash table
					data_ij = data.at(key_ij);
				}

				// If the data for the pair i-k have not been computed yet
				if (data.count(key_ik) == 0)
				{
					//1,4 -> ok

					// Compute data!
					GetFastCenter(edge_i, edge_k, data_ik);
					// Insert computed data in the hash table
					data.insert(pair<uint, EllipseData>(key_ik, data_ik));
				}
				else
				{
					// Otherwise, just lookup the data in the hash table
					data_ik = data.at(key_ik);
				}

				// INVALID CENTERS
				if (!data_ij.isValid || !data_ik.isValid)
				{
					continue;
				}

#ifndef DISCARD_CONSTRAINT_CENTER
				// Selection strategy - Step 3. See Sect [3.2.2] in the paper
				// The computed centers are not close enough
				if (ed2(data_ij.Cab, data_ik.Cab) > _fMaxCenterDistance2)
				{
					//discard
					continue;
				}
#endif
				// If all constraints of the selection strategy have been satisfied, 
				// we can start estimating the ellipse parameters

				// Find ellipse parameters

				// Get the coordinates of the center (xc, yc)
				Point2f center = GetCenterCoordinates(data_ij, data_ik);

				// Find remaining paramters (A,B,rho)
				FindEllipses(center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses);
			}
		}
	}
};



void CNEllipseDetector::Triplets231(VVP& pi,
									   VVP& pj,
									   VVP& pk,
									   unordered_map<uint, EllipseData>& data,
									   vector<Ellipse>& ellipses
									   )
{
	ushort sz_i = ushort(pi.size());
	ushort sz_j = ushort(pj.size());
	ushort sz_k = ushort(pk.size());

	// For each edge i
	for (ushort i = 0; i < sz_i; ++i)
	{
		VP& edge_i = pi[i];
		ushort sz_ei = ushort(edge_i.size());

		Point& pif = edge_i[0];
		Point& pim = edge_i[sz_ei/2];
		Point& pil = edge_i[sz_ei - 1];

		VP rev_i(edge_i.size());
		reverse_copy(edge_i.begin(), edge_i.end(), rev_i.begin());

		// For each edge j
		for (ushort j = 0; j < sz_j; ++j)
		{
			VP& edge_j = pj[j];
			ushort sz_ej = ushort(edge_j.size());

			Point& pjf = edge_j[0];
			Point& pjm = edge_j[sz_ej/2];
			Point& pjl = edge_j[sz_ej - 1];			

#ifndef DISCARD_CONSTRAINT_POSITION
			// CONSTRAINTS on position
			if (pjf.y < pif.y - _fThPosition)
			{
				//discard
				continue;
			}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
			//添加约束
			// 添加限定条件 2es se3 //pif,pim,pil,pjf,pjm,pjl pil,pim,pif,pjf,pjm,pjl
			if(myselect1&&fabs(value4SixPoints(T231)-1)>tCNC)
			{
				continue;
			}
#endif

			VP rev_j(edge_j.size());
			reverse_copy(edge_j.begin(), edge_j.end(), rev_j.begin());

			uint key_ij = GenerateKey(PAIR_23, i, j);

			// For each edge k
			for (ushort k = 0; k < sz_k; ++k)
			{
				VP& edge_k = pk[k];
				ushort sz_ek = ushort(edge_k.size());

				Point& pkf = edge_k[0];
				Point& pkm = edge_k[sz_ek/2];
				Point& pkl = edge_k[sz_ek - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
				// CONSTRAINTS on position
				if (pkf.x < pil.x - _fThPosition)
				{
					//discard
					continue;
				}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
				// 添加限定条件2
				if(myselect2&&fabs(value4SixPoints(pif,pim,pil,pkf,pkm,pkl)-1)>tCNC)
				{
					continue;
				}
				// 添加限定条件3
				if(myselect3&&fabs(value4SixPoints(pjf,pjm,pjl,pkf,pkm,pkl)-1)>tCNC)
				{
					continue;
				}
#endif
				uint key_ik = GenerateKey(PAIR_12, k, i);

				// Find centers

				EllipseData data_ij, data_ik;

				if (data.count(key_ij) == 0)
				{
					// 2,3 -> reverse 2,3

					GetFastCenter(rev_i, rev_j, data_ij);
					data.insert(pair<uint, EllipseData>(key_ij, data_ij));
				}
				else
				{
					data_ij = data.at(key_ij);
				}

				if (data.count(key_ik) == 0)
				{
					// 2,1 -> reverse 1
					VP rev_k(edge_k.size());
					reverse_copy(edge_k.begin(), edge_k.end(), rev_k.begin());

					GetFastCenter(edge_i, rev_k, data_ik);
					data.insert(pair<uint, EllipseData>(key_ik, data_ik));
				}
				else
				{
					data_ik = data.at(key_ik);
				}

				// INVALID CENTERS
				if (!data_ij.isValid || !data_ik.isValid)
				{
					continue;
				}

#ifndef DISCARD_CONSTRAINT_CENTER
				// CONSTRAINT ON CENTERS
				if (ed2(data_ij.Cab, data_ik.Cab) > _fMaxCenterDistance2)
				{
					//discard
					continue;
				}
#endif
				// Find ellipse parameters
				Point2f center = GetCenterCoordinates(data_ij, data_ik);

				FindEllipses(center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses);

			}
		}
	}
};


void CNEllipseDetector::Triplets342(VVP& pi,
									   VVP& pj,
									   VVP& pk,
									   unordered_map<uint, EllipseData>& data,
									   vector<Ellipse>& ellipses
									   )
{
	ushort sz_i = ushort(pi.size());
	ushort sz_j = ushort(pj.size());
	ushort sz_k = ushort(pk.size());

	// For each edge i
	for (ushort i = 0; i < sz_i; ++i)
	{
		VP& edge_i = pi[i];
		ushort sz_ei = ushort(edge_i.size());

		Point& pif = edge_i[0];
		Point& pim = edge_i[sz_ei/2];
		Point& pil = edge_i[sz_ei - 1];

		VP rev_i(edge_i.size());
		reverse_copy(edge_i.begin(), edge_i.end(), rev_i.begin());

		// For each edge j
		for (ushort j = 0; j < sz_j; ++j)
		{
			VP& edge_j = pj[j];
			ushort sz_ej = ushort(edge_j.size());

			Point& pjf = edge_j[0];
			Point& pjm = edge_j[sz_ej/2];
			Point& pjl = edge_j[sz_ej - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
			//CONSTRAINTS on position
			if (pjf.x < pil.x - _fThPosition) 		//is left
			{
				//discard
				continue;
			}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
			// 添加限定条件 3se se4 // pil,pim,pif,pjf,pjm,pjl pif,pim,pil,pjf,pjm,pjl
			if(myselect1&&fabs(value4SixPoints(T342)-1)>tCNC)
			{
				continue;
			}
#endif

			VP rev_j(edge_j.size());
			reverse_copy(edge_j.begin(), edge_j.end(), rev_j.begin());

			uint key_ij = GenerateKey(PAIR_34, i, j);

			// For each edge k
			for (ushort k = 0; k < sz_k; ++k)
			{
				VP& edge_k = pk[k];
				ushort sz_ek = ushort(edge_k.size());

				Point& pkf = edge_k[0];
				Point& pkm = edge_k[sz_ek/2];
				Point& pkl = edge_k[sz_ek - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
				//CONSTRAINTS on position
				if (pkf.y > pif.y + _fThPosition)
				{
					//discard
					continue;
				}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
				// 添加限定条件2
				if(myselect2&&fabs(value4SixPoints(pif,pim,pil,pkf,pkm,pkl)-1)>tCNC)
				{
					continue;
				}
				// 添加限定条件3
				if(myselect3&&fabs(value4SixPoints(pjf,pjm,pjl,pkf,pkm,pkl)-1)>tCNC)
				{
					continue;
				}
#endif
				uint key_ik = GenerateKey(PAIR_23, k, i);

				// Find centers

				EllipseData data_ij, data_ik;

				if (data.count(key_ij) == 0)
				{
					//3,4 -> reverse 4

					GetFastCenter(edge_i, rev_j, data_ij);
					data.insert(pair<uint, EllipseData>(key_ij, data_ij));
				}
				else
				{
					data_ij = data.at(key_ij);
				}

				if (data.count(key_ik) == 0)
				{
					//3,2 -> reverse 3,2

					VP rev_k(edge_k.size());
					reverse_copy(edge_k.begin(), edge_k.end(), rev_k.begin());

					GetFastCenter(rev_i, rev_k, data_ik);

					data.insert(pair<uint, EllipseData>(key_ik, data_ik));
				}
				else
				{
					data_ik = data.at(key_ik);
				}


				// INVALID CENTERS
				if (!data_ij.isValid || !data_ik.isValid)
				{
					continue;
				}

#ifndef DISCARD_CONSTRAINT_CENTER
				// CONSTRAINT ON CENTERS
				if (ed2(data_ij.Cab, data_ik.Cab) > _fMaxCenterDistance2)
				{
					//discard
					continue;
				}
#endif
				// Find ellipse parameters
				Point2f center = GetCenterCoordinates(data_ij, data_ik);
				FindEllipses(center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses);
			}
		}

	}
};



void CNEllipseDetector::Triplets413(VVP& pi,
									   VVP& pj,
									   VVP& pk,
									   unordered_map<uint, EllipseData>& data,
									   vector<Ellipse>& ellipses
									   )
{
	ushort sz_i = ushort(pi.size());
	ushort sz_j = ushort(pj.size());
	ushort sz_k = ushort(pk.size());

	// For each edge i
	for (ushort i = 0; i < sz_i; ++i)
	{
		VP& edge_i = pi[i];
		ushort sz_ei = ushort(edge_i.size());

		Point& pif = edge_i[0];
		Point& pim = edge_i[sz_ei/2];
		Point& pil = edge_i[sz_ei - 1];

		VP rev_i(edge_i.size());
		reverse_copy(edge_i.begin(), edge_i.end(), rev_i.begin());

		// For each edge j
		for (ushort j = 0; j < sz_j; ++j)
		{
			VP& edge_j = pj[j];
			ushort sz_ej = ushort(edge_j.size());

			Point& pjf = edge_j[0];
			Point& pjm = edge_j[sz_ej/2];
			Point& pjl = edge_j[sz_ej - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
			//CONSTRAINTS on position
			if (pjl.y > pil.y + _fThPosition)  		//is below
			{
				//discard
				continue;
			}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
			// 添加限定条件 4se es1//pif,pim,pil,pjf,pjm,pjl pil,pim,pif,pjl,pjm,pjf pif,pim,pil,pjl,pjm,pjf
			if(myselect1&&fabs(value4SixPoints(T413)-1)>tCNC)
			{
				continue;
			}
#endif

			uint key_ij = GenerateKey(PAIR_14, j, i);

			// For each edge k
			for (ushort k = 0; k < sz_k; ++k)
			{
				VP& edge_k = pk[k];
				ushort sz_ek = ushort(edge_k.size());

				Point& pkf = edge_k[0];
				Point& pkm = edge_k[sz_ek/2];
				Point& pkl = edge_k[sz_ek - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
				//CONSTRAINTS on position
				if (pkl.x > pif.x + _fThPosition)
				{
					//discard
					continue;
				}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
				// 添加限定条件2
				if(myselect2&&fabs(value4SixPoints(pif,pim,pil,pkf,pkm,pkl)-1)>tCNC)
				{
					continue;
				}
				// 添加限定条件3
				if(myselect3&&fabs(value4SixPoints(pjf,pjm,pjl,pkf,pkm,pkl)-1)>tCNC)
				{
					continue;
				}
#endif
				uint key_ik = GenerateKey(PAIR_34, k, i);

				// Find centers

				EllipseData data_ij, data_ik;

				if (data.count(key_ij) == 0)
				{
					// 4,1 -> OK
					GetFastCenter(edge_i, edge_j, data_ij);
					data.insert(pair<uint, EllipseData>(key_ij, data_ij));
				}
				else
				{
					data_ij = data.at(key_ij);
				}

				if (data.count(key_ik) == 0)
				{
					// 4,3 -> reverse 4
					GetFastCenter(rev_i, edge_k, data_ik);
					data.insert(pair<uint, EllipseData>(key_ik, data_ik));
				}
				else
				{
					data_ik = data.at(key_ik);
				}

				// INVALID CENTERS
				if (!data_ij.isValid || !data_ik.isValid)
				{
					continue;
				}

#ifndef DISCARD_CONSTRAINT_CENTER
				// CONSTRAIN ON CENTERS
				if (ed2(data_ij.Cab, data_ik.Cab) > _fMaxCenterDistance2)
				{
					//discard
					continue;
				}
#endif
				// Find ellipse parameters
				Point2f center = GetCenterCoordinates(data_ij, data_ik);

				FindEllipses(center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses);

			}
		}
	}
};

void CNEllipseDetector::Triplets124_rect(VVP& pi,
	VVP& pj,
	VVP& pk,
	unordered_map<uint, EllipseData>& data,
	vector<Ellipse>& ellipses,
	Rectangle rect
)
{

	//cout << "in 124" << endl;
	// get arcs length
	ushort sz_i = ushort(pi.size());
	ushort sz_j = ushort(pj.size());
	ushort sz_k = ushort(pk.size());
	cout << "sz_i "<<sz_i << endl;
	cout << "sz_j "<<sz_j << endl;
	cout <<"sz_k " <<sz_k << endl;
	// For each edge i
	for (ushort i = 0; i < sz_i; ++i)
	{
		VP& edge_i = pi[i];
		ushort sz_ei = ushort(edge_i.size());

		Point& pif = edge_i[0];
		Point& pim = edge_i[sz_ei / 2];
		Point& pil = edge_i[sz_ei - 1];

		// 1,2 -> reverse 1, swap
		VP rev_i(edge_i.size());
		reverse_copy(edge_i.begin(), edge_i.end(), rev_i.begin());

		// For each edge j
		for (ushort j = 0; j < sz_j; ++j)
		{
			VP& edge_j = pj[j];
			ushort sz_ej = ushort(edge_j.size());

			Point& pjf = edge_j[0];
			Point& pjm = edge_j[sz_ej / 2];
			Point& pjl = edge_j[sz_ej - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
			// CONSTRAINTS on position
			if (pjl.x > pif.x + _fThPosition) //is right	
			{
				//discard
				cout << "i=" << i << "j=" << j << "pjl.x > pif.x + _fThPosition discard" << endl;
				continue;
			}
#endif
#ifndef DISCARD_CONSTRAINT_CNC
			// 添加限定条件  2se se1//pil,pim,pif,pjf,pjm,pjl pjf,pjm,pjl,pif,pim,pil
			//double cnc1 = fabs(value4SixPoints(T124) - 1);
			
			if (myselect1 && fabs(value4SixPoints(T124) - 1) > tCNC)
			{
				//cout <<cnc1 << endl;
				cout << "i=" << i << "j=" << j << "myselect1&&fabs(value4SixPoints(T124)-1)>tCNC" << endl;
				continue;
			}
#endif
			uint key_ij = GenerateKey(PAIR_12, i, j);

			//for each edge k
			for (ushort k = 0; k < sz_k; ++k)
			{
				VP& edge_k = pk[k];
				ushort sz_ek = ushort(edge_k.size());

				Point& pkf = edge_k[0];
				Point& pkm = edge_k[sz_ek / 2];
				Point& pkl = edge_k[sz_ek - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
				//CONSTRAINTS on position
				if (pkl.y < pil.y - _fThPosition)
				{
					//discard
					cout << "i=" << i << "j=" << j <<"k="<<k<<"pkl.y < pil.y - _fThPosition"<< endl;
					continue;
				}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
				// 添加限定条件2

				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{

					cout << "i=" << i << "j=" << j << "k=" << k << "myselect2&&fabs(value4SixPoints(pif,pim,pil,pkf,pkm,pkl)-1)>tCNC" << endl;
					continue;
				}
				// 添加限定条件3
				if (myselect3&&fabs(value4SixPoints(pjf, pjm, pjl, pkf, pkm, pkl) - 1) > tCNC)
				{
					cout << "i=" << i << "j=" << j << "k=" << k << "myselect3&&fabs(value4SixPoints(pjf,pjm,pjl,pkf,pkm,pkl)-1)>tCNC" << endl;
					continue;
				}
#endif

				uint key_ik = GenerateKey(PAIR_14, i, k);
				cout << "i=" << i << "j=" << j << "k=" << k << "pass this selection" << endl;

				// Find centers

				EllipseData data_ij, data_ik;

				// If the data for the pair i-j have not been computed yet
				if (data.count(key_ij) == 0)
				{
					//1,2 -> reverse 1, swap

					// Compute data!
					GetFastCenter(edge_j, rev_i, data_ij);
					// Insert computed data in the hash table
					data.insert(pair<uint, EllipseData>(key_ij, data_ij));
				}
				else
				{
					// Otherwise, just lookup the data in the hash table
					data_ij = data.at(key_ij);
				}

				// If the data for the pair i-k have not been computed yet
				if (data.count(key_ik) == 0)
				{
					//1,4 -> ok

					// Compute data!
					GetFastCenter(edge_i, edge_k, data_ik);
					// Insert computed data in the hash table
					data.insert(pair<uint, EllipseData>(key_ik, data_ik));
				}
				else
				{
					// Otherwise, just lookup the data in the hash table
					data_ik = data.at(key_ik);
				}

				// INVALID CENTERS
				if (!data_ij.isValid || !data_ik.isValid)
				{
					cout << "i=" << i << "j=" << j << "k=" << k << "is dropped cuz invalid center" << endl;
					continue;
				}

#ifndef DISCARD_CONSTRAINT_CENTER
				// Selection strategy - Step 3. See Sect [3.2.2] in the paper
				// The computed centers are not close enough
				if (ed2(data_ij.Cab, data_ik.Cab) > _fMaxCenterDistance2)
				{
					//discard
					cout << "ed2(data_ij.Cab, data_ik.Cab)" << ed2(data_ij.Cab, data_ik.Cab) << endl;
					cout << "_fMaxCenterDistance2" << _fMaxCenterDistance2 << endl;
					cout << "i=" << i << "j=" << j << "k=" << k << "cuz not closed enough" << endl;
					cout << "i=" << i << "j=" << j << "k=" << k << "is dropped cuz centers are not closed enough" << endl;
					continue;
				}
#endif
				// If all constraints of the selection strategy have been satisfied, 
				// we can start estimating the ellipse parameters

				// Find ellipse parameters

				// Get the coordinates of the center (xc, yc)
				Point2f center = GetCenterCoordinates(data_ij, data_ik);
				//cout << "in 124" << endl;
				cout << "i=" << i << "j=" << j << "k=" << k << "pass all the tests, now begins to find the params" << endl;
				// Find remaining paramters (A,B,rho)
				//cout << "in 124 one ell pass the CNC detect" << endl;
				FindEllipses(center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses, rect);


			}
		}
	}
};



void CNEllipseDetector::Triplets231_rect(VVP& pi,
	VVP& pj,
	VVP& pk,
	unordered_map<uint, EllipseData>& data,
	vector<Ellipse>& ellipses,
	Rectangle rect
)
{
	//cout << "in 231" << endl;
	ushort sz_i = ushort(pi.size());
	ushort sz_j = ushort(pj.size());
	ushort sz_k = ushort(pk.size());
	// For each edge i
	for (ushort i = 0; i < sz_i; ++i)
	{
		VP& edge_i = pi[i];
		ushort sz_ei = ushort(edge_i.size());

		Point& pif = edge_i[0];
		Point& pim = edge_i[sz_ei / 2];
		Point& pil = edge_i[sz_ei - 1];

		VP rev_i(edge_i.size());
		reverse_copy(edge_i.begin(), edge_i.end(), rev_i.begin());

		// For each edge j
		for (ushort j = 0; j < sz_j; ++j)
		{
			VP& edge_j = pj[j];
			ushort sz_ej = ushort(edge_j.size());

			Point& pjf = edge_j[0];
			Point& pjm = edge_j[sz_ej / 2];
			Point& pjl = edge_j[sz_ej - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
			// CONSTRAINTS on position
			if (pjf.y < pif.y - _fThPosition)
			{
				//discard
				continue;
			}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
			//添加约束
			// 添加限定条件 2es se3 //pif,pim,pil,pjf,pjm,pjl pil,pim,pif,pjf,pjm,pjl
			if (myselect1&&fabs(value4SixPoints(T231) - 1) > tCNC)
			{
				continue;
			}
#endif

			VP rev_j(edge_j.size());
			reverse_copy(edge_j.begin(), edge_j.end(), rev_j.begin());

			uint key_ij = GenerateKey(PAIR_23, i, j);

			// For each edge k
			for (ushort k = 0; k < sz_k; ++k)
			{
				VP& edge_k = pk[k];
				ushort sz_ek = ushort(edge_k.size());

				Point& pkf = edge_k[0];
				Point& pkm = edge_k[sz_ek / 2];
				Point& pkl = edge_k[sz_ek - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
				// CONSTRAINTS on position
				if (pkf.x < pil.x - _fThPosition)
				{
					//discard
					continue;
				}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
				// 添加限定条件2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
				// 添加限定条件3
				if (myselect3&&fabs(value4SixPoints(pjf, pjm, pjl, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
#endif
				uint key_ik = GenerateKey(PAIR_12, k, i);

				// Find centers

				EllipseData data_ij, data_ik;

				if (data.count(key_ij) == 0)
				{
					// 2,3 -> reverse 2,3

					GetFastCenter(rev_i, rev_j, data_ij);
					data.insert(pair<uint, EllipseData>(key_ij, data_ij));
				}
				else
				{
					data_ij = data.at(key_ij);
				}

				if (data.count(key_ik) == 0)
				{
					// 2,1 -> reverse 1
					VP rev_k(edge_k.size());
					reverse_copy(edge_k.begin(), edge_k.end(), rev_k.begin());

					GetFastCenter(edge_i, rev_k, data_ik);
					data.insert(pair<uint, EllipseData>(key_ik, data_ik));
				}
				else
				{
					data_ik = data.at(key_ik);
				}

				// INVALID CENTERS
				if (!data_ij.isValid || !data_ik.isValid)
				{
					continue;
				}

#ifndef DISCARD_CONSTRAINT_CENTER
				// CONSTRAINT ON CENTERS
				if (ed2(data_ij.Cab, data_ik.Cab) > _fMaxCenterDistance2)
				{
					//discard
					//cout << "ed2(data_ij.Cab, data_ik.Cab)" << ed2(data_ij.Cab, data_ik.Cab) << endl;
					//cout << "_fMaxCenterDistance2" << _fMaxCenterDistance2 << endl;
					//cout << "i=" << i << "j=" << j << "k=" << k << "cuz not closed enough" << endl;
					continue;
				}
#endif
				// Find ellipse parameters
				Point2f center = GetCenterCoordinates(data_ij, data_ik);
				//cout << "in 231 one ell pass the CNC detect" << endl;
				FindEllipses(center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses, rect);

			}
		}
	}
};


void CNEllipseDetector::Triplets342_rect(VVP& pi,
	VVP& pj,
	VVP& pk,
	unordered_map<uint, EllipseData>& data,
	vector<Ellipse>& ellipses,
	Rectangle rect
)
{
	//cout << "in 342" << endl;
	ushort sz_i = ushort(pi.size());
	ushort sz_j = ushort(pj.size());
	ushort sz_k = ushort(pk.size());

	// For each edge i
	for (ushort i = 0; i < sz_i; ++i)
	{
		VP& edge_i = pi[i];
		ushort sz_ei = ushort(edge_i.size());

		Point& pif = edge_i[0];
		Point& pim = edge_i[sz_ei / 2];
		Point& pil = edge_i[sz_ei - 1];

		VP rev_i(edge_i.size());
		reverse_copy(edge_i.begin(), edge_i.end(), rev_i.begin());

		// For each edge j
		for (ushort j = 0; j < sz_j; ++j)
		{
			VP& edge_j = pj[j];
			ushort sz_ej = ushort(edge_j.size());

			Point& pjf = edge_j[0];
			Point& pjm = edge_j[sz_ej / 2];
			Point& pjl = edge_j[sz_ej - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
			//CONSTRAINTS on position
			if (pjf.x < pil.x - _fThPosition) 		//is left
			{
				//discard
				continue;
			}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
			// 添加限定条件 3se se4 // pil,pim,pif,pjf,pjm,pjl pif,pim,pil,pjf,pjm,pjl
			if (myselect1&&fabs(value4SixPoints(T342) - 1) > tCNC)
			{
				continue;
			}
#endif

			VP rev_j(edge_j.size());
			reverse_copy(edge_j.begin(), edge_j.end(), rev_j.begin());

			uint key_ij = GenerateKey(PAIR_34, i, j);

			// For each edge k
			for (ushort k = 0; k < sz_k; ++k)
			{
				VP& edge_k = pk[k];
				ushort sz_ek = ushort(edge_k.size());

				Point& pkf = edge_k[0];
				Point& pkm = edge_k[sz_ek / 2];
				Point& pkl = edge_k[sz_ek - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
				//CONSTRAINTS on position
				if (pkf.y > pif.y + _fThPosition)
				{
					//discard
					continue;
				}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
				// 添加限定条件2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
				// 添加限定条件3
				if (myselect3&&fabs(value4SixPoints(pjf, pjm, pjl, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
#endif
				uint key_ik = GenerateKey(PAIR_23, k, i);

				// Find centers

				EllipseData data_ij, data_ik;

				if (data.count(key_ij) == 0)
				{
					//3,4 -> reverse 4

					GetFastCenter(edge_i, rev_j, data_ij);
					data.insert(pair<uint, EllipseData>(key_ij, data_ij));
				}
				else
				{
					data_ij = data.at(key_ij);
				}

				if (data.count(key_ik) == 0)
				{
					//3,2 -> reverse 3,2

					VP rev_k(edge_k.size());
					reverse_copy(edge_k.begin(), edge_k.end(), rev_k.begin());

					GetFastCenter(rev_i, rev_k, data_ik);

					data.insert(pair<uint, EllipseData>(key_ik, data_ik));
				}
				else
				{
					data_ik = data.at(key_ik);
				}


				// INVALID CENTERS
				if (!data_ij.isValid || !data_ik.isValid)
				{
					continue;
				}

#ifndef DISCARD_CONSTRAINT_CENTER
				// CONSTRAINT ON CENTERS
				if (ed2(data_ij.Cab, data_ik.Cab) > _fMaxCenterDistance2)
				{
					//discard
					//cout << "ed2(data_ij.Cab, data_ik.Cab)" << ed2(data_ij.Cab, data_ik.Cab) << endl;
					//cout << "_fMaxCenterDistance2" << _fMaxCenterDistance2 << endl;
					//cout << "i=" << i << "j=" << j << "k=" << k << "cuz not closed enough" << endl;
					continue;
				}
#endif
				// Find ellipse parameters
				Point2f center = GetCenterCoordinates(data_ij, data_ik);
				//cout << "in 342 one ell pass the CNC detect" << endl;
				FindEllipses(center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses, rect);
			}
		}

	}
};



void CNEllipseDetector::Triplets413_rect(VVP& pi,
	VVP& pj,
	VVP& pk,
	unordered_map<uint, EllipseData>& data,
	vector<Ellipse>& ellipses,
	Rectangle rect
)
{
	//cout << "in 413" << endl;
	//cout << "_fMaxCenterDistance2" <<_fMaxCenterDistance2 << endl;
	ushort sz_i = ushort(pi.size());
	ushort sz_j = ushort(pj.size());
	ushort sz_k = ushort(pk.size());
	// For each edge i
	for (ushort i = 0; i < sz_i; ++i)
	{
		VP& edge_i = pi[i];
		ushort sz_ei = ushort(edge_i.size());

		Point& pif = edge_i[0];
		Point& pim = edge_i[sz_ei / 2];
		Point& pil = edge_i[sz_ei - 1];

		VP rev_i(edge_i.size());
		reverse_copy(edge_i.begin(), edge_i.end(), rev_i.begin());

		// For each edge j
		for (ushort j = 0; j < sz_j; ++j)
		{
			VP& edge_j = pj[j];
			ushort sz_ej = ushort(edge_j.size());

			Point& pjf = edge_j[0];
			Point& pjm = edge_j[sz_ej / 2];
			Point& pjl = edge_j[sz_ej - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
			//CONSTRAINTS on position
			if (pjl.y > pil.y + _fThPosition)  		//is below
			{
				//discard
				//cout << "i=" << i << "j=" << j << "pjl.y > pil.y + _fThPosition" << endl;
				continue;
			}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
			// 添加限定条件 4se es1//pif,pim,pil,pjf,pjm,pjl pil,pim,pif,pjl,pjm,pjf pif,pim,pil,pjl,pjm,pjf
			double cnc1 = fabs(value4SixPoints(T413) - 1);
			if (myselect1&&cnc1 > tCNC)
			{
				//cout << "cnc1" <<cnc1<< endl;
				//cout << "i=" << i << "j=" << j << "myselect1&&fabs(value4SixPoints(T413)-1)>tCNC" << endl;
				continue;
			}
#endif

			uint key_ij = GenerateKey(PAIR_14, j, i);

			// For each edge k
			for (ushort k = 0; k < sz_k; ++k)
			{
				VP& edge_k = pk[k];
				ushort sz_ek = ushort(edge_k.size());

				Point& pkf = edge_k[0];
				Point& pkm = edge_k[sz_ek / 2];
				Point& pkl = edge_k[sz_ek - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
				//CONSTRAINTS on position
				if (pkl.x > pif.x + _fThPosition)
				{
					//discard
					//cout << "i=" << i << "j=" << j << "k=" << k << "pkl.x > pif.x + _fThPosition" << endl;
					continue;
				}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
				// 添加限定条件2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{
					//cout << "i=" << i << "j=" << j << "k=" << k << "myselect2&&fabs(value4SixPoints(pif,pim,pil,pkf,pkm,pkl)-1)>tCNC" << endl;
					continue;
				}
				// 添加限定条件3
				if (myselect3&&fabs(value4SixPoints(pjf, pjm, pjl, pkf, pkm, pkl) - 1) > tCNC)
				{
					//cout << "i=" << i << "j=" << j << "k=" << k << "myselect3&&fabs(value4SixPoints(pjf,pjm,pjl,pkf,pkm,pkl)-1)" << endl;
					continue;
				}
#endif
				uint key_ik = GenerateKey(PAIR_34, k, i);

				// Find centers

				EllipseData data_ij, data_ik;

				if (data.count(key_ij) == 0)
				{
					// 4,1 -> OK
					GetFastCenter(edge_i, edge_j, data_ij);
					data.insert(pair<uint, EllipseData>(key_ij, data_ij));
				}
				else
				{
					data_ij = data.at(key_ij);
				}

				if (data.count(key_ik) == 0)
				{
					// 4,3 -> reverse 4
					GetFastCenter(rev_i, edge_k, data_ik);
					data.insert(pair<uint, EllipseData>(key_ik, data_ik));
				}
				else
				{
					data_ik = data.at(key_ik);
				}

				// INVALID CENTERS
				if (!data_ij.isValid || !data_ik.isValid)
				{
					//cout << "i=" << i << "j=" << j << "k=" << k << "invalid centers" << endl;
					continue;
				}

#ifndef DISCARD_CONSTRAINT_CENTER
				// CONSTRAIN ON CENTERS

				//cout << "i=" << i << "j=" << j << "k=" << k << "center distance" << ed2(data_ij.Cab, data_ik.Cab) << endl;
				//cout << "_fMaxCenterDistance2" << _fMaxCenterDistance2 << endl;
				if (ed2(data_ij.Cab, data_ik.Cab) > _fMaxCenterDistance2)
				{
					//discard
					//cout << "ed2(data_ij.Cab, data_ik.Cab)" << ed2(data_ij.Cab, data_ik.Cab) << endl;
					//cout << "_fMaxCenterDistance2" << _fMaxCenterDistance2 << endl;
					//cout << "i=" << i << "j=" << j << "k=" << k << "cuz not closed enough" << endl;
					continue;
				}
#endif
				// Find ellipse parameters
				Point2f center = GetCenterCoordinates(data_ij, data_ik);
				//cout << "in 413" << endl;
				//cout << "i=" << i << "j=" << j << "k=" << k << "pass all the tests, now begins to find the params" << endl;
				//cout << "in 413 one ell pass the CNC detect" << endl;
				FindEllipses(center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses, rect);


			}
		}
	}
};



void CNEllipseDetector::Triplets124_ArcsIn(VVP& pi,
	VVP& pj,
	VVP& pk,
	unordered_map<uint, EllipseData>& data,
	vector<Ellipse>& ellipses
)
{
	myselect1 = true;
	myselect2 = false;
	myselect3 = false;
	// get arcs length
	ushort sz_i = ushort(pi.size());
	ushort sz_j = ushort(pj.size());
	ushort sz_k = ushort(pk.size());

	// For each edge i
	for (ushort i = 0; i < sz_i; ++i)
	{
		VP& edge_i = pi[i];
		ushort sz_ei = ushort(edge_i.size());

		Point& pif = edge_i[0];
		Point& pim = edge_i[sz_ei / 2];
		Point& pil = edge_i[sz_ei - 1];

		// 1,2 -> reverse 1, swap
		VP rev_i(edge_i.size());
		reverse_copy(edge_i.begin(), edge_i.end(), rev_i.begin());

		// For each edge j
		for (ushort j = 0; j < sz_j; ++j)
		{
			VP& edge_j = pj[j];
			ushort sz_ej = ushort(edge_j.size());

			Point& pjf = edge_j[0];
			Point& pjm = edge_j[sz_ej / 2];
			Point& pjl = edge_j[sz_ej - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
			// CONSTRAINTS on position
			if (pjl.x > pif.x + _fThPosition) //is right	
			{
				//discard
				continue;
			}
#endif
#ifndef DISCARD_CONSTRAINT_CNC
			// 添加限定条件  2se se1//pil,pim,pif,pjf,pjm,pjl pjf,pjm,pjl,pif,pim,pil
			double cnc1 = fabs(value4SixPoints(T124) - 1);
			//cout << "i=" << i << "j=" << j << endl;
			//cout <<"cncT124"<<cnc1<< endl;
			if (myselect1&&fabs(value4SixPoints(T124) - 1) > tCNC)
			{
				continue;
			}
#endif
			uint key_ij = GenerateKey(PAIR_12, i, j);

			//for each edge k
			for (ushort k = 0; k < sz_k; ++k)
			{
				VP& edge_k = pk[k];
				ushort sz_ek = ushort(edge_k.size());

				Point& pkf = edge_k[0];
				Point& pkm = edge_k[sz_ek / 2];
				Point& pkl = edge_k[sz_ek - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
				//CONSTRAINTS on position
				if (pkl.y < pil.y - _fThPosition)
				{
					//discard
					continue;
				}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
				// 添加限定条件2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
				// 添加限定条件3
				if (myselect3&&fabs(value4SixPoints(pjf, pjm, pjl, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}

#endif

				uint key_ik = GenerateKey(PAIR_14, i, k);

				// Find centers

				EllipseData data_ij, data_ik;

				// If the data for the pair i-j have not been computed yet
				if (data.count(key_ij) == 0)
				{
					//1,2 -> reverse 1, swap

					// Compute data!
					GetFastCenter(edge_j, rev_i, data_ij);
					// Insert computed data in the hash table
					data.insert(pair<uint, EllipseData>(key_ij, data_ij));
				}
				else
				{
					// Otherwise, just lookup the data in the hash table
					data_ij = data.at(key_ij);
				}

				// If the data for the pair i-k have not been computed yet
				if (data.count(key_ik) == 0)
				{
					//1,4 -> ok

					// Compute data!
					GetFastCenter(edge_i, edge_k, data_ik);
					// Insert computed data in the hash table
					data.insert(pair<uint, EllipseData>(key_ik, data_ik));
				}
				else
				{
					// Otherwise, just lookup the data in the hash table
					data_ik = data.at(key_ik);
				}

				// INVALID CENTERS
				if (!data_ij.isValid || !data_ik.isValid)
				{
					continue;
				}

#ifndef DISCARD_CONSTRAINT_CENTER
				// Selection strategy - Step 3. See Sect [3.2.2] in the paper
				// The computed centers are not close enough
				if (ed2(data_ij.Cab, data_ik.Cab) > _fMaxCenterDistance2)
				{
					//discard
					continue;
				}
#endif
				// If all constraints of the selection strategy have been satisfied, 
				// we can start estimating the ellipse parameters

				// Find ellipse parameters

				// Get the coordinates of the center (xc, yc)
				Point2f center = GetCenterCoordinates(data_ij, data_ik);

				// Find remaining paramters (A,B,rho)
				//cout << "i=" << i << "j=" << j << "k=" << k << "pass all the detection" << endl;
				bool flag = FindEllipses(1, 2, 4, i, j, k, center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses);
			}
		}
	}
};

void CNEllipseDetector::Triplets231_ArcsIn(VVP& pi,
	VVP& pj,
	VVP& pk,
	unordered_map<uint, EllipseData>& data,
	vector<Ellipse>& ellipses
)
{
	myselect1 = true;
	myselect2 = false;
	myselect3 = false;
	ushort sz_i = ushort(pi.size());
	ushort sz_j = ushort(pj.size());
	ushort sz_k = ushort(pk.size());

	// For each edge i
	for (ushort i = 0; i < sz_i; ++i)
	{
		VP& edge_i = pi[i];
		ushort sz_ei = ushort(edge_i.size());

		Point& pif = edge_i[0];
		Point& pim = edge_i[sz_ei / 2];
		Point& pil = edge_i[sz_ei - 1];

		VP rev_i(edge_i.size());
		reverse_copy(edge_i.begin(), edge_i.end(), rev_i.begin());

		// For each edge j
		for (ushort j = 0; j < sz_j; ++j)
		{
			VP& edge_j = pj[j];
			ushort sz_ej = ushort(edge_j.size());

			Point& pjf = edge_j[0];
			Point& pjm = edge_j[sz_ej / 2];
			Point& pjl = edge_j[sz_ej - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
			// CONSTRAINTS on position
			if (pjf.y < pif.y - _fThPosition)
			{
				//discard
				continue;
			}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
			//添加约束
			// 添加限定条件 2es se3 //pif,pim,pil,pjf,pjm,pjl pil,pim,pif,pjf,pjm,pjl
			if (myselect1&&fabs(value4SixPoints(T231) - 1) > tCNC)
			{
				continue;
			}
#endif

			VP rev_j(edge_j.size());
			reverse_copy(edge_j.begin(), edge_j.end(), rev_j.begin());

			uint key_ij = GenerateKey(PAIR_23, i, j);

			// For each edge k
			for (ushort k = 0; k < sz_k; ++k)
			{
				VP& edge_k = pk[k];
				ushort sz_ek = ushort(edge_k.size());

				Point& pkf = edge_k[0];
				Point& pkm = edge_k[sz_ek / 2];
				Point& pkl = edge_k[sz_ek - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
				// CONSTRAINTS on position
				if (pkf.x < pil.x - _fThPosition)
				{
					//discard
					continue;
				}
#endif

#ifndef DISCARD_CONSTRAINT_CNC

				// 添加限定条件2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
				// 添加限定条件3
				if (myselect3&&fabs(value4SixPoints(pjf, pjm, pjl, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}

#endif
				uint key_ik = GenerateKey(PAIR_12, k, i);

				// Find centers

				EllipseData data_ij, data_ik;

				if (data.count(key_ij) == 0)
				{
					// 2,3 -> reverse 2,3

					GetFastCenter(rev_i, rev_j, data_ij);
					data.insert(pair<uint, EllipseData>(key_ij, data_ij));
				}
				else
				{
					data_ij = data.at(key_ij);
				}

				if (data.count(key_ik) == 0)
				{
					// 2,1 -> reverse 1
					VP rev_k(edge_k.size());
					reverse_copy(edge_k.begin(), edge_k.end(), rev_k.begin());

					GetFastCenter(edge_i, rev_k, data_ik);
					data.insert(pair<uint, EllipseData>(key_ik, data_ik));
				}
				else
				{
					data_ik = data.at(key_ik);
				}

				// INVALID CENTERS
				if (!data_ij.isValid || !data_ik.isValid)
				{
					continue;
				}

#ifndef DISCARD_CONSTRAINT_CENTER
				// CONSTRAINT ON CENTERS
				if (ed2(data_ij.Cab, data_ik.Cab) > _fMaxCenterDistance2)
				{
					//discard
					continue;
				}
#endif
				// Find ellipse parameters
				Point2f center = GetCenterCoordinates(data_ij, data_ik);
				//cout << "in 231" << endl;
				//cout << "i=" << i << "j=" << j << "k=" << k << "pass all the detection" << endl;
				FindEllipses(2, 3, 1, i, j, k, center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses);

			}
		}
	}
};

void CNEllipseDetector::Triplets342_ArcsIn(VVP& pi,
	VVP& pj,
	VVP& pk,
	unordered_map<uint, EllipseData>& data,
	vector<Ellipse>& ellipses
)
{
	myselect1 = true;
	myselect2 = false;
	myselect3 = false;
	ushort sz_i = ushort(pi.size());
	ushort sz_j = ushort(pj.size());
	ushort sz_k = ushort(pk.size());

	// For each edge i
	for (ushort i = 0; i < sz_i; ++i)
	{
		VP& edge_i = pi[i];
		ushort sz_ei = ushort(edge_i.size());

		Point& pif = edge_i[0];
		Point& pim = edge_i[sz_ei / 2];
		Point& pil = edge_i[sz_ei - 1];

		VP rev_i(edge_i.size());
		reverse_copy(edge_i.begin(), edge_i.end(), rev_i.begin());

		// For each edge j
		for (ushort j = 0; j < sz_j; ++j)
		{
			VP& edge_j = pj[j];
			ushort sz_ej = ushort(edge_j.size());

			Point& pjf = edge_j[0];
			Point& pjm = edge_j[sz_ej / 2];
			Point& pjl = edge_j[sz_ej - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
			//CONSTRAINTS on position
			if (pjf.x < pil.x - _fThPosition) 		//is left
			{
				//discard
				continue;
			}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
			// 添加限定条件 3se se4 // pil,pim,pif,pjf,pjm,pjl pif,pim,pil,pjf,pjm,pjl
			if (myselect1&&fabs(value4SixPoints(T342) - 1) > tCNC)
			{
				continue;
			}
#endif

			VP rev_j(edge_j.size());
			reverse_copy(edge_j.begin(), edge_j.end(), rev_j.begin());

			uint key_ij = GenerateKey(PAIR_34, i, j);

			// For each edge k
			for (ushort k = 0; k < sz_k; ++k)
			{
				VP& edge_k = pk[k];
				ushort sz_ek = ushort(edge_k.size());

				Point& pkf = edge_k[0];
				Point& pkm = edge_k[sz_ek / 2];
				Point& pkl = edge_k[sz_ek - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
				//CONSTRAINTS on position
				if (pkf.y > pif.y + _fThPosition)
				{
					//discard
					continue;
				}
#endif

#ifndef DISCARD_CONSTRAINT_CNC

				// 添加限定条件2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
				// 添加限定条件3
				if (myselect3&&fabs(value4SixPoints(pjf, pjm, pjl, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}

#endif
				uint key_ik = GenerateKey(PAIR_23, k, i);

				// Find centers

				EllipseData data_ij, data_ik;

				if (data.count(key_ij) == 0)
				{
					//3,4 -> reverse 4

					GetFastCenter(edge_i, rev_j, data_ij);
					data.insert(pair<uint, EllipseData>(key_ij, data_ij));
				}
				else
				{
					data_ij = data.at(key_ij);
				}

				if (data.count(key_ik) == 0)
				{
					//3,2 -> reverse 3,2

					VP rev_k(edge_k.size());
					reverse_copy(edge_k.begin(), edge_k.end(), rev_k.begin());

					GetFastCenter(rev_i, rev_k, data_ik);

					data.insert(pair<uint, EllipseData>(key_ik, data_ik));
				}
				else
				{
					data_ik = data.at(key_ik);
				}


				// INVALID CENTERS
				if (!data_ij.isValid || !data_ik.isValid)
				{
					continue;
				}

#ifndef DISCARD_CONSTRAINT_CENTER
				// CONSTRAINT ON CENTERS
				if (ed2(data_ij.Cab, data_ik.Cab) > _fMaxCenterDistance2)
				{
					//discard
					continue;
				}
#endif
				// Find ellipse parameters
				Point2f center = GetCenterCoordinates(data_ij, data_ik);
				//cout << "in 342" << endl;
				//cout << "i=" << i << "j=" << j << "k=" << k << "pass all the detection" << endl;
				FindEllipses(3, 4, 2, i, j, k, center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses);


			}
		}

	}
};


void CNEllipseDetector::Triplets413_ArcsIn(VVP& pi,
	VVP& pj,
	VVP& pk,
	unordered_map<uint, EllipseData>& data,
	vector<Ellipse>& ellipses
)
{
	myselect1 = true;
	myselect2 = false;
	myselect3 = false;
	ushort sz_i = ushort(pi.size());
	ushort sz_j = ushort(pj.size());
	ushort sz_k = ushort(pk.size());

	// For each edge i
	for (ushort i = 0; i < sz_i; ++i)
	{
		VP& edge_i = pi[i];
		ushort sz_ei = ushort(edge_i.size());

		Point& pif = edge_i[0];
		Point& pim = edge_i[sz_ei / 2];
		Point& pil = edge_i[sz_ei - 1];

		VP rev_i(edge_i.size());
		reverse_copy(edge_i.begin(), edge_i.end(), rev_i.begin());

		// For each edge j
		for (ushort j = 0; j < sz_j; ++j)
		{
			VP& edge_j = pj[j];
			ushort sz_ej = ushort(edge_j.size());

			Point& pjf = edge_j[0];
			Point& pjm = edge_j[sz_ej / 2];
			Point& pjl = edge_j[sz_ej - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
			//CONSTRAINTS on position
			if (pjl.y > pil.y + _fThPosition)  		//is below
			{
				//discard
				continue;
			}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
			// 添加限定条件 4se es1//pif,pim,pil,pjf,pjm,pjl pil,pim,pif,pjl,pjm,pjf pif,pim,pil,pjl,pjm,pjf
			if (myselect1&&fabs(value4SixPoints(T413) - 1) > tCNC)
			{
				continue;
			}
#endif


			uint key_ij = GenerateKey(PAIR_14, j, i);

			// For each edge k
			for (ushort k = 0; k < sz_k; ++k)
			{
				VP& edge_k = pk[k];
				ushort sz_ek = ushort(edge_k.size());

				Point& pkf = edge_k[0];
				Point& pkm = edge_k[sz_ek / 2];
				Point& pkl = edge_k[sz_ek - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
				//CONSTRAINTS on position
				if (pkl.x > pif.x + _fThPosition)
				{
					//discard
					continue;
				}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
				// 添加限定条件2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
				// 添加限定条件3
				if (myselect3&&fabs(value4SixPoints(pjf, pjm, pjl, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}

#endif
				uint key_ik = GenerateKey(PAIR_34, k, i);

				// Find centers

				EllipseData data_ij, data_ik;

				if (data.count(key_ij) == 0)
				{
					// 4,1 -> OK
					GetFastCenter(edge_i, edge_j, data_ij);
					data.insert(pair<uint, EllipseData>(key_ij, data_ij));
				}
				else
				{
					data_ij = data.at(key_ij);
				}

				if (data.count(key_ik) == 0)
				{
					// 4,3 -> reverse 4
					GetFastCenter(rev_i, edge_k, data_ik);
					data.insert(pair<uint, EllipseData>(key_ik, data_ik));
				}
				else
				{
					data_ik = data.at(key_ik);
				}

				// INVALID CENTERS
				if (!data_ij.isValid || !data_ik.isValid)
				{
					continue;
				}

#ifndef DISCARD_CONSTRAINT_CENTER
				// CONSTRAIN ON CENTERS
				//cout << "i=" << i << "j=" << j << "k=" << k <<"center distance"<<ed2(data_ij.Cab, data_ik.Cab) << endl;
				//cout << "_fMaxCenterDistance2" << _fMaxCenterDistance2 << endl;
				if (ed2(data_ij.Cab, data_ik.Cab) > _fMaxCenterDistance2)
				{
					//discard
					continue;
				}
#endif
				// Find ellipse parameters
				Point2f center = GetCenterCoordinates(data_ij, data_ik);
				//cout << "in 413" << endl;
				//cout << "i=" << i << "j=" << j << "k=" << k << "pass all the detection" << endl;
				FindEllipses(4, 1, 3, i, j, k, center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses);

			}
		}
	}
};

void CNEllipseDetector::Triplets124_3(VVP& pi,
	VVP& pj,
	VVP& pk,
	unordered_map<uint, EllipseData>& data,
	vector<Ellipse>& ellipses
)
{
	// get arcs length
	myselect1 = true;
	myselect2 = true;
	myselect3 = true;
	ushort sz_i = ushort(pi.size());
	ushort sz_j = ushort(pj.size());
	ushort sz_k = ushort(pk.size());

	// For each edge i
	for (ushort i = 0; i < sz_i; ++i)
	{
		VP& edge_i = pi[i];
		ushort sz_ei = ushort(edge_i.size());

		Point& pif = edge_i[0];
		Point& pim = edge_i[sz_ei / 2];
		Point& pil = edge_i[sz_ei - 1];

		// 1,2 -> reverse 1, swap
		VP rev_i(edge_i.size());
		reverse_copy(edge_i.begin(), edge_i.end(), rev_i.begin());

		// For each edge j
		for (ushort j = 0; j < sz_j; ++j)
		{
			VP& edge_j = pj[j];
			ushort sz_ej = ushort(edge_j.size());

			Point& pjf = edge_j[0];
			Point& pjm = edge_j[sz_ej / 2];
			Point& pjl = edge_j[sz_ej - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
			// CONSTRAINTS on position
			if (pjl.x > pif.x + _fThPosition) //is right	
			{
				//discard
				continue;
			}
#endif
#ifndef DISCARD_CONSTRAINT_CNC
			// 添加限定条件  2se se1//pil,pim,pif,pjf,pjm,pjl pjf,pjm,pjl,pif,pim,pil
			double cnc1 = fabs(value4SixPoints(T124) - 1);
			//cout << "i=" << i << "j=" << j << endl;
			//cout <<"cncT124"<<cnc1<< endl;
			if (myselect1&&fabs(value4SixPoints(T124) - 1) > tCNC)
			{
				continue;
			}
#endif
			uint key_ij = GenerateKey(PAIR_12, i, j);

			//for each edge k
			for (ushort k = 0; k < sz_k; ++k)
			{
				VP& edge_k = pk[k];
				ushort sz_ek = ushort(edge_k.size());

				Point& pkf = edge_k[0];
				Point& pkm = edge_k[sz_ek / 2];
				Point& pkl = edge_k[sz_ek - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
				//CONSTRAINTS on position
				if (pkl.y < pil.y - _fThPosition)
				{
					//discard
					continue;
				}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
				// 添加限定条件2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
				// 添加限定条件3
				if (myselect3&&fabs(value4SixPoints(pjf, pjm, pjl, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}

#endif

				uint key_ik = GenerateKey(PAIR_14, i, k);

				// Find centers

				EllipseData data_ij, data_ik;

				// If the data for the pair i-j have not been computed yet
				if (data.count(key_ij) == 0)
				{
					//1,2 -> reverse 1, swap

					// Compute data!
					GetFastCenter(edge_j, rev_i, data_ij);
					// Insert computed data in the hash table
					data.insert(pair<uint, EllipseData>(key_ij, data_ij));
				}
				else
				{
					// Otherwise, just lookup the data in the hash table
					data_ij = data.at(key_ij);
				}

				// If the data for the pair i-k have not been computed yet
				if (data.count(key_ik) == 0)
				{
					//1,4 -> ok

					// Compute data!
					GetFastCenter(edge_i, edge_k, data_ik);
					// Insert computed data in the hash table
					data.insert(pair<uint, EllipseData>(key_ik, data_ik));
				}
				else
				{
					// Otherwise, just lookup the data in the hash table
					data_ik = data.at(key_ik);
				}

				// INVALID CENTERS
				if (!data_ij.isValid || !data_ik.isValid)
				{
					continue;
				}

#ifndef DISCARD_CONSTRAINT_CENTER
				// Selection strategy - Step 3. See Sect [3.2.2] in the paper
				// The computed centers are not close enough
				if (ed2(data_ij.Cab, data_ik.Cab) > _fMaxCenterDistance2)
				{
					//discard
					continue;
				}
#endif
				// If all constraints of the selection strategy have been satisfied, 
				// we can start estimating the ellipse parameters

				// Find ellipse parameters

				// Get the coordinates of the center (xc, yc)
				Point2f center = GetCenterCoordinates(data_ij, data_ik);

				// Find remaining paramters (A,B,rho)
				FindEllipses(center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses);
			}
		}
	}
};

void CNEllipseDetector::Triplets231_3(VVP& pi,
	VVP& pj,
	VVP& pk,
	unordered_map<uint, EllipseData>& data,
	vector<Ellipse>& ellipses
)
{
	myselect1 = true;
	myselect2 = true;
	myselect3 = true;
	ushort sz_i = ushort(pi.size());
	ushort sz_j = ushort(pj.size());
	ushort sz_k = ushort(pk.size());

	// For each edge i
	for (ushort i = 0; i < sz_i; ++i)
	{
		VP& edge_i = pi[i];
		ushort sz_ei = ushort(edge_i.size());

		Point& pif = edge_i[0];
		Point& pim = edge_i[sz_ei / 2];
		Point& pil = edge_i[sz_ei - 1];

		VP rev_i(edge_i.size());
		reverse_copy(edge_i.begin(), edge_i.end(), rev_i.begin());

		// For each edge j
		for (ushort j = 0; j < sz_j; ++j)
		{
			VP& edge_j = pj[j];
			ushort sz_ej = ushort(edge_j.size());

			Point& pjf = edge_j[0];
			Point& pjm = edge_j[sz_ej / 2];
			Point& pjl = edge_j[sz_ej - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
			// CONSTRAINTS on position
			if (pjf.y < pif.y - _fThPosition)
			{
				//discard
				continue;
			}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
			//添加约束
			// 添加限定条件 2es se3 //pif,pim,pil,pjf,pjm,pjl pil,pim,pif,pjf,pjm,pjl
			if (myselect1&&fabs(value4SixPoints(T231) - 1) > tCNC)
			{
				continue;
			}
#endif

			VP rev_j(edge_j.size());
			reverse_copy(edge_j.begin(), edge_j.end(), rev_j.begin());

			uint key_ij = GenerateKey(PAIR_23, i, j);

			// For each edge k
			for (ushort k = 0; k < sz_k; ++k)
			{
				VP& edge_k = pk[k];
				ushort sz_ek = ushort(edge_k.size());

				Point& pkf = edge_k[0];
				Point& pkm = edge_k[sz_ek / 2];
				Point& pkl = edge_k[sz_ek - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
				// CONSTRAINTS on position
				if (pkf.x < pil.x - _fThPosition)
				{
					//discard
					continue;
				}
#endif

#ifndef DISCARD_CONSTRAINT_CNC

				// 添加限定条件2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
				// 添加限定条件3
				if (myselect3&&fabs(value4SixPoints(pjf, pjm, pjl, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}

#endif
				uint key_ik = GenerateKey(PAIR_12, k, i);

				// Find centers

				EllipseData data_ij, data_ik;

				if (data.count(key_ij) == 0)
				{
					// 2,3 -> reverse 2,3

					GetFastCenter(rev_i, rev_j, data_ij);
					data.insert(pair<uint, EllipseData>(key_ij, data_ij));
				}
				else
				{
					data_ij = data.at(key_ij);
				}

				if (data.count(key_ik) == 0)
				{
					// 2,1 -> reverse 1
					VP rev_k(edge_k.size());
					reverse_copy(edge_k.begin(), edge_k.end(), rev_k.begin());

					GetFastCenter(edge_i, rev_k, data_ik);
					data.insert(pair<uint, EllipseData>(key_ik, data_ik));
				}
				else
				{
					data_ik = data.at(key_ik);
				}

				// INVALID CENTERS
				if (!data_ij.isValid || !data_ik.isValid)
				{
					continue;
				}

#ifndef DISCARD_CONSTRAINT_CENTER
				// CONSTRAINT ON CENTERS
				if (ed2(data_ij.Cab, data_ik.Cab) > _fMaxCenterDistance2)
				{
					//discard
					continue;
				}
#endif
				// Find ellipse parameters
				Point2f center = GetCenterCoordinates(data_ij, data_ik);
				//cout << "in 231" << endl;
				//cout << "i=" << i << "j=" << j << "k=" << k << "pass all the detection" << endl;
				FindEllipses(center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses);

			}
		}
	}
};

void CNEllipseDetector::Triplets342_3(VVP& pi,
	VVP& pj,
	VVP& pk,
	unordered_map<uint, EllipseData>& data,
	vector<Ellipse>& ellipses
)
{
	myselect1 = true;
	myselect2 = true;
	myselect3 = true;
	ushort sz_i = ushort(pi.size());
	ushort sz_j = ushort(pj.size());
	ushort sz_k = ushort(pk.size());

	// For each edge i
	for (ushort i = 0; i < sz_i; ++i)
	{
		VP& edge_i = pi[i];
		ushort sz_ei = ushort(edge_i.size());

		Point& pif = edge_i[0];
		Point& pim = edge_i[sz_ei / 2];
		Point& pil = edge_i[sz_ei - 1];

		VP rev_i(edge_i.size());
		reverse_copy(edge_i.begin(), edge_i.end(), rev_i.begin());

		// For each edge j
		for (ushort j = 0; j < sz_j; ++j)
		{
			VP& edge_j = pj[j];
			ushort sz_ej = ushort(edge_j.size());

			Point& pjf = edge_j[0];
			Point& pjm = edge_j[sz_ej / 2];
			Point& pjl = edge_j[sz_ej - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
			//CONSTRAINTS on position
			if (pjf.x < pil.x - _fThPosition) 		//is left
			{
				//discard
				continue;
			}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
			// 添加限定条件 3se se4 // pil,pim,pif,pjf,pjm,pjl pif,pim,pil,pjf,pjm,pjl
			if (myselect1&&fabs(value4SixPoints(T342) - 1) > tCNC)
			{
				continue;
			}
#endif

			VP rev_j(edge_j.size());
			reverse_copy(edge_j.begin(), edge_j.end(), rev_j.begin());

			uint key_ij = GenerateKey(PAIR_34, i, j);

			// For each edge k
			for (ushort k = 0; k < sz_k; ++k)
			{
				VP& edge_k = pk[k];
				ushort sz_ek = ushort(edge_k.size());

				Point& pkf = edge_k[0];
				Point& pkm = edge_k[sz_ek / 2];
				Point& pkl = edge_k[sz_ek - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
				//CONSTRAINTS on position
				if (pkf.y > pif.y + _fThPosition)
				{
					//discard
					continue;
				}
#endif

#ifndef DISCARD_CONSTRAINT_CNC

				// 添加限定条件2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
				// 添加限定条件3
				if (myselect3&&fabs(value4SixPoints(pjf, pjm, pjl, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}

#endif
				uint key_ik = GenerateKey(PAIR_23, k, i);

				// Find centers

				EllipseData data_ij, data_ik;

				if (data.count(key_ij) == 0)
				{
					//3,4 -> reverse 4

					GetFastCenter(edge_i, rev_j, data_ij);
					data.insert(pair<uint, EllipseData>(key_ij, data_ij));
				}
				else
				{
					data_ij = data.at(key_ij);
				}

				if (data.count(key_ik) == 0)
				{
					//3,2 -> reverse 3,2

					VP rev_k(edge_k.size());
					reverse_copy(edge_k.begin(), edge_k.end(), rev_k.begin());

					GetFastCenter(rev_i, rev_k, data_ik);

					data.insert(pair<uint, EllipseData>(key_ik, data_ik));
				}
				else
				{
					data_ik = data.at(key_ik);
				}


				// INVALID CENTERS
				if (!data_ij.isValid || !data_ik.isValid)
				{
					continue;
				}

#ifndef DISCARD_CONSTRAINT_CENTER
				// CONSTRAINT ON CENTERS
				if (ed2(data_ij.Cab, data_ik.Cab) > _fMaxCenterDistance2)
				{
					//discard
					continue;
				}
#endif
				// Find ellipse parameters
				Point2f center = GetCenterCoordinates(data_ij, data_ik);
				//cout << "in 342" << endl;
				//cout << "i=" << i << "j=" << j << "k=" << k << "pass all the detection" << endl;
				FindEllipses(center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses);


			}
		}

	}
};

void CNEllipseDetector::Triplets413_3(VVP& pi,
	VVP& pj,
	VVP& pk,
	unordered_map<uint, EllipseData>& data,
	vector<Ellipse>& ellipses
)
{
	myselect1 = true;
	myselect2 = true;
	myselect3 = true;
	ushort sz_i = ushort(pi.size());
	ushort sz_j = ushort(pj.size());
	ushort sz_k = ushort(pk.size());

	// For each edge i
	for (ushort i = 0; i < sz_i; ++i)
	{
		VP& edge_i = pi[i];
		ushort sz_ei = ushort(edge_i.size());

		Point& pif = edge_i[0];
		Point& pim = edge_i[sz_ei / 2];
		Point& pil = edge_i[sz_ei - 1];

		VP rev_i(edge_i.size());
		reverse_copy(edge_i.begin(), edge_i.end(), rev_i.begin());

		// For each edge j
		for (ushort j = 0; j < sz_j; ++j)
		{
			VP& edge_j = pj[j];
			ushort sz_ej = ushort(edge_j.size());

			Point& pjf = edge_j[0];
			Point& pjm = edge_j[sz_ej / 2];
			Point& pjl = edge_j[sz_ej - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
			//CONSTRAINTS on position
			if (pjl.y > pil.y + _fThPosition)  		//is below
			{
				//discard
				continue;
			}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
			// 添加限定条件 4se es1//pif,pim,pil,pjf,pjm,pjl pil,pim,pif,pjl,pjm,pjf pif,pim,pil,pjl,pjm,pjf
			if (myselect1&&fabs(value4SixPoints(T413) - 1) > tCNC)
			{
				continue;
			}
#endif


			uint key_ij = GenerateKey(PAIR_14, j, i);

			// For each edge k
			for (ushort k = 0; k < sz_k; ++k)
			{
				VP& edge_k = pk[k];
				ushort sz_ek = ushort(edge_k.size());

				Point& pkf = edge_k[0];
				Point& pkm = edge_k[sz_ek / 2];
				Point& pkl = edge_k[sz_ek - 1];

#ifndef DISCARD_CONSTRAINT_POSITION
				//CONSTRAINTS on position
				if (pkl.x > pif.x + _fThPosition)
				{
					//discard
					continue;
				}
#endif

#ifndef DISCARD_CONSTRAINT_CNC
				// 添加限定条件2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
				// 添加限定条件3
				if (myselect3&&fabs(value4SixPoints(pjf, pjm, pjl, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}

#endif
				uint key_ik = GenerateKey(PAIR_34, k, i);

				// Find centers

				EllipseData data_ij, data_ik;

				if (data.count(key_ij) == 0)
				{
					// 4,1 -> OK
					GetFastCenter(edge_i, edge_j, data_ij);
					data.insert(pair<uint, EllipseData>(key_ij, data_ij));
				}
				else
				{
					data_ij = data.at(key_ij);
				}

				if (data.count(key_ik) == 0)
				{
					// 4,3 -> reverse 4
					GetFastCenter(rev_i, edge_k, data_ik);
					data.insert(pair<uint, EllipseData>(key_ik, data_ik));
				}
				else
				{
					data_ik = data.at(key_ik);
				}

				// INVALID CENTERS
				if (!data_ij.isValid || !data_ik.isValid)
				{
					continue;
				}

#ifndef DISCARD_CONSTRAINT_CENTER
				// CONSTRAIN ON CENTERS
				//cout << "i=" << i << "j=" << j << "k=" << k <<"center distance"<<ed2(data_ij.Cab, data_ik.Cab) << endl;
				//cout << "_fMaxCenterDistance2" << _fMaxCenterDistance2 << endl;
				if (ed2(data_ij.Cab, data_ik.Cab) > _fMaxCenterDistance2)
				{
					//discard
					continue;
				}
#endif
				// Find ellipse parameters
				Point2f center = GetCenterCoordinates(data_ij, data_ik);
				//cout << "in 413" << endl;
				//cout << "i=" << i << "j=" << j << "k=" << k << "pass all the detection" << endl;
				FindEllipses(center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses);

			}
		}
	}
};


void CNEllipseDetector::RemoveShortEdges(Mat1b& edges, Mat1b& clean)
{
	VVP contours;

	// Labeling and contraints on length
	Labeling(edges, contours, _iMinEdgeLength);

	int iContoursSize = contours.size();
	for (int i = 0; i < iContoursSize; ++i)
	{
		VP& edge = contours[i];
		unsigned szEdge = edge.size();

		// Constraint on axes aspect ratio
		RotatedRect oriented = minAreaRect(edge);
		if (oriented.size.width < _fMinOrientedRectSide ||
			oriented.size.height < _fMinOrientedRectSide ||
			oriented.size.width > oriented.size.height * _fMaxRectAxesRatio ||
			oriented.size.height > oriented.size.width * _fMaxRectAxesRatio)
		{
			continue;
		}

		for (unsigned j = 0; j < szEdge; ++j)
		{
			clean(edge[j]) = (uchar)255;
		}
	}
}



void CNEllipseDetector::PrePeocessing(Mat1b& I,
										 Mat1b& DP,
										 Mat1b& DN
										 )
{

	Tic(0); //edge detection
	// Mid smooth
	//medianBlur(I,I,3);
	// Smooth image
	GaussianBlur(I, I, _szPreProcessingGaussKernelSize, _dPreProcessingGaussSigma);

	// Temp variables
	Mat1b E;				//edge mask
	Mat1s DX, DY;			//sobel derivatives

	// Detect edges
	Canny3(I, E, DX, DY, 3, false);

	Toc(0); //edge detection

	Tac(1); //preprocessing

	// For each edge points, compute the edge direction
	for (int i = 0; i<_szImg.height; ++i)
	{
		short* _dx = DX.ptr<short>(i);
		short* _dy = DY.ptr<short>(i);
		uchar* _e = E.ptr<uchar>(i);
		uchar* _dp = DP.ptr<uchar>(i);
		uchar* _dn = DN.ptr<uchar>(i);

		for (int j = 0; j<_szImg.width; ++j)
		{
			if (!((_e[j] <= 0) || (_dx[j] == 0) || (_dy[j] == 0)))
			{
				// Angle of the tangent
				float phi = -(float(_dx[j]) / float(_dy[j]));

				// Along positive or negative diagonal
				if (phi > 0)	_dp[j] = (uchar)255;
				else if (phi < 0)	_dn[j] = (uchar)255;
			}
		}
	}
	Toc(1); //preprocessing
};


void CNEllipseDetector::DetectAfterPreProcessing(vector<Ellipse>& ellipses, Mat1b& E, Mat1f& PHI)
{
	// Set the image size
	_szImg = E.size();

	// Initialize temporary data structures
	Mat1b DP = Mat1b::zeros(_szImg);		// arcs along positive diagonal
	Mat1b DN = Mat1b::zeros(_szImg);		// arcs along negative diagonal

	// For each edge points, compute the edge direction
	for (int i = 0; i<_szImg.height; ++i)
	{
		float* _phi = PHI.ptr<float>(i);
		uchar* _e = E.ptr<uchar>(i);
		uchar* _dp = DP.ptr<uchar>(i);
		uchar* _dn = DN.ptr<uchar>(i);

		for (int j = 0; j<_szImg.width; ++j)
		{
			if ((_e[j] > 0) && (_phi[j] != 0))
			{
				// Angle

				// along positive or negative diagonal
				if (_phi[j] > 0)	_dp[j] = (uchar)255;
				else if (_phi[j] < 0)	_dn[j] = (uchar)255;
			}
		}
	}

	// Initialize accumulator dimensions
	ACC_N_SIZE = 101;
	ACC_R_SIZE = 180;
	ACC_A_SIZE = max(_szImg.height, _szImg.width);

	// Allocate accumulators
	accN = new int[ACC_N_SIZE];
	accR = new int[ACC_R_SIZE];
	accA = new int[ACC_A_SIZE];

	// Other temporary 
	VVP points_1, points_2, points_3, points_4;		//vector of points, one for each convexity class
	unordered_map<uint, EllipseData> centers;		//hash map for reusing already computed EllipseData

	// Detect edges and find convexities
	DetectEdges13(DP, points_1, points_3);
	DetectEdges24(DN, points_2, points_4);

	// Find triplets
	Triplets124(points_1, points_2, points_4, centers, ellipses);
	Triplets231(points_2, points_3, points_1, centers, ellipses);
	Triplets342(points_3, points_4, points_2, centers, ellipses);
	Triplets413(points_4, points_1, points_3, centers, ellipses);

	// Sort detected ellipses with respect to score
	sort(ellipses.begin(), ellipses.end());

	//free accumulator memory
	delete[] accN;
	delete[] accR;
	delete[] accA;

	//cluster detections
	//ClusterEllipses(ellipses);
};


void CNEllipseDetector::Detect_2(Mat1b& I, vector<Ellipse>& ellipses)
{
	countsOfFindEllipse=0;
	countsOfGetFastCenter=0;
	Tic(1); //prepare data structure

	// Set the image size
	_szImg = I.size();

	// Initialize temporary data structures
	Mat1b DP = Mat1b::zeros(_szImg);		// arcs along positive diagonal
	Mat1b DN = Mat1b::zeros(_szImg);		// arcs along negative diagonal

	// Initialize accumulator dimensions
	ACC_N_SIZE = 101;
	ACC_R_SIZE = 180;
	ACC_A_SIZE = max(_szImg.height, _szImg.width);

	// Allocate accumulators
	accN = new int[ACC_N_SIZE];
	accR = new int[ACC_R_SIZE];
	accA = new int[ACC_A_SIZE];

	// Other temporary 
	VVP points_1, points_2, points_3, points_4;		//vector of points, one for each convexity class
	unordered_map<uint, EllipseData> centers;		//hash map for reusing already computed EllipseData

	Toc(1); //prepare data structure

	// Preprocessing
	// From input image I, find edge point with coarse convexity along positive (DP) or negative (DN) diagonal
	PrePeocessing(I, DP, DN);

	Tac(1); //preprocessing
	// Detect edges and find convexities
	DetectEdges13(DP, points_1, points_3);
	DetectEdges24(DN, points_2, points_4);

	Toc(1); //preprocessing


	//// DEBUG
	//Mat3b out(I.rows, I.cols, Vec3b(0,0,0));
	//for(unsigned i=0; i<points_1.size(); ++i)
	//{
	//	//Vec3b color(rand()%255, 128+rand()%127, 128+rand()%127);
	//	Vec3b color(255,0,0);
	//	for(unsigned j=0; j<points_1[i].size(); ++j)
	//		out(points_1[i][j]) = color;
	//}

	//for(unsigned i=0; i<points_2.size(); ++i)
	//{
	//	//Vec3b color(rand()%255, 128+rand()%127, 128+rand()%127);
	//	Vec3b color(0,255,0);
	//	for(unsigned j=0; j<points_2[i].size(); ++j)
	//		out(points_2[i][j]) = color;
	//}
	//for(unsigned i=0; i<points_3.size(); ++i)
	//{
	//	//Vec3b color(rand()%255, 128+rand()%127, 128+rand()%127);
	//	Vec3b color(0,0,255);
	//	for(unsigned j=0; j<points_3[i].size(); ++j)
	//		out(points_3[i][j]) = color;
	//}

	//for(unsigned i=0; i<points_4.size(); ++i)
	//{
	//	//Vec3b color(rand()%255, 128+rand()%127, 128+rand()%127);
	//	Vec3b color(255,0,255);
	//	for(unsigned j=0; j<points_4[i].size(); ++j)
	//		out(points_4[i][j]) = color;
	//}


	// time estimation, validation  inside

	Tic(2); //grouping
	//find triplets
	Triplets124(points_1, points_2, points_4, centers, ellipses);
	Triplets231(points_2, points_3, points_1, centers, ellipses);
	Triplets342(points_3, points_4, points_2, centers, ellipses);
	Triplets413(points_4, points_1, points_3, centers, ellipses);
	Toc(2); //grouping	
	// time estimation, validation inside
	_times[2] -= (_times[3] + _times[4]);//不够准确

	Tac(4); //validation
	// Sort detected ellipses with respect to score
	sort(ellipses.begin(), ellipses.end());//按照分数排序 分数相同 越接近圆越小
	Toc(4); //validation


	// Free accumulator memory
	delete[] accN;
	delete[] accR;
	delete[] accA;

	Tic(5);
	// Cluster detections
	ClusterEllipses(ellipses);
	Toc(5);

};
void CNEllipseDetector::Detect_rect(Mat1b& I, vector<Ellipse>& ellipses, Rectangle rect)
{
	countsOfFindEllipse = 0;
	countsOfGetFastCenter = 0;
	Tic(1); //prepare data structure

	// Set the image size
	_szImg = I.size();

	// Initialize temporary data structures
	Mat1b DP = Mat1b::zeros(_szImg);		// arcs along positive diagonal
	Mat1b DN = Mat1b::zeros(_szImg);		// arcs along negative diagonal

	// Initialize accumulator dimensions
	ACC_N_SIZE = 101;
	ACC_R_SIZE = 180;
	ACC_A_SIZE = max(_szImg.height, _szImg.width);

	// Allocate accumulators
	accN = new int[ACC_N_SIZE];
	accR = new int[ACC_R_SIZE];
	accA = new int[ACC_A_SIZE];

	// Other temporary 
	VVP points_1, points_2, points_3, points_4;		//vector of points, one for each convexity class
	unordered_map<uint, EllipseData> centers;		//hash map for reusing already computed EllipseData

	Toc(1); //prepare data structure

	// Preprocessing
	// From input image I, find edge point with coarse convexity along positive (DP) or negative (DN) diagonal
	PrePeocessing(I, DP, DN);

	namedWindow("PrePeocessing", CV_WINDOW_AUTOSIZE);
	imshow("PrePeocessing", DN + DP);

	Tac(1); //preprocessing
	// Detect edges and find convexities
	DetectEdges13_zrk(DP, points_1, points_3);
	DetectEdges24_zrk(DN, points_2, points_4);

	namedWindow("detectEdges_13", CV_WINDOW_AUTOSIZE);
	imshow("detectEdges_13", DP);

	namedWindow("detectEdges_24", CV_WINDOW_AUTOSIZE);
	imshow("detectEdges_24", DN);

	Toc(1); //preprocessing


	//// DEBUG
	//Mat3b out(I.rows, I.cols, Vec3b(0,0,0));
	//for(unsigned i=0; i<points_1.size(); ++i)
	//{
	//	//Vec3b color(rand()%255, 128+rand()%127, 128+rand()%127);
	//	Vec3b color(255,0,0);
	//	for(unsigned j=0; j<points_1[i].size(); ++j)
	//		out(points_1[i][j]) = color;
	//}

	//for(unsigned i=0; i<points_2.size(); ++i)
	//{
	//	//Vec3b color(rand()%255, 128+rand()%127, 128+rand()%127);
	//	Vec3b color(0,255,0);
	//	for(unsigned j=0; j<points_2[i].size(); ++j)
	//		out(points_2[i][j]) = color;
	//}
	//for(unsigned i=0; i<points_3.size(); ++i)
	//{
	//	//Vec3b color(rand()%255, 128+rand()%127, 128+rand()%127);
	//	Vec3b color(0,0,255);
	//	for(unsigned j=0; j<points_3[i].size(); ++j)
	//		out(points_3[i][j]) = color;
	//}

	//for(unsigned i=0; i<points_4.size(); ++i)
	//{
	//	//Vec3b color(rand()%255, 128+rand()%127, 128+rand()%127);
	//	Vec3b color(255,0,255);
	//	for(unsigned j=0; j<points_4[i].size(); ++j)
	//		out(points_4[i][j]) = color;
	//}


	// time estimation, validation  inside

	Tic(2); //grouping
	//find triplets
	/*
	for (int i = 0; i < points_1.size(); i++)
	{
		VP temp = points_1[i];
		cout << "第" << i << "条第1象限的弧" << endl;
		for (int j = 0; j < temp.size(); j++)
		{
			cout << "(" << temp[j].x + rect.x << "," << temp[j].y + rect.y << ")" << " ";
		}
		cout << endl;

	}
	cout << "----------------" << endl;
	for (int i = 0; i < points_2.size(); i++)
	{
		VP temp = points_2[i];
		cout << "第" << i << "条第2象限的弧" << endl;
		for (int j = 0; j < temp.size(); j++)
		{
			cout << "(" << temp[j].x + rect.x << "," << temp[j].y + rect.y << ")" << " ";
		}
		cout << endl;

	}
	cout << "----------------" << endl;
	for (int i = 0; i < points_3.size(); i++)
	{
		VP temp = points_3[i];
		cout << "第" << i << "条第3象限的弧" << endl;
		for (int j = 0; j < temp.size(); j++)
		{
			cout << "(" << temp[j].x + rect.x << "," << temp[j].y+rect.y << ")" << " ";
		}
		cout << endl;

	}
	cout << "----------------" << endl;

	for (int i = 0; i < points_4.size(); i++)
	{
		VP temp = points_4[i];
		cout << "第" << i << "条第4象限的弧" << endl;
		for (int j = 0; j < temp.size(); j++)
		{
			cout << "(" << temp[j].x + rect.x << "," << temp[j].y + rect.y << ")" << " ";
		}
		cout << endl;

	}
	cout << "----------------" << endl;
	//张睿珂在瞎搞
	//points_1[1].pop_back();
	//VP temp = points_2[1];
	//temp.insert(temp.begin(), Point(345 - rect.x ,111 - rect.y));
	//temp.insert(temp.begin(), Point(344 - rect.x, 112 - rect.y));
	//temp.insert(temp.begin(), Point(344 - rect.x, 113 - rect.y));
	//temp.insert(temp.begin(), Point(343 - rect.x, 113 - rect.y));
	//temp.insert(temp.begin(), Point(343 - rect.x, 114 - rect.y));
	//temp.insert(temp.begin(), Point(343 - rect.x, 115 - rect.y));

	for (int i = 0; i < points_2.size(); i++)
	{
		VP temp = points_2[i];
		cout << "第" << i << "条第2象限的弧" << endl;
		for (int j = 0; j < temp.size(); j++)
		{
			cout << "(" << temp[j].x + rect.x << "," << temp[j].y + rect.y << ")" << " ";
		}
		cout << endl;

	}*/
	cout << "1st " << points_1.size() << endl;
	cout << "2nd " << points_2.size() << endl;
	cout << "3rd " << points_3.size() << endl;
	cout << "4th " << points_4.size() << endl;
	Triplets124_rect(points_1, points_2, points_4, centers, ellipses, rect);
	Triplets231_rect(points_2, points_3, points_1, centers, ellipses, rect);
	Triplets342_rect(points_3, points_4, points_2, centers, ellipses, rect);
	Triplets413_rect(points_4, points_1, points_3, centers, ellipses, rect);

	Toc(2); //grouping	
	// time estimation, validation inside
	_times[2] -= (_times[3] + _times[4]);//不够准确

	Tac(4); //validation
	// Sort detected ellipses with respect to score
	sort(ellipses.begin(), ellipses.end());//按照分数排序 分数相同 越接近圆越小
	//cout << "after sorting there are " << ellipses.size()<< "ellipses" << endl;
	Toc(4); //validation


	// Free accumulator memory
	delete[] accN;
	delete[] accR;
	delete[] accA;

	Tic(5);
	// Cluster detections
	ClusterEllipses(ellipses);
	Toc(5);
	//cout << "after clustering there are " << ellipses.size() << "ellipses" << endl;
};



void CNEllipseDetector::Detect_3(Mat1b& I, vector<Ellipse>& ellipses)
{
	countsOfFindEllipse = 0;
	countsOfGetFastCenter = 0;
	Tic(1); //prepare data structure

	// Set the image size
	_szImg = I.size();

	// Initialize temporary data structures
	Mat1b DP = Mat1b::zeros(_szImg);		// arcs along positive diagonal
	Mat1b DN = Mat1b::zeros(_szImg);		// arcs along negative diagonal

	// Initialize accumulator dimensions
	ACC_N_SIZE = 101;
	ACC_R_SIZE = 180;
	ACC_A_SIZE = max(_szImg.height, _szImg.width);

	// Allocate accumulators
	accN = new int[ACC_N_SIZE];
	accR = new int[ACC_R_SIZE];
	accA = new int[ACC_A_SIZE];

	// Other temporary 
	VVP points_1, points_2, points_3, points_4;		//vector of points, one for each convexity class
	unordered_map<uint, EllipseData> centers;		//hash map for reusing already computed EllipseData

	Toc(1); //prepare data structure

	// Preprocessing
	// From input image I, find edge point with coarse convexity along positive (DP) or negative (DN) diagonal
	PrePeocessing(I, DP, DN);

	namedWindow("PrePeocessing", CV_WINDOW_AUTOSIZE);
	imshow("PrePeocessing", DN + DP);

	Tac(1); //preprocessing
	// Detect edges and find convexities
	DetectEdges13(DP, points_1, points_3);
	DetectEdges24(DN, points_2, points_4);

	namedWindow("detectEdges_13", CV_WINDOW_AUTOSIZE);
	imshow("detectEdges_13", DP);

	namedWindow("detectEdges_24", CV_WINDOW_AUTOSIZE);
	imshow("detectEdges_24", DN);
	Toc(1); //preprocessing


	//// DEBUG
	//Mat3b out(I.rows, I.cols, Vec3b(0,0,0));
	//for(unsigned i=0; i<points_1.size(); ++i)
	//{
	//	//Vec3b color(rand()%255, 128+rand()%127, 128+rand()%127);
	//	Vec3b color(255,0,0);
	//	for(unsigned j=0; j<points_1[i].size(); ++j)
	//		out(points_1[i][j]) = color;
	//}

	//for(unsigned i=0; i<points_2.size(); ++i)
	//{
	//	//Vec3b color(rand()%255, 128+rand()%127, 128+rand()%127);
	//	Vec3b color(0,255,0);
	//	for(unsigned j=0; j<points_2[i].size(); ++j)
	//		out(points_2[i][j]) = color;
	//}
	//for(unsigned i=0; i<points_3.size(); ++i)
	//{
	//	//Vec3b color(rand()%255, 128+rand()%127, 128+rand()%127);
	//	Vec3b color(0,0,255);
	//	for(unsigned j=0; j<points_3[i].size(); ++j)
	//		out(points_3[i][j]) = color;
	//}

	//for(unsigned i=0; i<points_4.size(); ++i)
	//{
	//	//Vec3b color(rand()%255, 128+rand()%127, 128+rand()%127);
	//	Vec3b color(255,0,255);
	//	for(unsigned j=0; j<points_4[i].size(); ++j)
	//		out(points_4[i][j]) = color;
	//}


	// time estimation, validation  inside
	/*
	for (int i = 0; i < points_1.size(); i++)
	{
		VP temp = points_1[i];
		cout << "第" << i << "条第1象限的弧" << endl;
		for (int j = 0; j < temp.size(); j++)
		{
			cout << "(" << temp[j].x  << "," << temp[j].y  << ")" << " ";
		}
		cout << endl;

	}
	cout << "----------------" << endl;
	for (int i = 0; i < points_2.size(); i++)
	{
		VP temp = points_2[i];
		cout << "第" << i << "条第2象限的弧" << endl;
		for (int j = 0; j < temp.size(); j++)
		{
			cout << "(" << temp[j].x  << "," << temp[j].y  << ")" << " ";
		}
		cout << endl;

	}
	cout << "----------------" << endl;
	for (int i = 0; i < points_3.size(); i++)
	{
		VP temp = points_3[i];
		cout << "第" << i << "条第3象限的弧" << endl;
		for (int j = 0; j < temp.size(); j++)
		{
			cout << "(" << temp[j].x  << "," << temp[j].y  << ")" << " ";
		}
		cout << endl;

	}
	cout << "----------------" << endl;

	for (int i = 0; i < points_4.size(); i++)
	{
		VP temp = points_4[i];
		cout << "第" << i << "条第4象限的弧" << endl;
		for (int j = 0; j < temp.size(); j++)
		{
			cout << "(" << temp[j].x << "," << temp[j].y << ")" << " ";
		}
		cout << endl;

	}
	cout << "----------------" << endl;*/
	Tic(2); //grouping
	//find triplets
	Triplets124_3(points_1, points_2, points_4, centers, ellipses);
	Triplets231_3(points_2, points_3, points_1, centers, ellipses);
	Triplets342_3(points_3, points_4, points_2, centers, ellipses);
	Triplets413_3(points_4, points_1, points_3, centers, ellipses);
	Toc(2); //grouping	
	// time estimation, validation inside
	_times[2] -= (_times[3] + _times[4]);//不够准确

	Tac(4); //validation
	// Sort detected ellipses with respect to score
	sort(ellipses.begin(), ellipses.end());//按照分数排序 分数相同 越接近圆越小
	Toc(4); //validation


	// Free accumulator memory
	delete[] accN;
	delete[] accR;
	delete[] accA;

	Tic(5);
	// Cluster detections
	ClusterEllipses(ellipses);
	Toc(5);

};

void CNEllipseDetector::Detect_ArcsIn(Mat1b& I, vector<Ellipse>& ellipses)
{
	countsOfFindEllipse = 0;
	countsOfGetFastCenter = 0;
	Tic(1); //prepare data structure

	// Set the image size
	_szImg = I.size();

	// Initialize temporary data structures
	Mat1b DP = Mat1b::zeros(_szImg);		// arcs along positive diagonal
	Mat1b DN = Mat1b::zeros(_szImg);		// arcs along negative diagonal

	// Initialize accumulator dimensions
	ACC_N_SIZE = 101;
	ACC_R_SIZE = 180;
	ACC_A_SIZE = max(_szImg.height, _szImg.width);

	// Allocate accumulators
	accN = new int[ACC_N_SIZE];
	accR = new int[ACC_R_SIZE];
	accA = new int[ACC_A_SIZE];

	// Other temporary 
	VVP points_1, points_2, points_3, points_4;		//vector of points, one for each convexity class
	unordered_map<uint, EllipseData> centers;		//hash map for reusing already computed EllipseData

	Toc(1); //prepare data structure

	// Preprocessing
	// From input image I, find edge point with coarse convexity along positive (DP) or negative (DN) diagonal
	PrePeocessing(I, DP, DN);

	namedWindow("PrePeocessing", CV_WINDOW_AUTOSIZE);
	imshow("PrePeocessing", DN + DP);

	Tac(1); //preprocessing
	// Detect edges and find convexities
	DetectEdges13(DP, points_1, points_3);
	DetectEdges24(DN, points_2, points_4);

	namedWindow("detectEdges_13", CV_WINDOW_AUTOSIZE);
	imshow("detectEdges_13", DP);

	namedWindow("detectEdges_24", CV_WINDOW_AUTOSIZE);
	imshow("detectEdges_24", DN);
	Toc(1); //preprocessing


	//// DEBUG
	//Mat3b out(I.rows, I.cols, Vec3b(0,0,0));
	//for(unsigned i=0; i<points_1.size(); ++i)
	//{
	//	//Vec3b color(rand()%255, 128+rand()%127, 128+rand()%127);
	//	Vec3b color(255,0,0);
	//	for(unsigned j=0; j<points_1[i].size(); ++j)
	//		out(points_1[i][j]) = color;
	//}

	//for(unsigned i=0; i<points_2.size(); ++i)
	//{
	//	//Vec3b color(rand()%255, 128+rand()%127, 128+rand()%127);
	//	Vec3b color(0,255,0);
	//	for(unsigned j=0; j<points_2[i].size(); ++j)
	//		out(points_2[i][j]) = color;
	//}
	//for(unsigned i=0; i<points_3.size(); ++i)
	//{
	//	//Vec3b color(rand()%255, 128+rand()%127, 128+rand()%127);
	//	Vec3b color(0,0,255);
	//	for(unsigned j=0; j<points_3[i].size(); ++j)
	//		out(points_3[i][j]) = color;
	//}

	//for(unsigned i=0; i<points_4.size(); ++i)
	//{
	//	//Vec3b color(rand()%255, 128+rand()%127, 128+rand()%127);
	//	Vec3b color(255,0,255);
	//	for(unsigned j=0; j<points_4[i].size(); ++j)
	//		out(points_4[i][j]) = color;
	//}


	// time estimation, validation  inside
	/*
	for (int i = 0; i < points_1.size(); i++)
	{
		VP temp = points_1[i];
		cout << "第" << i << "条第1象限的弧" << endl;
		for (int j = 0; j < temp.size(); j++)
		{
			cout << "(" << temp[j].x  << "," << temp[j].y  << ")" << " ";
		}
		cout << endl;

	}
	cout << "----------------" << endl;
	for (int i = 0; i < points_2.size(); i++)
	{
		VP temp = points_2[i];
		cout << "第" << i << "条第2象限的弧" << endl;
		for (int j = 0; j < temp.size(); j++)
		{
			cout << "(" << temp[j].x  << "," << temp[j].y  << ")" << " ";
		}
		cout << endl;

	}
	cout << "----------------" << endl;
	for (int i = 0; i < points_3.size(); i++)
	{
		VP temp = points_3[i];
		cout << "第" << i << "条第3象限的弧" << endl;
		for (int j = 0; j < temp.size(); j++)
		{
			cout << "(" << temp[j].x  << "," << temp[j].y  << ")" << " ";
		}
		cout << endl;

	}
	cout << "----------------" << endl;

	for (int i = 0; i < points_4.size(); i++)
	{
		VP temp = points_4[i];
		cout << "第" << i << "条第4象限的弧" << endl;
		for (int j = 0; j < temp.size(); j++)
		{
			cout << "(" << temp[j].x << "," << temp[j].y << ")" << " ";
		}
		cout << endl;

	}
	cout << "----------------" << endl;*/
	Tic(2); //grouping
	//find triplets
	Triplets124_ArcsIn(points_1, points_2, points_4, centers, ellipses);
	Triplets231_ArcsIn(points_2, points_3, points_1, centers, ellipses);
	Triplets342_ArcsIn(points_3, points_4, points_2, centers, ellipses);
	Triplets413_ArcsIn(points_4, points_1, points_3, centers, ellipses);
	Toc(2); //grouping	
	// time estimation, validation inside
	_times[2] -= (_times[3] + _times[4]);//不够准确

	Tac(4); //validation
	// Sort detected ellipses with respect to score
	sort(ellipses.begin(), ellipses.end());//按照分数排序 分数相同 越接近圆越小
	Toc(4); //validation


	// Free accumulator memory
	delete[] accN;
	delete[] accR;
	delete[] accA;

	Tic(5);
	// Cluster detections
	ClusterEllipses(ellipses);
	Toc(5);

};



// Ellipse clustering procedure. See Sect [3.3.2] in the paper.
void CNEllipseDetector::ClusterEllipses(vector<Ellipse>& ellipses)
{
	float th_Da = 0.1f;
	float th_Db = 0.1f;
	float th_Dr = 0.1f;

	float th_Dc_ratio = 0.1f;
	float th_Dr_circle = 0.9f;

	int iNumOfEllipses = int(ellipses.size());
	if (iNumOfEllipses == 0) return;

	// The first ellipse is assigned to a cluster
	vector<Ellipse> clusters;
	clusters.push_back(ellipses[0]);

	bool bFoundCluster = false;

	for (int i = 1; i<iNumOfEllipses; ++i)
	{
		Ellipse& e1 = ellipses[i];

		int sz_clusters = int(clusters.size());

		float ba_e1 = e1._b / e1._a;
		float Decc1 = e1._b / e1._a;

		bool bFoundCluster = false;
		for (int j = 0; j<sz_clusters; ++j)
		{
			Ellipse& e2 = clusters[j];

			float ba_e2 = e2._b / e2._a;
			float th_Dc = min(e1._b, e2._b) * th_Dc_ratio;
			th_Dc *= th_Dc;

			// Centers
			float Dc = ((e1._xc - e2._xc)*(e1._xc - e2._xc) + (e1._yc - e2._yc)*(e1._yc - e2._yc));
			if (Dc > th_Dc)
			{
				//not same cluster
				continue;
			}

			// a
			float Da = abs(e1._a - e2._a) / max(e1._a, e2._a);
			if (Da > th_Da)
			{
				//not same cluster
				continue;
			}

			// b
			float Db = abs(e1._b - e2._b) / min(e1._b, e2._b);
			if (Db > th_Db)
			{
				//not same cluster
				continue;
			}

			// angle
			float Dr = GetMinAnglePI(e1._rad, e2._rad) / float(CV_PI);
			if ((Dr > th_Dr) && (ba_e1 < th_Dr_circle) && (ba_e2 < th_Dr_circle))
			{
				//not same cluster
				continue;
			}

			// Same cluster as e2
			bFoundCluster = true;//
			// Discard, no need to create a new cluster
			break;
		}

		if (!bFoundCluster)
		{
			// Create a new cluster			
			clusters.push_back(e1);
		}
	}

	clusters.swap(ellipses);
};



//Draw at most iTopN detected ellipses.
void CNEllipseDetector::DrawDetectedEllipses(Mat3b& output, vector<Ellipse>& ellipses, int iTopN, int thickness)
{
	int sz_ell = int(ellipses.size());
	int n = (iTopN == 0) ? sz_ell : min(iTopN, sz_ell);
	for (int i = 0; i < n; ++i)
	{
		Ellipse& e = ellipses[n - i - 1];
		int g = cvRound(e._score * 255.f);
		Scalar color(0, g, 0);
		ellipse(output, Point(cvRound(e._xc), cvRound(e._yc)), Size(cvRound(e._a), cvRound(e._b)), e._rad*180.0 / CV_PI, 0.0, 360.0, color, thickness);
	}
}

// add 2015年12月24日
void CNEllipseDetector:: showEdgeInPic(Mat1b& I){
	_szImg = I.size();
	// initialize temporary data structures
	Mat1b DP	= Mat1b::zeros(_szImg);		// arcs along positive diagonal
	Mat1b DN	= Mat1b::zeros(_szImg);		// arcs along negative diagonal

	//other temporary 
	vector<vector<Point>> points_1, points_2, points_3, points_4;		//vector of points, one for each convexity class

	//Preprocessing
	//From input image I, find edge point with coarse convexity along positive (DP) or negative (DN) diagonal


	PrePeocessing(I, DP, DN);

	// show DP
	Mat1b DTMP	= Mat1b::ones(_szImg)*255;	
	Mat1b tmp=DTMP-DP;
	Mat1b tmp2=DTMP-DP;
	/*cvShowImage("DP",&(IplImage)tmp );
	cvSaveImage("DP.jpg",&IplImage(tmp));
	tmp=DTMP-DN;
	cvShowImage("DN",&(IplImage)tmp );
	cvSaveImage("DN.jpg",&IplImage(tmp));*/
	tmp = DP+DN;
	tmp2 = DTMP-tmp;
	cvShowImage("DNP",&(IplImage)tmp2 );
	cvSaveImage("DPN.jpg",&IplImage(tmp2));
	//end show DpDn
	float t_fMinOrientedRectSide = _fMinOrientedRectSide;
	_fMinOrientedRectSide=4;
	//detect edges and find convexities// 4个象限
	DetectEdges13(DP, points_1, points_3);
	DetectEdges24(DN, points_2, points_4);
	_fMinOrientedRectSide=t_fMinOrientedRectSide;
	//显示选出的边
	Mat picture(_szImg, CV_8UC3,Scalar(255,255,255));
	IplImage *image= cvCreateImage(_szImg,IPL_DEPTH_8U,3);
	// Mat picture(_szImg,CV_8UC3);
	/*for(int ih=0;ih<_szImg.height;ih++){
	for(int iw=0;iw<_szImg.width;iw++){
	picture.at<Vec3b>(ih,iw)[0]=I.at<uchar>(ih,iw);
	picture.at<Vec3b>(ih,iw)[1]=I.at<uchar>(ih,iw);
	picture.at<Vec3b>(ih,iw)[2]=I.at<uchar>(ih,iw);

	}
	}*/
	vector<Mat> imgs(4);
	Mat picture1=picture.clone();
	Mat picture2=picture.clone();
	Mat picture3=picture.clone();
	Mat picture4=picture.clone();
	showEdge(points_1,picture1);
	showEdge(points_2,picture2);
	showEdge(points_3,picture3);
	showEdge(points_4,picture4);
	/*
	imshow("第1象限", picture1);
	imshow("第2象限", picture2);
	imshow("第3象限", picture3);
	imshow("第4象限", picture4);*/
	
	Mat picture5=picture.clone();
	showEdge(points_1,picture5);
	showEdge(points_2,picture5);
	showEdge(points_3,picture5);
	showEdge(points_4,picture5);

	imgs[0] = picture2;
	imgs[1] = picture1;
	imgs[2] = picture3;
	imgs[3] = picture4; 
	
	cvNamedWindow("all arcs", CV_WINDOW_NORMAL);  
	cvShowImage("all arcs",&(IplImage)picture5 );
	cvSaveImage("arcs.jpg",&IplImage(picture5));

	//imshow("all arcs", picture5);

	MultiImage_OneWin("象限中的弧", imgs, cvSize(2, 2), cvSize(400,800)); 
	//cvWaitKey(0);  
	//cvDestroyWindow("all arcs");
}
void CNEllipseDetector:: showAllEdgeInPic(Mat1b& I){
	_szImg = I.size();
	// initialize temporary data structures
	Mat1b DP	= Mat1b::zeros(_szImg);		// arcs along positive diagonal
	Mat1b DN	= Mat1b::zeros(_szImg);		// arcs along negative diagonal

	//other temporary 
	vector<vector<Point>> points_1, points_2, points_3, points_4;		//vector of points, one for each convexity class

	//Preprocessing
	//From input image I, find edge point with coarse convexity along positive (DP) or negative (DN) diagonal

	//this->_fMaxRectAxesRatio=0;
	//this->_iMinEdgeLength=0;
	//this->_uNs=0;

	PrePeocessing(I, DP, DN);
	cvNamedWindow("aaas", CV_WINDOW_NORMAL);
	Mat1b tmp	= Mat1b::zeros(_szImg);
	

	cvShowImage("aaas",&(IplImage)DP );
	//detect edges and find convexities// 4个象限
	//DetectEdges13(DP, points_1, points_3);
	//DetectEdges24(DN, points_2, points_4);
}

int CNEllipseDetector:: showEdgeInPic(Mat1b& I,bool showedge){
	_szImg = I.size();
	// initialize temporary data structures
	Mat1b DP	= Mat1b::zeros(_szImg);		// arcs along positive diagonal
	Mat1b DN	= Mat1b::zeros(_szImg);		// arcs along negative diagonal

	//other temporary 
	vector<vector<Point>> points_1, points_2, points_3, points_4;		//vector of points, one for each convexity class

	//Preprocessing
	//From input image I, find edge point with coarse convexity along positive (DP) or negative (DN) diagonal


	PrePeocessing(I, DP, DN);
	//detect edges and find convexities// 4个象限
	DetectEdges13(DP, points_1, points_3);
	DetectEdges24(DN, points_2, points_4);
	int EdgesNumber=points_1.size()+points_2.size()+points_3.size()+points_4.size();
	//显示选出的边
	return EdgesNumber;
	//cvWaitKey(0);  
	//cvDestroyWindow("all arcs");
}