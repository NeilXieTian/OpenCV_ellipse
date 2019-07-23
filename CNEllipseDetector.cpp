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


//zh
#include <opencv2/opencv.hpp>


//
////zh���ӵ�˫������ֵ
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
using namespace cv;
using namespace std;
#define PI 3.14159265
float BiCubicPoly(float x);
void MyScaleBiCubicInter(Mat& src, Mat& dst, float TransMat[3][3]);




void getEllipsePara(RotatedRect & ellipsemege, EllipsePara& EP_t)
{


	EP_t.theta = ellipsemege.angle * CV_PI / 180.0;
	EP_t.a = max(ellipsemege.size.width, ellipsemege.size.height) / 2.0;
	EP_t.b = min(ellipsemege.size.width, ellipsemege.size.height) / 2.0;


	EP_t.c.x = ellipsemege.center.x;
	EP_t.c.y = ellipsemege.center.y;

	/*//����Բһ�㷽���еĲ���
	EP_t.A = a * a * sin(theta) * sin(theta) + b * b * cos(theta) * cos(theta);
	EP_t.B = (-2.0) * (a * a - b * b) * sin(theta) * cos(theta);
	EP_t.C = a * a * cos(theta) * cos(theta) + b * b * sin(theta) * sin(theta);
	EP_t.F = (-1.0) * a * a * b * b;
	*/

	cout << "theta: " << EP_t.theta << " x: " << EP_t.c.x << " y: " << EP_t.c.y << " a: " << EP_t.a << " b: " << EP_t.b << endl;
	cout << "�Ƕȣ�" << ellipsemege.angle << endl;

}

////12.12Xt������С���˷������Բ
Point2f LeastSquareFittingEllipse(vector<Point> temp_coordinates)//QR�ֽ�
{
	float x1 = 0;
	float x2 = 0;
	float x3 = 0;
	float x4 = 0;
	float y1 = 0;
	float y2 = 0;
	float y3 = 0;
	float y4 = 0;
	float x1y1 = 0;
	float x1y2 = 0;
	float x1y3 = 0;
	float x2y1 = 0;
	float x2y2 = 0;
	float x3y1 = 0;
	int num;
	vector<Point>::iterator k;
	Point tempcircle;
	num = temp_coordinates.size();
	for (k = temp_coordinates.begin(); k != temp_coordinates.end(); k++)
	{
		x1 = x1 + (*k).x;
		x2 = x2 + pow((*k).x, 2);
		x3 = x3 + pow((*k).x, 3);
		x4 = x4 + pow((*k).x, 4);
		y1 = y1 + (*k).y;
		y2 = y2 + pow((*k).y, 2);
		y3 = y3 + pow((*k).y, 3);
		y4 = y4 + pow((*k).y, 4);
		x1y1 = x1y1 + (*k).x * (*k).y;
		x1y2 = x1y2 + (*k).x * pow((*k).y, 2);
		x1y3 = x1y3 + (*k).x * pow((*k).y, 3);
		x2y1 = x2y1 + pow((*k).x, 2) * (*k).y;
		x2y2 = x2y2 + pow((*k).x, 2) * pow((*k).y, 2);
		x3y1 = x3y1 + pow((*k).x, 3) * (*k).y;
	}
	Mat left_matrix = (Mat_<float>(6, 5) << x3y1, x2y2 - x4, x3, x2y1, x2,
		x2y2, x1y3 - x3y1, x2y1, x1y2, x1y1,
		x1y3, y4 - x2y2, x1y2, y3, y2,
		x2y1, x1y2 - x3, x2, x1y1, x1,
		x1y2, y3 - x2y1, x1y1, y2, y1,
		x1y1, y2 - x2, x1, y1, num);
	//cout << "left_matrix[6,5]=" << left_matrix << endl;
	Mat right_matrix = (Mat_<float>(6, 1) << -x4, -x3y1, -x2y2, -x3, -x2y1, -x2);
	//cout << "right_matrix[6,1]=" << right_matrix << endl;
	Mat ellipse_solution(5, 1, CV_32F);
	double t = getTickCount();
	solve(left_matrix, right_matrix, ellipse_solution, DECOMP_QR);
	t = getTickCount() - t;
	double ConsumeTime = t / (getTickFrequency());
	//cout << "time =" << ConsumeTime << endl;
	//cout << "QR_ellipse_solution[5,1]=" << ellipse_solution << endl;
	float A, B, C, D, E, F;
	A = 1 - ellipse_solution.at<float>(1);
	B = ellipse_solution.at<float>(0);
	C = ellipse_solution.at<float>(1);
	D = ellipse_solution.at<float>(2);
	E = ellipse_solution.at<float>(3);
	F = ellipse_solution.at<float>(4);
	tempcircle.x = (B*E - 2 * C*D) / (4 * A*C - B*B);
	tempcircle.y = (B*D - 2 * A*E) / (4 * A*C - B*B);
	cout << "QR:" << tempcircle.x << "		" << tempcircle.y << endl;
	return tempcircle;
}
Point2f LeastSquareFittingEllipse_LU(VP temp_coordinates, double *FittingResult)//LU�ֽ�
{
	float x1 = 0;
	float x2 = 0;
	float x3 = 0;
	//float x4 = 0;
	float y1 = 0;
	float y2 = 0;
	float y3 = 0;
	float y4 = 0;
	float x1y1 = 0;
	float x1y2 = 0;
	float x1y3 = 0;
	float x2y1 = 0;
	float x2y2 = 0;
	float x3y1 = 0;
	int num;
	VP::iterator k;
	Point2f tempcircle;
	num = temp_coordinates.size();
	for (k = temp_coordinates.begin(); k != temp_coordinates.end(); k++)
	{
		x1 = x1 + (*k).x;
		x2 = x2 + pow((*k).x, 2);
		x3 = x3 + pow((*k).x, 3);
		//x4 = x4 + pow((*k).x, 4);
		y1 = y1 + (*k).y;
		y2 = y2 + pow((*k).y, 2);
		y3 = y3 + pow((*k).y, 3);
		y4 = y4 + pow((*k).y, 4);
		x1y1 = x1y1 + (*k).x * (*k).y;
		x1y2 = x1y2 + (*k).x * pow((*k).y, 2);
		x1y3 = x1y3 + (*k).x * pow((*k).y, 3);
		x2y1 = x2y1 + pow((*k).x, 2) * (*k).y;
		x2y2 = x2y2 + pow((*k).x, 2) * pow((*k).y, 2);
		x3y1 = x3y1 + pow((*k).x, 3) * (*k).y;
	}
	Mat left_matrix = (Mat_<float>(5, 5) << x2y2, x1y3, x2y1, x1y2, x1y1,
		x1y3, y4, x1y2, y3, y2,
		x2y1, x1y2, x2, x1y1, x1,
		x1y2, y3, x1y1, y2, y1,
		x1y1, y2, x1, y1, num);
	//cout << "left_matrix[5,5]=" << left_matrix << endl;
	Mat right_matrix = (Mat_<float>(5, 1) << -x3y1, -x2y2, -x3, -x2y1, -x2);
	//cout << "right_matrix[5,1]=" << right_matrix << endl;
	Mat ellipse_solution(5, 1, CV_32F);
	solve(left_matrix, right_matrix, ellipse_solution, DECOMP_LU);
	//cout << "ellipse_solution[5,1]=" << ellipse_solution << endl;
	float A, B, C, D, E, F;
	A = ellipse_solution.at<float>(0);
	B = ellipse_solution.at<float>(1);
	C = ellipse_solution.at<float>(2);
	D = ellipse_solution.at<float>(3);
	E = ellipse_solution.at<float>(4);

	tempcircle.x = (2 * B*C - A*D) / (A*A - 4 * B);
	tempcircle.y = (2 * D - A*D) / (A*A - 4 * B);
	double fenzi = 2 * (A*C*D - B*C*C - D*D + 4*B*E - A*A*E);
	double fenmu = (A*A - 4*B)*(B-sqrt(A*A + (1-B)*(1-B))+1);
	double fenmu2 = (A*A - 4*B)*(B + sqrt(A*A + (1 - B)*(1 - B))+1);
	double Long = sqrt(fabs(fenzi / fenmu));
	double Short = sqrt(fabs(fenzi / fenmu2));
	double theta = atan(sqrt((Long * Long - Short * Short * B) / (Long * Long * B - Short * Short))) * 180 / CV_PI;
	cout << "��Ϻ�Ľ�����£�" << endl;
	cout << "a = " << Long << endl;
	cout << "b = " << Short << endl;
	cout << "x = " << tempcircle.x << endl;
	cout << "y = " << tempcircle.y << endl;
	cout << "theta = " << theta << endl;
	
	FittingResult[0] = tempcircle.x;
	FittingResult[1] = tempcircle.y;
	FittingResult[2] = Long;
	FittingResult[3] = Short;
	FittingResult[4] = theta;

	return tempcircle;
}
/////////
//Xt 2.2�������δ�����7�������㡢mean-shift����

//1.������һ������б��
double Ptangent(VP arc, int index_p0) {
	double R = 3;
	int sz = arc.size();
	int p = index_p0;
	Point p0 = arc[index_p0];
	int i = p + 1;
	for (; i++; i < sz)		//Ѱ��p1
	{
		Point current = arc[i];
		if (abs(sqrt((current.x - p0.x)*(current.x - p0.x) + (current.y - p0.y) * (current.y - p0.y)) - R) < 1 / sqrt(2))
			break;
	}
	Point p1 = arc[i];
	if(p1.x == arc[i].x)
	cout << "+++++++++++����Ѱ�ҵ��ĵ��ǵ� "<< i << "����������(" << p1.x << "," << p1.y << ")" << endl;

	i = p - 1;			
	for (; i--; i > 0)		//Ѱ��p2
	{
		Point current = arc[i];
		if (abs(sqrt((current.x - p0.x)*(current.x - p0.x) + (current.y - p0.y) * (current.y - p0.y)) - R) < 1 / sqrt(2))
			break;
	}
	Point p2 = arc[i];
	if (p2.x == arc[i].x)
	cout << "-----------����Ѱ�ҵ��ĵ��ǵ� " << i << "����������(" << p2.x << "," << p2.y << ")" << endl;
	cout << "б�ʣ�" << (p1.y - p2.y) / ((double)p1.x - p2.x) << endl;
	return (p1.y - p2.y) / ((double)p1.x - p2.x);
}

//������ֱ�ߵĽ���
struct LINE
{
	double slope;
	Point pStart;
	LINE(double k, Point point){
		slope = k;
		pStart = point;
	}
};
Point getCrossPoint(LINE& line1, LINE& line2)
{
	Point pt;
	pt.x = (line1.slope * line1.pStart.x - line2.slope * line2.pStart.x + line2.pStart.y - line1.pStart.y) / ((double)line1.slope - line2.slope);
	pt.y = line1.slope * (pt.x - line1.pStart.x) + line1.pStart.y;
	return pt;
}

//mean-shift�㷨
Point MeanShift(VP PointCloud)
{
	Point center = PointCloud[0];
	Point move = Point(0, 0);
	int time = 0;
	Point center_temp;
	do
	{
		for (int i = 0; i < PointCloud.size(); ++i)
		{
			move.x = move.x + PointCloud[i].x - center.x;
			move.y = move.y + PointCloud[i].y - center.y;
			cout << i << ": (" << PointCloud[i].x << "," << PointCloud[i].y << ")" << endl;
			cout << "		(" << move.x << "," << move.y << ")" << endl;

		}
		cout << "��Ҫ�ƶ���" << move.x / 7 << "," << move.y / 7 << "��" << endl;
		center.x += move.x / 7;
		center.y += move.y / 7;
		cout << "��ǰ�����ĵ��Ƶ��ˣ�" << center.x << "  " << center.y << endl;
		center_temp = center;
		time++;
	} while (time != 1 && sqrt(center.x * center.x + center.y * center.y) - sqrt(center_temp.x * center_temp.x + center_temp.y * center_temp.y) < 5);
	return center;
}

//׼��7��������
VP FindPointCloud(VP arc1, VP arc2, double& N, double& K, double& rho)
{
	//���C12��ȫ���̣�����ֱ��L1��L2
	VP PointCloud;
	int sz_e1 = arc1.size();
	cout << "����ĵ�һ�λ���sz_e1 = " << sz_e1 << endl;

	Point p1 = arc1[sz_e1 / 5];
	Point p12 = arc1[ sz_e1 / 2];
	Point p2 = arc1[sz_e1 - sz_e1/5];
	double k12 = Ptangent(arc1, sz_e1 / 2);
	double k1 = Ptangent(arc1, sz_e1 / 5);
	double k2 = Ptangent(arc1,sz_e1 - sz_e1 / 5);
	LINE Line1(k1, p1), Line2(k2, p2), Line12(k12, p12);
	cout << "P1a�������뿴���" << k1 << "	" << k12 << endl;
	
	Point P1a = getCrossPoint(Line1, Line12);
	cout << "P1a: " << P1a.x << ", " << P1a.y << endl;
	Point P2a = getCrossPoint(Line12, Line2);
	cout << "P2a: " << P2a.x << ", " << P2a.y << endl;
	Point m1 = Point((p1.x + p12.x) / 2, (p1.y + p12.y) / 2);
	Point m2 = Point((p12.x + p2.x) / 2, (p12.y + p2.y) / 2);

	cout << "p1: " << p1.x << ", " << p1.y << endl;
	cout << "p12: " << p12.x << ", " << p12.y << endl;
	cout << "m1: " << m1.x << ", " << m1.y << endl;

	double kl1 = (m1.y - P1a.y) / ((double)m1.x - P1a.x);
	double kl2 = (m2.y - P2a.y) / ((double)m2.x - P2a.x);
	cout << "kl1  kl2 �ֱ�Ϊ��" << kl1 << kl2 << endl;
	LINE L1(kl1, m1), L2(kl2, m2);
	Point C12 = getCrossPoint(L1, L2);		//�õ���һ������C12
	cout << "*C12: " << C12.x << ", " << C12.y << endl;
	PointCloud.push_back(C12);

	//���C34��ȫ���̣�����ֱ��L3��L4
	int sz_e2 = arc2.size();
	Point p3 = arc2[4];
	Point p34 = arc2[sz_e2 / 2];
	Point p4 = arc2[sz_e2 - 6];
	double k34 = Ptangent(arc2, sz_e2 / 2);
	double k3 = Ptangent(arc2, 4);
	double k4 = Ptangent(arc2, sz_e2 - 6);
	LINE Line3(k3, p3), Line4(k4, p4), Line34(k34, p34);

	Point p3a = getCrossPoint(Line1, Line12);
	Point p4a = getCrossPoint(Line12, Line2);
	Point m3 = Point((p3.x + p34.x) / 2, (p3.y + p34.y) / 2);
	Point m4 = Point((p34.x + p4.x) / 2, (p34.y + p4.y) / 2);
	double kl3 = (m3.y - p3a.y) / ((double)m3.x - p3a.x);
	double kl4 = (m4.y - p4a.y) / ((double)m4.x - p4a.x);
	LINE L3(kl3, m3), L4(kl4, m4);
	Point C34 = getCrossPoint(L3, L4);		//�õ��ڶ�������C34
	PointCloud.push_back(C34);

	//��������C12��C34��ƽ��ֵ
	PointCloud.push_back(Point((C12.x + C34.x) / 2, (C12.y + C34.y) / 2));
	//L1��L3�Ľ���
	Point C13 = getCrossPoint(L1, L3);
	PointCloud.push_back(C13);

	//L1��L4�Ľ���
	Point C14 = getCrossPoint(L1, L4);
	PointCloud.push_back(C14);

	//L2��L3�Ľ���
	Point C23 = getCrossPoint(L2, L3);
	PointCloud.push_back(C23);

	//L2��L4�Ľ���
	Point C24 = getCrossPoint(L2, L4);
	PointCloud.push_back(C24);

	//������Բ���������Ĳ���
	Point center = MeanShift(PointCloud);
	double s1 = (p12.y - p1.y) / ((double)p12.x - p1.x), s2 = (center.y - (p1.y+p12.y)/2) / ((double)center.x - (p1.x+p12.x)/2);
	double s3 = (p2.y - p12.y) / ((double)p2.x - p12.x), s4 = (center.y - (p2.y + p12.y) / 2) / ((double)center.x - (p2.x + p12.x) / 2);

	cout << "s1: " << s1 << "  s2:" << s2 << endl;
	double alpha = s1 * s2 - s3 * s4;
	double beta = s2 * s4 * (s3 - s1) + s1 * s3 * (s4 - s2) + (s1 + s2 - s3 - s4);
	cout << "alpha:  " << alpha << endl;
	cout << "beta:   " << beta << endl;
	cout << beta / alpha << endl;
	double K_p = sqrt(1 - beta / alpha), K_n = -K_p;		//K = tan�ѣ�������Բ����ʱ�뷽�����ƫת��

	rho = atan(K_p);
	cout << "�����ǻ������p1��" << p1.x << "," << p1.y << "��" << endl;
	cout << "�����ǻ����е�p12��" << p12.x << "," << p12.y << "��" << endl;
	cout << "�����ǻ����յ�p2��" << p2.x << "," << p2.y << "��" << endl;
	
	cout <<"K:" << K_p << "   s1:" << s1 << "    s2:" << s2 << endl;
	N = sqrt(-(s1 - K_p) * (s2 - K_p) / ((1 + s1 * K_p) * (1 + s2 * K_p)));   //N = b/a, �̰���ͳ�����֮��
	K = K_p;

	return PointCloud;

}


/////////


//Xt���ӵĵ��������Բ����
void drawEllipseFromOneArc(VVP contours, Rectangle rect, int quadrant, VP& result) {
	ushort sz = ushort(contours.size());
	
	float *scores = new float[sz];
	int rect_up = rect.y;
	int rect_down = rect.y + rect.height;
	int rect_left = rect.x;
	int rect_right = rect.x + rect.width;

	cout << "��ı߽���Ϣ" << endl;
	cout << "up = " << rect_up << endl;
	cout << "down = " << rect_down << endl;
	cout << "left = " << rect_left << endl;
	cout << "right = " << rect_right << endl;
	cout << "�ڵ�" << quadrant << "���޹���" << sz << "����,�÷����£�" << endl;
	for (int i = 0; i < sz; i++)
		for (int j = 0; j < contours[i].size(); j++)
		{
			contours[i][j].x += rect.x;
			contours[i][j].y += rect.y;
		}
	//
	switch (quadrant)
	{
	case 1:
		for (int i = 0; i < sz; i++)
		{
			VP& edge_i = contours[i];
			ushort sz_ei = ushort(edge_i.size());		//�ó��������������ٸ���
			Point& pif = edge_i[0];				//�������
			Point& pim = edge_i[sz_ei / 2];		//�м��
			Point& pil = edge_i[sz_ei - 1];		//�յ�
			cout << "����Ϊ" << sz_ei << endl;
			if(sz_ei >= 16 &&
			   pif.x > rect_left && pif.x < rect_right && pif.y > rect_up && pif.y < rect_down &&
			   pim.x > rect_left && pim.x < rect_right && pim.y > rect_up && pim.y < rect_down &&
			   pil.x > rect_left && pil.x < rect_right && pil.y > rect_up && pil.y < rect_down)
			{
				float ave_distance = (min(rect_right - pif.x, pif.y - rect_up) + min(rect_right - pim.x, pim.y - rect_up) + min(rect_right - pil.x, pil.y - rect_up)) / 3.0;
				scores[i] = ave_distance * 0.7 + 10/sz_ei * 0.3;
			}
			else scores[i] = 1000;
			cout << scores[i] << endl;
		}
		break;
	case 2:
		for (int i = 0; i < sz; i++)
		{
			VP& edge_i = contours[i];
			ushort sz_ei = ushort(edge_i.size());
			Point& pif = edge_i[0];				
			Point& pim = edge_i[sz_ei / 2];		
			Point& pil = edge_i[sz_ei - 1];		
			float ave_distance = (min(pif.x-rect_left, pif.y-rect_up) + min(pim.x-rect_left, pim.y-rect_up) + min(pil.x - rect_left, pil.y-rect_up)) / 3.0;
			if (sz_ei >= 16 &&
				pif.x > rect_left && pif.x < rect_right && pif.y > rect_up && pif.y < rect_down &&
				pim.x > rect_left && pim.x < rect_right && pim.y > rect_up && pim.y < rect_down &&
				pil.x > rect_left && pil.x < rect_right && pil.y > rect_up && pil.y < rect_down)
			{
				scores[i] = ave_distance * 0.7 + 10/sz_ei * 0.3;
			}
			else scores[i] = 1000;
			cout << scores[i] << endl;
		}
		break;
	case 3:
		for (int i = 0; i < sz; i++)
		{
			VP& edge_i = contours[i];
			ushort sz_ei = ushort(edge_i.size());
			Point& pif = edge_i[0];				
			Point& pim = edge_i[sz_ei / 2];		
			Point& pil = edge_i[sz_ei - 1];		
			float ave_distance = (min(pif.x - rect_left, rect_down-pif.y) + min(pim.x - rect_left, rect_down-pim.y) + min(pil.x - rect_left, rect_down-pil.y)) / 3.0;
			if (sz_ei >= 16 &&
				pif.x > rect_left && pif.x < rect_right && pif.y > rect_up && pif.y < rect_down &&
				pim.x > rect_left && pim.x < rect_right && pim.y > rect_up && pim.y < rect_down &&
				pil.x > rect_left && pil.x < rect_right && pil.y > rect_up && pil.y < rect_down)
			{
				scores[i] = ave_distance * 0.7 + 10/sz_ei * 0.3;
			}
			else scores[i] = 1000;
			cout << scores[i] << endl;
		}
		break;
	case 4:
		for (int i = 0; i < sz; i++)
		{
			VP& edge_i = contours[i];
			ushort sz_ei = ushort(edge_i.size());	
			Point& pif = edge_i[0];				//�������
			Point& pim = edge_i[sz_ei / 2];		//�м��
			Point& pil = edge_i[sz_ei - 1];		//�յ�
			float ave_distance = (min(rect_right - pif.x, rect_down - pif.y) + min(rect_right - pim.x, rect_down - pim.y) + min(rect_right - pil.x, rect_down - pil.y)) / 3.0;
			if (sz_ei >= 16 &&
				pif.x > rect_left && pif.x < rect_right && pif.y > rect_up && pif.y < rect_down &&
				pim.x > rect_left && pim.x < rect_right && pim.y > rect_up && pim.y < rect_down &&
				pil.x > rect_left && pil.x < rect_right && pil.y > rect_up && pil.y < rect_down)
			{
				scores[i] = ave_distance * 0.7 + 10/sz_ei * 0.3;
			}
			else scores[i] = 1000;
			cout << scores[i] << endl;
		}
		break;
	default:
		break;
	}
	int flag = -1;
	float temp = 100;
	for (int i = 0; i < sz; i++)
	{
		if (scores[i] < temp)
		{
			flag = i;
			temp = scores[i];
		}
	}
	if (flag >= 0)
	{
		cout << "�ҵ������Ż���Ϊ" << contours[flag].size() << ",λ�ڵ�" << quadrant << "����,�ǵ�" << flag << "�����÷�Ϊ��" << scores[flag] << endl;
		//cout << "�������ϵĵ��������£�" << endl;
		//for (int i = 0; i < contours[flag].size(); i++)
		//	cout << contours[flag][i].x << "	" << contours[flag][i].y << endl;
		result = contours[flag];
	}
	else cout << "��������û�з��������Ļ��Σ�" << endl;
}



////////////////
//4.17xt��������ȷ��6�㲢������ϳ�����
void Find_Para_by2(VP arc1, VP arc2, VP arc3, int rand1, int rand2, int rand3, float* EllPara) {
	ushort sz_1 = ushort(arc1.size());
	ushort sz_2 = ushort(arc2.size());
	ushort sz_3 = ushort(arc3.size());
	//��֤arc1
	Point& p11 = arc1[rand1];				    
	Point& p12 = arc1[rand1 + sz_1 / 2];		

	Point& p21 = arc2[rand2];
	Point& p22 = arc2[rand2 + sz_2 / 2];
	
	Point& p31 = arc3[rand3];
	Point& p32 = arc3[rand3 + sz_3 / 2];

	VP test1;
	test1.push_back(p11);
	test1.push_back(p12);
	test1.push_back(p21);
	test1.push_back(p22);
	test1.push_back(p31);
	test1.push_back(p32);

	RotatedRect box = fitEllipse(test1);
	EllipsePara EP_t;
	getEllipsePara(box, EP_t);
	EllPara[0] = EP_t.c.x;
	EllPara[1] = EP_t.c.y;
	EllPara[2] = EP_t.a;
	EllPara[3] = EP_t.b;
	EllPara[4] = EP_t.theta;

}

void choose_6points(VP arc1, VP arc2, VP arc3) {
	ushort sz_1 = ushort(arc1.size());
	ushort sz_2 = ushort(arc2.size());
	ushort sz_3 = ushort(arc3.size());
	
	int rand1 = rand() % (sz_1 / 2);
	int rand2 = rand() % (sz_2 / 2);
	int rand3 = rand() % (sz_3 / 2);

	float testR1[5];
	float testR2[5];
	Find_Para_by2(arc1, arc2, arc3, rand1, rand2, rand3, testR1);
	bool need_another_test = true;
	int loop_times = 0;

	while (need_another_test) {
		int rand11 = rand() % (sz_1 / 2);
		if (rand11 == rand1)
			continue;
		Find_Para_by2(arc1, arc2, arc3, rand11, rand2, rand3, testR2);
		for (int i = 0; i < 5; i++)
		{
			if ((testR1[i] - testR2[i]) / testR1[i] < 0.1 || (testR1[i] - testR2[i]) / testR2[i] < 0.1)
				need_another_test = false;
			else {
				need_another_test = true;
				break;
			}
		}
		if(need_another_test)
			loop_times++;
		if (true) {

		}
	}
}


/////////////////


////////////////
//5.1Xt����Ransac
//����[0,1]֮����Ͼ��ȷֲ�����
double uniformRandom(void)
{
	return (double)rand() / (double)RAND_MAX;
}

//5.12����  ��д����С���˷��������
/*************************************************************************
��С���˷����ֱ�ߣ�y = a*x + b; n������; r-���ϵ��[-1,1],fabs(r)->1,˵��x,y֮�����Թ�ϵ�ã�fabs(r)->0��x,y֮�������Թ�ϵ�����������
a = (n*C - B*D) / (n*A - B*B)
b = (A*D - B*C) / (n*A - B*B)
r = E / F
���У�
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

	for (int i = 0; i<data_n; i++)
	{
		ptsF[i].x;
		ptsF[i].y;
		A += ptsF[i].x * ptsF[i].x;
		B += ptsF[i].x;
		C += ptsF[i].x * ptsF[i].y;
		D += ptsF[i].y;
	}

	// ����б��a�ͽؾ�b
	float a, b, temp = 0;
	if (temp = (data_n*A - B*B))// �жϷ�ĸ��Ϊ0
	{
		a = (data_n*C - B*D) / temp;
		b = (A*D - B*C) / temp;
	}
	else
	{
		a = 1;
		b = 0;
	}

	// �������ϵ��r
	float Xmean, Ymean;
	Xmean = B / data_n;
	Ymean = D / data_n;

	float tempSumXX = 0.0, tempSumYY = 0.0;
	for (int i = 0; i<data_n; i++)
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

//���ݵ㼯���ֱ��ax+by+c=0��resΪ�в�
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

//�õ�ֱ���������������ֱ�߲����㼯�����ѡ2����
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

//ֱ���������������λ�ò���̫��
bool verifyComposition(const vector<Point2d> pts)
{
	cv::Point2d pt1 = pts[0];
	cv::Point2d pt2 = pts[1];
	if (abs(pt1.x - pt2.x) < 5 && abs(pt1.y - pt2.y) < 5)
		return false;

	return true;
}

//RANSACֱ�����
void fitLineRANSAC(vector<Point2f> ptSet, double &a, double &b, double &c, vector<bool> &inlierFlag)
{
	double residual_error = 2.99; //�ڵ���ֵ

	bool stop_loop = false;
	int maximum = 0;  //����ڵ���

	//�����ڵ��ʶ����в�
	inlierFlag = vector<bool>(ptSet.size(), false);
	vector<double> resids_(ptSet.size(), 3);
	int sample_count = 0;
	int N = 500;

	double res = 0;

	// RANSAC
	srand((unsigned int)time(NULL)); //�������������
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

		// ����ֱ�߷���
		calcLinePara(pt_sam, a, b, c, res);
		//�ڵ����
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
		// �ҵ�������ֱ��
		if (inlier_count >= maximum)
		{
			maximum = inlier_count;
			resids_ = residualstemp;
			inlierFlag = inlierstemp;
		}
		// ����RANSAC�����������Լ��ڵ����
		if (inlier_count == 0)
		{
			N = 500;
		}
		else
		{
			double epsilon = 1.0 - double(inlier_count) / (double)ptSet.size(); //Ұֵ�����
			double p = 0.99; //���������д���1���������ĸ���
			double s = 2.0;
			N = int(log(1.0 - p) / log(1.0 - pow((1.0 - epsilon), s)));
		}
		++sample_count;
	}

	//���������ڵ��������ֱ��
	vector<Point2d> pset;
	for (unsigned int i = 0; i < ptSet.size(); i++)
	{
		if (inlierFlag[i])
			pset.push_back(ptSet[i]);
	}

	calcLinePara(pset, a, b, c, res);
}

////////////////


bool myselect1 = true;
bool myselect2 = false;
bool myselect3 = false;
float tCNC = 0.3f;
float tTCNl = 0.005f;//����ֱ̫������

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
	_fMinScore = 0.7f;
	_fMinReliability = 0.5;
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

//����ֵ�� GetFastCenter(
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

		unsigned minPoints = (_uNs < hsize_2) ? _uNs : hsize_2;//ƽ���Ҹ���

		vector<uint> indexes(minPoints);
		if (_uNs < hsize_2)
		{//��������������趨������ƽ������
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
			iota(indexes.begin(), indexes.end(), hsize_2);//ת��unsigned
		}
		for (uint ii = 0; ii<minPoints; ++ii)
		{//��2��ÿ�����ҵ���Ӧƽ���ҵ�
			uint i = indexes[ii];

			float x1 = float(e2[i].x);
			float y1 = float(e2[i].y);

			uint begin = 0;
			uint end = size_1 - 1;
			//һ�����뻡1�ϵĿ�ʼ��������ƽ��������
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
			//�����һ��
			if ((sign_begin + sign_end) != 0)
			{//������û��ƽ�е���
				continue;
			}

			//���ַ�����ƽ����
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
			//���ַ�����ƽ���� end �д�
			med.push_back(Point2f((e1[j].x + x1)* 0.5f, (e1[j].y + y1)* 0.5f));
		}

		if (med.size() < 2)
		{
			data.isValid = false;
			return;
		}

		q2 = GetMedianSlope(med, M12, data.Sa);//�õ�Sa ta=q2 Ma
		cout << "************************************************" << endl;
		cout << "1-2ԭ�����õ����е�����б��Ϊ" << q2 << endl;
		double A, B, C;
		vector<bool> inliers;
		fitLineRANSAC(med, A, B, C, inliers);
		cout << "RANSAC�õ���A��BΪ��" << A << "	" << B << endl;
		cout << "б����" << -A / B << endl;
		cout << "************************************************" << endl;

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
		cout << "************************************************" << endl;
		cout << "2-1ԭ�����õ����е�����б��Ϊ" << q4 << endl;
		double A, B, C;
		vector<bool> inliers;
		fitLineRANSAC(med, A, B, C, inliers);
		cout << "RANSAC�õ���A��BΪ��" << A << "	" << B << endl;
		cout << "б����" << -A / B << endl;
		cout << "************************************************" << endl;
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

//��13���޵Ļ���  ***DetectAfterPreProcessing Detect showEdgeInPic  showEdgeInPic
#define	DISCARD_TCN1
void CNEllipseDetector::DetectEdges13(Mat1b& DP, VVP& points_1, VVP& points_3)
{
	// Vector of connected edge points
	VVP contours;
	int countedges = 0;
	// Labeling 8-connected edge points, discarding edge too small
	Labeling(DP, contours, _iMinEdgeLength);//����һ�����ĵ�ŵ�һ�����򼯺���
	int iContoursSize = int(contours.size());

	// For each edge
	for (int i = 0; i < iContoursSize; ++i)
	{
		VP& edgeSegment = contours[i];
#ifndef DISCARD_CONSTRAINT_OBOX //����Լ��OBOX

		// Selection strategy - Step 1 - See Sect [3.1.2] of the paper
		// Constraint on axes aspect ratio ��չ��Լ��
		RotatedRect oriented = minAreaRect(edgeSegment);//������
		float o_min = min(oriented.size.width, oriented.size.height);

		if (o_min < _fMinOrientedRectSide)
		{
			countedges++;
			continue;
		}
#endif
		// Order edge points of the same arc ��ͬһ���ı�Ե���������
		//�������Ͻ�-���½�
		sort(edgeSegment.begin(), edgeSegment.end(), SortTopLeft2BottomRight);//�������ӵ�
		int iEdgeSegmentSize = unsigned(edgeSegment.size());

		// Get extrema of the arc �õ����ļ�ֵ
		Point& left = edgeSegment[0];
		Point& right = edgeSegment[iEdgeSegmentSize - 1];
#ifndef DISCARD_TCN
#ifndef DISCARD_TCN2
		int flag = 0;
		for (int j = 0; j<iEdgeSegmentSize; j++) {
			Point& mid = edgeSegment[j];
			float data[] = { left.x, left.y, 1,mid.x, mid.y, 1,right.x, right.y, 1 };
			Mat threePoints(Size(3, 3), CV_32FC1, data);
			double ans = determinant(threePoints);

			float dx = 1.0f*(left.x - right.x);
			float dy = 1.0f*(left.y - right.y);
			float edgelength2 = dx*dx + dy*dy;
			//double TCNl=ans/edgelength2;
			double TCNl = ans / (2 * sqrt(edgelength2));
			if (abs(TCNl)>tTCNl) {
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
		float edgelength2 = dx*dx + dy*dy;
		double TCNl = ans / edgelength2;
		//double TCNl=ans/(2*sqrt(edgelength2));
		if (abs(TCNl)<tTCNl) {
			countedges++;
			continue;
		}
#endif
#endif
		// Find convexity - See Sect [3.1.3] of the paper ��͹
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
		//namedWindow("detectEdges", CV_WINDOW_AUTOSIZE);
		//imshow("detectEdges", DP);
		//namedWindow("detectEdges3", CV_WINDOW_AUTOSIZE);
		//imshow("detectEdges3", points_3);
	}
};

//Ƕ����13���޻��� ***Detect showEdgeInPic_zrk showEdgeInPic_zrk
//***Detect OnImage_zrk main_Rect main_allRects case11 case12
#define	DISCARD_TCN1
void CNEllipseDetector::DetectEdges13_zrk(Mat1b& DP, VVP& points_1, VVP& points_3)
{
	// Vector of connected edge points
	VVP contours;
	int countedges = 0;
	// Labeling 8-connected edge points, discarding edge too small  ���ڿ�
	Labeling_zrk(DP, contours, _iMinEdgeLength);//����һ�����ĵ�ŵ�һ�����򼯺���
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
		sort(edgeSegment.begin(), edgeSegment.end(), SortTopLeft2BottomRight);//�������ӵ�
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
		// Find convexity - See Sect [3.1.3] of the paper �ҵ�͹��
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

//��24���޵Ļ���  ***DetectAfterPreProcessing  Detect2563 showEdgeInPic showEdgeInPic
void CNEllipseDetector::DetectEdges24(Mat1b& DN, VVP& points_2, VVP& points_4)
{
	// Vector of connected edge points
	VVP contours;
	int countedges = 0;
	///Labeling 8-connected edge points, discarding edge too small
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
		int flag = 0;
		for (int j = 0; j<iEdgeSegmentSize; j++) {
			Point& mid = edgeSegment[j];
			float data[] = { left.x, left.y, 1,mid.x, mid.y, 1,right.x, right.y, 1 };
			Mat threePoints(Size(3, 3), CV_32FC1, data);
			double ans = determinant(threePoints);

			float dx = 1.0f*(left.x - right.x);
			float dy = 1.0f*(left.y - right.y);
			float edgelength2 = dx*dx + dy*dy;
			//double TCNl=ans/edgelength2;
			double TCNl = ans / (2 * sqrt(edgelength2));
			if (abs(TCNl)>tTCNl) {
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
		float edgelength2 = dx*dx + dy*dy;
		double TCNl = ans / edgelength2;
		//double TCNl=ans/(2*sqrt(edgelength2));
		if (abs(TCNl)<tTCNl) {
			countedges++;
			continue;
		}
#endif
#endif
		// Find convexity - See Sect [3.1.3] of the paper
		//�жϰ�͹��
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

//Ƕ����24���޻���  ***Detect3456 showEdgeInPic_zrk showEdgeInPic_zrk
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
		//�жϰ�͹��
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

//xt
//void fromVectorToMat(vector<Ellipse>& v, Mat &pts)
//{
//	int nbrOfPoints = (int)v.size();
//
//	if (pts.empty())
//		pts.create(2, nbrOfPoints, CV_32F);
//
//	for (int i = 0; i < nbrOfPoints; ++i)
//	{
//		pts.at<float>(0, i) = v[i]._xc;
//		pts.at<float>(1, i) = v[i]._yc;
//	}
//}

//zh
/****************** vectorתMat *********************/
//template<typename Ellipse>
//cv::Mat convertVector2Mat(vector<Ellipse> v, int channels, int rows)
//{
//	cv::Mat mat = cv::Mat(v);//��vector��ɵ��е�mat
//	cv::Mat dest = mat.reshape(channels, rows).clone();//PS������clone()һ�ݣ����򷵻س���
//	return dest;
//}

//zh�������Ҷ���Բ
Mat1b showT(Mat1b image, CNEllipseDetector cned, vector<Ellipse> ellsCned) {
	//Mat image(ori_image, Rect(104, 125, 42, 49));
	Mat1b resultImage = image.clone();
	//Mat3b ori_image = imread(filename);
	//Mat image(ori_image, Rect(104, 125, 42, 49)); //10.19

	// Draw GT ellipses
	/*for (unsigned i = 0; i < ellsCned.size(); ++i)
	{
		Ellipse& e = ellsCned[i];
		Scalar color(0, 0, 1);
		ellipse(resultImage, Point(cvRound(e._xc), cvRound(e._yc)), Size(cvRound(e._a), cvRound(e._b)), e._rad*180.0 / CV_PI, 0.0, 360.0, color, 3);
	}*/

	cned.DrawDetectedEllipses(resultImage, ellsCned);
	return resultImage;
}

// Most important function for detecting ellipses. See Sect[3.2.3] of the paper
//�������λ�������Բ������������õ�����Բ������׼ȷ����
//***Triplets124 Detect  OnImage_zrk main_Rect main_allRects case11 case12
void CNEllipseDetector::FindEllipses(Point2f& center,
	VP& edge_i, VP& edge_j, VP& edge_k,
	EllipseData& data_ij, EllipseData& data_ik,
	vector<Ellipse>& ellipses,
	Rectangle rect
)
{//�������λ�������Բ������������õ�����Բ������׼ȷ����
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
			float q2 = data_ij.Sa[ij1];//�б�Ҫ������

			float q1xq2 = q1*q2;
			// ij1 and ik1
			for (int ik1 = 0; ik1 < sz_ik1; ++ik1)
			{
				float q4 = data_ik.Sa[ik1];//�б�Ҫ������

				float q3xq4 = q3*q4;

				// See Eq. [13-18] in the paper

				float a = (q1xq2 - q3xq4);//gama
				float b = (q3xq4 + 1)*(q1 + q2) - (q1xq2 + 1)*(q3 + q4);//beta
				float Kp = (-b + sqrt(b*b + 4 * a*a)) / (2 * a);//K+
				float zplus = ((q1 - Kp)*(q2 - Kp)) / ((1 + q1*Kp)*(1 + q2*Kp));
				// �ж���Ч��  zplus ��K�������Թ�ϵ
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
					rhoDeg = cvRound((rho * 180 / CV_PI) + 90) % 180; // [0,180)//rho �Ƕȱ�ʾ����һ��
				}

				int iNp = cvRound(Np * 100); // [0, 100]ȡ��

				if (0 <= iNp	&& iNp < ACC_N_SIZE &&
					0 <= rhoDeg	&& rhoDeg < ACC_R_SIZE
					)
				{//Ϊʲô��������ֻ����һ��������  ��Ϊ zplus ��K�������Թ�ϵ
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
	// ���λ��ϵĵ㶼��������A
	//�����Ż�   
	for (ushort l = 0; l < sz_ei; ++l)
	{
		Point& pp = edge_i[l];
		float sk = 1.f / sqrt(Kp*Kp + 1.f);//cos rho
		float x0 = ((pp.x - a0) * sk) + (((pp.y - b0)*Kp) * sk);//���Ż�
		float y0 = -(((pp.x - a0) * Kp) * sk) + ((pp.y - b0) * sk);//���Ż�
		float Ax = sqrt((x0*x0*Np*Np + y0*y0) / ((Np*Np)*(1.f + Kp*Kp)));
		int A = cvRound(abs(Ax / cos(rho)));//�����Ż�
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

	if (rel < _fMinReliability)
	{
		Toc(4); //validation
		return;
	}

	// Assign the new score!
	ell._score = (score + rel) * 0.5f;//need to change

									  // The tentative detection has been confirmed. Save it!
									  //zrk add rect into the ell
	ell._rect = rect;
	ellipses.push_back(ell);

	Toc(4); // Validation
};



//�������λ�������Բ������������õ�����Բ������׼ȷ����
//***Triplets124 Detect  OnImage_zrk main_Rect main_allRects case1 
void CNEllipseDetector::FindEllipses(Point2f& center,
	VP& edge_i, VP& edge_j, VP& edge_k,
	EllipseData& data_ij, EllipseData& data_ik,
	vector<Ellipse>& ellipses
)
{//�������λ�������Բ������������õ�����Բ������׼ȷ����
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
			float q2 = data_ij.Sa[ij1];//�б�Ҫ������

			float q1xq2 = q1 * q2;
			// ij1 and ik1
			for (int ik1 = 0; ik1 < sz_ik1; ++ik1)
			{
				float q4 = data_ik.Sa[ik1];//�б�Ҫ������

				float q3xq4 = q3 * q4;

				// See Eq. [13-18] in the paper

				float a = (q1xq2 - q3xq4);//gama
				float b = (q3xq4 + 1)*(q1 + q2) - (q1xq2 + 1)*(q3 + q4);//beta
				float Kp = (-b + sqrt(b*b + 4 * a*a)) / (2 * a);//K+
				float zplus = ((q1 - Kp)*(q2 - Kp)) / ((1 + q1 * Kp)*(1 + q2 * Kp));
				// �ж���Ч��  zplus ��K�������Թ�ϵ
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
					rhoDeg = cvRound((rho * 180 / CV_PI) + 90) % 180; // [0,180)//rho �Ƕȱ�ʾ����һ��
				}

				int iNp = cvRound(Np * 100); // [0, 100]ȡ��

				if (0 <= iNp && iNp < ACC_N_SIZE &&
					0 <= rhoDeg && rhoDeg < ACC_R_SIZE
					)
				{//Ϊʲô��������ֻ����һ��������  ��Ϊ zplus ��K�������Թ�ϵ
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

	// Find peak in N and K accumulator �ҳ�N��K�ۼ����еķ�ֵ
	int iN = FindMaxN(accN);
	int iK = FindMaxK(accR);

	// Recover real values
	float fK = float(iK);
	float Np = float(iN) * 0.01f;
	float rho = fK * float(CV_PI) / 180.f;	//deg 2 rad
	float Kp = tan(rho);

	// Estimate A. See Eq. [19 - 22] in Sect [3.2.3] of the paper  
	// ���λ��ϵĵ㶼��������A
	//�����Ż�   
	for (ushort l = 0; l < sz_ei; ++l)
	{
		Point& pp = edge_i[l];
		float sk = 1.f / sqrt(Kp*Kp + 1.f);//cos rho
		float x0 = ((pp.x - a0) * sk) + (((pp.y - b0)*Kp) * sk);//���Ż�
		float y0 = -(((pp.x - a0) * Kp) * sk) + ((pp.y - b0) * sk);//���Ż�
		float Ax = sqrt((x0*x0*Np*Np + y0 * y0) / ((Np*Np)*(1.f + Kp * Kp)));
		int A = cvRound(abs(Ax / cos(rho)));//�����Ż�
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
		if (abs(h - 1.f) < _fDistanceToEllipseContour)//DistanceToEllipseContour ����Բ�����ľ���
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
		cout << "û��ͨ���ĵ÷��ǣ�" << counter_on_perimeter << endl;
		return;
	}
	cout << "ͨ����validation1.1" << endl;


	// Compute score
	float score = float(counter_on_perimeter) * invNofPoints;
	if (score < _fMinScore)
	{
		Toc(4); //validation
		cout << "û��ͨ���ĵ÷��ǣ�" << score << endl;
		return;
	}
	cout << "ͨ����validation1.2" << endl;

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

	// This allows to get rid of thick edges ȥӵ������
	float rel = min(1.f, ((di + dj + dk) / (3 * (ell._a + ell._b))));

	if (rel < _fMinReliability)
	{
		Toc(4); //validation
		cout << "δͨ����relֵ�ǣ�" << rel << endl;
		return;
	}
	cout << "ͨ����validation��" << endl;

	// Assign the new score!
	ell._score = (score + rel) * 0.5f;//need to change

									  // The tentative detection has been confirmed. Save it!
									  //zrk add rect into the ell
	ellipses.push_back(ell);



	//CNEllipseDetector cned;
	//Mat3b image1;
	////image1 = Mat(ellipses, true);
	//fromVectorToMat(ellipses, image1);

	//cned.DrawDetectedEllipses(image1, ellipses);//////////////////////////////////////////////////////////////////////////////////////////////
	////cout << "image(opencvĬ�Ϸ��) = " << image1<< ";" << endl << endl;
	//if (true) {
	//	imshow("finding", image1);
	//	//cvSaveImage("C:/Users/zhangruike/Desktop/ell/CNC-Ellipse/image/result/result.jpg",&IplImage(resultImage));
	//	//imwrite("C:/Users/zhangruike/Desktop/ell/CNC-Ellipse/image/result/033_0022_normal.jpg", resultImage);
	//}
	//namedWindow("finding", CV_WINDOW_AUTOSIZE);
	//imshow("finding", image1);

	Toc(4); // Validation
};




// Get the coordinates of the center, given the intersection of the estimated lines. See Fig. [8] in Sect [3.2.3] in the paper.
//***Triplets124 Detect  OnImage_zrk main_Rect main_allRects case11 case12
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

#define T124 pjf,pjm,pjl,pif,pim,pil //ԭ����
#define T231 pil,pim,pif,pjf,pjm,pjl
#define T342 pif,pim,pil,pjf,pjm,pjl
#define T413 pif,pim,pil,pjl,pjm,pjf


//��֤���Խ�Բ����͹��:i=1, j=2, k=4
// Verify triplets of arcs with convexity: i=1, j=2, k=4
void CNEllipseDetector::Triplets124(VVP& pi,
	VVP& pj,
	VVP& pk,
	unordered_map<uint, EllipseData>& data,
	vector<Ellipse>& ellipses,
	Rectangle rect
)
{
	// get arcs length  �ó�ÿһ�����ж�������
	ushort sz_i = ushort(pi.size());
	ushort sz_j = ushort(pj.size());
	ushort sz_k = ushort(pk.size());

	// For each edge i  �ڵ�һ�����е�ϵ�л�������
	for (ushort i = 0; i < sz_i; ++i)
	{
		VP& edge_i = pi[i];
		ushort sz_ei = ushort(edge_i.size());		//�ó��������������ٸ���

		Point& pif = edge_i[0];				//�������
		Point& pim = edge_i[sz_ei / 2];		//�м��
		Point& pil = edge_i[sz_ei - 1];		//�յ�

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
			// ����޶�����  2se se1//pil,pim,pif,pjf,pjm,pjl pjf,pjm,pjl,pif,pim,pil
			if (myselect1&&fabs(value4SixPoints(T124) - 1)>tCNC)
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
				// ����޶�����2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1)>tCNC)
				{
					continue;
				}
				// ����޶�����3
				if (myselect3&&fabs(value4SixPoints(pjf, pjm, pjl, pkf, pkm, pkl) - 1)>tCNC)
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
				FindEllipses(center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses, rect);
			}
		}
	}
};


//***Detect  OnImage_zrk main_Rect main_allRects case11 case12
void CNEllipseDetector::Triplets231(VVP& pi,
	VVP& pj,
	VVP& pk,
	unordered_map<uint, EllipseData>& data,
	vector<Ellipse>& ellipses,
	Rectangle rect
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
//����_Լ��
#ifndef DISCARD_CONSTRAINT_CNC
			//���Լ��
			// ����޶����� 2es se3 //pif,pim,pil,pjf,pjm,pjl pil,pim,pif,pjf,pjm,pjl
			if (myselect1&&fabs(value4SixPoints(T231) - 1)>tCNC)
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
				// ����޶�����2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1)>tCNC)
				{
					continue;
				}
				// ����޶�����3
				if (myselect3&&fabs(value4SixPoints(pjf, pjm, pjl, pkf, pkm, pkl) - 1)>tCNC)
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

				FindEllipses(center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses, rect);

			}
		}
	}
};


void CNEllipseDetector::Triplets342(VVP& pi,
	VVP& pj,
	VVP& pk,
	unordered_map<uint, EllipseData>& data,
	vector<Ellipse>& ellipses,
	Rectangle rect
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
			// ����޶����� 3se se4 // pil,pim,pif,pjf,pjm,pjl pif,pim,pil,pjf,pjm,pjl
			if (myselect1&&fabs(value4SixPoints(T342) - 1)>tCNC)
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
				// ����޶�����2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1)>tCNC)
				{
					continue;
				}
				// ����޶�����3
				if (myselect3&&fabs(value4SixPoints(pjf, pjm, pjl, pkf, pkm, pkl) - 1)>tCNC)
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
				FindEllipses(center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses, rect);
			}
		}

	}
};



void CNEllipseDetector::Triplets413(VVP& pi,
	VVP& pj,
	VVP& pk,
	unordered_map<uint, EllipseData>& data,
	vector<Ellipse>& ellipses,
	Rectangle rect
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
			// ����޶����� 4se es1//pif,pim,pil,pjf,pjm,pjl pil,pim,pif,pjl,pjm,pjf pif,pim,pil,pjl,pjm,pjf
			if (myselect1&&fabs(value4SixPoints(T413) - 1)>tCNC)
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
				// ����޶�����2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1)>tCNC)
				{
					continue;
				}
				// ����޶�����3
				if (myselect3&&fabs(value4SixPoints(pjf, pjm, pjl, pkf, pkm, pkl) - 1)>tCNC)
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

				FindEllipses(center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses, rect);

			}
		}
	}
};
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
		Point& pim = edge_i[sz_ei / 2];
		Point& pil = edge_i[sz_ei - 1];

		// 1,2 -> reverse 1, swap ����1������
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
			// ����޶�����  2se se1//pil,pim,pif,pjf,pjm,pjl pjf,pjm,pjl,pif,pim,pil
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


				//CNEllipseDetector cned;
				//Mat3b image1;
				////cned.DrawDetectedEllipses(image1, VP);//////////////////////////////////////////////////////////////////////////////////////////////
				////namedWindow("CONSTRAINTS1", CV_WINDOW_AUTOSIZE);
				////imshow("CONSTRAINTS1", image1);



#ifndef DISCARD_CONSTRAINT_CNC
				// ����޶�����2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}

				//cned.DrawDetectedEllipses(image1, ellipses);//////////////////////////////////////////////////////////////////////////////////////////////
				//namedWindow("CONSTRAINTS2", CV_WINDOW_AUTOSIZE);
				//imshow("CONSTRAINTS2", image1);

				// ����޶�����3
				if (myselect3&&fabs(value4SixPoints(pjf, pjm, pjl, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
#endif
				//cned.DrawDetectedEllipses(image1, ellipses);//////////////////////////////////////////////////////////////////////////////////////////////
				//namedWindow("CONSTRAINTS3", CV_WINDOW_AUTOSIZE);
				//imshow("CONSTRAINTS3", image1);

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

				//cned.DrawDetectedEllipses(image1, ellipses);//////////////////////////////////////////////////////////////////////////////////////////////
				//namedWindow("CONSTRAINTS3", CV_WINDOW_AUTOSIZE);
				//imshow("CONSTRAINTS3", image1);

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
				

				//cv::Mat image1 = convertVector2Mat<uchar>(ellipses, 1, 6);//������תΪ1ͨ����4�е�Mat����
				//CNEllipseDetector cned;
				//Mat3b image1;
				//////image1 = Mat(ellipses, true);
				////fromVectorToMat(ellipses, image1);

				//cned.DrawDetectedEllipses(image1, ellipses);//////////////////////////////////////////////////////////////////////////////////////////////
				//////cout << "image(opencvĬ�Ϸ��) = " << image1<< ";" << endl << endl;
				////if (true) {
				////	imshow("finding", image1);
				////	//cvSaveImage("C:/Users/zhangruike/Desktop/ell/CNC-Ellipse/image/result/result.jpg",&IplImage(resultImage));
				////	//imwrite("C:/Users/zhangruike/Desktop/ell/CNC-Ellipse/image/result/033_0022_normal.jpg", resultImage);
				////}
				//namedWindow("finding", CV_WINDOW_AUTOSIZE);
				//imshow("finding", image1);
				
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
			//���Լ��
			// ����޶����� 2es se3 //pif,pim,pil,pjf,pjm,pjl pil,pim,pif,pjf,pjm,pjl
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
				// ����޶�����2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
				// ����޶�����3
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
			// ����޶����� 3se se4 // pil,pim,pif,pjf,pjm,pjl pif,pim,pil,pjf,pjm,pjl
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
				// ����޶�����2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
				// ����޶�����3
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
			// ����޶����� 4se es1//pif,pim,pil,pjf,pjm,pjl pil,pim,pif,pjl,pjm,pjf pif,pim,pil,pjl,pjm,pjf
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
				// ����޶�����2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
				// ����޶�����3
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


//zhչʾ������
void CNEllipseDetector::Triplets124(Mat1b& I, VVP& pi,
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
		Point& pim = edge_i[sz_ei / 2];
		Point& pil = edge_i[sz_ei - 1];

		// 1,2 -> reverse 1, swap ����1������
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
			CNEllipseDetector cned;
			Mat1b image1 = showT(I, cned, ellipses);
			namedWindow("Triplets124_0", CV_WINDOW_AUTOSIZE);
			imshow("Triplets124_0", image1);

#ifndef DISCARD_CONSTRAINT_CNC
			// ����޶�����  2se se1//pil,pim,pif,pjf,pjm,pjl pjf,pjm,pjl,pif,pim,pil
			if (myselect1&&fabs(value4SixPoints(T124) - 1) > tCNC)
			{
				continue;
			}
#endif

			image1 = showT(I, cned, ellipses);
			namedWindow("Triplets124_1", CV_WINDOW_AUTOSIZE);
			imshow("Triplets124_1", image1);

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


				//CNEllipseDetector cned;
				//Mat3b image1;
				////cned.DrawDetectedEllipses(image1, VP);//////////////////////////////////////////////////////////////////////////////////////////////
				////namedWindow("CONSTRAINTS1", CV_WINDOW_AUTOSIZE);
				////imshow("CONSTRAINTS1", image1);



#ifndef DISCARD_CONSTRAINT_CNC
				// ����޶�����2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}


				// ����޶�����3
				if (myselect3&&fabs(value4SixPoints(pjf, pjm, pjl, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
#endif

				image1 = showT(I, cned, ellipses);
				namedWindow("Triplets124_2", CV_WINDOW_AUTOSIZE);
				imshow("Triplets124_2", image1);

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

				//cned.DrawDetectedEllipses(image1, ellipses);//////////////////////////////////////////////////////////////////////////////////////////////
				//namedWindow("CONSTRAINTS3", CV_WINDOW_AUTOSIZE);
				//imshow("CONSTRAINTS3", image1);

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


				
				image1 = showT(I, cned, ellipses);
				namedWindow("Triplets124_3", CV_WINDOW_AUTOSIZE);
				imshow("Triplets124_3", image1);

			}
		}
	}
};



void CNEllipseDetector::Triplets231(Mat1b& I, VVP& pi,
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
			CNEllipseDetector cned;
			Mat1b image1 = showT(I, cned, ellipses);
			namedWindow("Triplets231_0", CV_WINDOW_AUTOSIZE);
			imshow("Triplets231_0", image1);

#ifndef DISCARD_CONSTRAINT_CNC
			//���Լ��
			// ����޶����� 2es se3 //pif,pim,pil,pjf,pjm,pjl pil,pim,pif,pjf,pjm,pjl
			if (myselect1&&fabs(value4SixPoints(T231) - 1) > tCNC)
			{
				continue;
			}
#endif

			image1 = showT(I, cned, ellipses);
			namedWindow("Triplets231_1", CV_WINDOW_AUTOSIZE);
			imshow("Triplets231_1", image1);

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
				// ����޶�����2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
				// ����޶�����3
				if (myselect3&&fabs(value4SixPoints(pjf, pjm, pjl, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
#endif

				image1 = showT(I, cned, ellipses);
				namedWindow("Triplets231_2", CV_WINDOW_AUTOSIZE);
				imshow("Triplets231_2", image1);

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


				image1 = showT(I, cned, ellipses);
				namedWindow("Triplets231_3", CV_WINDOW_AUTOSIZE);
				imshow("Triplets231_3", image1);
			}
		}
	}
};


void CNEllipseDetector::Triplets342(Mat1b& I, VVP& pi,
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
			// ����޶����� 3se se4 // pil,pim,pif,pjf,pjm,pjl pif,pim,pil,pjf,pjm,pjl
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
				// ����޶�����2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
				// ����޶�����3
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
				FindEllipses(center, edge_i, edge_j, edge_k, data_ij, data_ik, ellipses);
			}
		}

	}
};



void CNEllipseDetector::Triplets413(Mat1b& I, VVP& pi,
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
			// ����޶����� 4se es1//pif,pim,pil,pjf,pjm,pjl pil,pim,pif,pjl,pjm,pjf pif,pim,pil,pjl,pjm,pjf
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
				// ����޶�����2
				if (myselect2&&fabs(value4SixPoints(pif, pim, pil, pkf, pkm, pkl) - 1) > tCNC)
				{
					continue;
				}
				// ����޶�����3
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



//ȥ��С����
//�޵���
void CNEllipseDetector::RemoveShortEdges(Mat1b& edges, Mat1b& clean)
{
	VVP contours;

	// Labeling and contraints on length
	Labeling_zrk(edges, contours, _iMinEdgeLength);

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

////Xt���ӵ��񻯺���
//void sharpen(const cv::Mat &image, cv::Mat &result) {
//	// allocate if necessary
//	result.create(image.size(), image.type());
//	for (int j = 1; j<image.rows - 1; j++) { // for all rows
//											 // (except first and last)
//		const uchar* previous =
//			image.ptr<const uchar>(j - 1); // previous row
//		const uchar* current =
//			image.ptr<const uchar>(j); // current row
//		const uchar* next =
//			image.ptr<const uchar>(j + 1); // next row
//		uchar* output = result.ptr<uchar>(j); // output row
//		for (int i = 1; i<image.cols - 1; i++) {
//			*output++ = cv::saturate_cast<uchar>(
//				5 * current[i] - current[i - 1]
//				- current[i + 1] - previous[i] - next[i]);
//		}
//	}
//}

//Xt���ӵ����ź���
void imageResize(Mat1b& I, int scale)
{
	IplImage* scr = &IplImage(I);
	IplImage* dst = &IplImage(I);
	CvSize dst_cvsize;
	dst_cvsize.width = (int)(scr->width*scale);
	dst_cvsize.height = (int)(scr->height*scale);
	dst = cvCreateImage(dst_cvsize, scr->depth, scr->nChannels);
	cvResize(scr, dst, CV_INTER_NN);
	I = cvarrToMat(dst);
}


////zh���ӵ�˫������ֵ
float BiCubicPoly(float x)
{
	float abs_x = abs(x);
	float a = -0.5;
	if (abs_x <= 1.0)
	{
		return (a + 2)*pow(abs_x, 3) - (a + 3)*pow(abs_x, 2) + 1;
	}
	else if (abs_x < 2.0)
	{
		return a*pow(abs_x, 3) - 5 * a*pow(abs_x, 2) + 8 * a*abs_x - 4 * a;
	}
	else
		return 0.0;
}

//�޵���
void MyScaleBiCubicInter(Mat& src, Mat& dst, float TransMat[3][3])
{
	CV_Assert(src.data);
	CV_Assert(src.depth() != sizeof(uchar));

	// calculate margin point of dst image
	float left = 0;
	float right = 0;
	float top = 0;
	float down = 0;

	float x = src.cols * 1.0f;
	float y = 0.0f;
	float u1 = x * TransMat[0][0] + y * TransMat[0][1];
	float v1 = x * TransMat[1][0] + y * TransMat[1][1];
	x = src.cols * 1.0f;
	y = src.rows * 1.0f;
	float u2 = x * TransMat[0][0] + y * TransMat[0][1];
	float v2 = x * TransMat[1][0] + y * TransMat[1][1];
	x = 0.0f;
	y = src.rows * 1.0f;
	float u3 = x * TransMat[0][0] + y * TransMat[0][1];
	float v3 = x * TransMat[1][0] + y * TransMat[1][1];

	left = min(min(min(0.0f, u1), u2), u3);
	right = max(max(max(0.0f, u1), u2), u3);
	top = min(min(min(0.0f, v1), v2), v3);
	down = max(max(max(0.0f, v1), v2), v3);

	// create dst image
	dst.create(int(abs(right - left)), int(abs(down - top)), src.type());

	CV_Assert(dst.channels() == src.channels());
	int channels = dst.channels();

	int i, j;
	uchar* p;
	uchar* q0;
	uchar* q1;
	uchar* q2;
	uchar* q3;
	for (i = 0; i < dst.rows; ++i)
	{
		p = dst.ptr<uchar>(i);
		for (j = 0; j < dst.cols; ++j)
		{
			// 
			x = (j + left) / TransMat[0][0];
			y = (i + top) / TransMat[1][1];

			int x0 = int(x) - 1;
			int y0 = int(y) - 1;
			int x1 = int(x);
			int y1 = int(y);
			int x2 = int(x) + 1;
			int y2 = int(y) + 1;
			int x3 = int(x) + 2;
			int y3 = int(y) + 2;

			if ((x0 >= 0) && (x3 < src.cols) && (y0 >= 0) && (y3 < src.rows))
			{
				q0 = src.ptr<uchar>(y0);
				q1 = src.ptr<uchar>(y1);
				q2 = src.ptr<uchar>(y2);
				q3 = src.ptr<uchar>(y3);

				float dist_x0 = BiCubicPoly(x - x0);
				float dist_x1 = BiCubicPoly(x - x1);
				float dist_x2 = BiCubicPoly(x - x2);
				float dist_x3 = BiCubicPoly(x - x3);
				float dist_y0 = BiCubicPoly(y - y0);
				float dist_y1 = BiCubicPoly(y - y1);
				float dist_y2 = BiCubicPoly(y - y2);
				float dist_y3 = BiCubicPoly(y - y3);

				float dist_x0y0 = dist_x0 * dist_y0;
				float dist_x0y1 = dist_x0 * dist_y1;
				float dist_x0y2 = dist_x0 * dist_y2;
				float dist_x0y3 = dist_x0 * dist_y3;
				float dist_x1y0 = dist_x1 * dist_y0;
				float dist_x1y1 = dist_x1 * dist_y1;
				float dist_x1y2 = dist_x1 * dist_y2;
				float dist_x1y3 = dist_x1 * dist_y3;
				float dist_x2y0 = dist_x2 * dist_y0;
				float dist_x2y1 = dist_x2 * dist_y1;
				float dist_x2y2 = dist_x2 * dist_y2;
				float dist_x2y3 = dist_x2 * dist_y3;
				float dist_x3y0 = dist_x3 * dist_y0;
				float dist_x3y1 = dist_x3 * dist_y1;
				float dist_x3y2 = dist_x3 * dist_y2;
				float dist_x3y3 = dist_x3 * dist_y3;

				switch (channels)
				{
				case 1:
				{
					break;
				}
				case 3:
				{
					p[3 * j] = (uchar)(q0[3 * x0] * dist_x0y0 +
						q1[3 * x0] * dist_x0y1 +
						q2[3 * x0] * dist_x0y2 +
						q3[3 * x0] * dist_x0y3 +
						q0[3 * x1] * dist_x1y0 +
						q1[3 * x1] * dist_x1y1 +
						q2[3 * x1] * dist_x1y2 +
						q3[3 * x1] * dist_x1y3 +
						q0[3 * x2] * dist_x2y0 +
						q1[3 * x2] * dist_x2y1 +
						q2[3 * x2] * dist_x2y2 +
						q3[3 * x2] * dist_x2y3 +
						q0[3 * x3] * dist_x3y0 +
						q1[3 * x3] * dist_x3y1 +
						q2[3 * x3] * dist_x3y2 +
						q3[3 * x3] * dist_x3y3);

					p[3 * j + 1] = (uchar)(q0[3 * x0 + 1] * dist_x0y0 +
						q1[3 * x0 + 1] * dist_x0y1 +
						q2[3 * x0 + 1] * dist_x0y2 +
						q3[3 * x0 + 1] * dist_x0y3 +
						q0[3 * x1 + 1] * dist_x1y0 +
						q1[3 * x1 + 1] * dist_x1y1 +
						q2[3 * x1 + 1] * dist_x1y2 +
						q3[3 * x1 + 1] * dist_x1y3 +
						q0[3 * x2 + 1] * dist_x2y0 +
						q1[3 * x2 + 1] * dist_x2y1 +
						q2[3 * x2 + 1] * dist_x2y2 +
						q3[3 * x2 + 1] * dist_x2y3 +
						q0[3 * x3 + 1] * dist_x3y0 +
						q1[3 * x3 + 1] * dist_x3y1 +
						q2[3 * x3 + 1] * dist_x3y2 +
						q3[3 * x3 + 1] * dist_x3y3);

					p[3 * j + 2] = (uchar)(q0[3 * x0 + 2] * dist_x0y0 +
						q1[3 * x0 + 2] * dist_x0y1 +
						q2[3 * x0 + 2] * dist_x0y2 +
						q3[3 * x0 + 2] * dist_x0y3 +
						q0[3 * x1 + 2] * dist_x1y0 +
						q1[3 * x1 + 2] * dist_x1y1 +
						q2[3 * x1 + 2] * dist_x1y2 +
						q3[3 * x1 + 2] * dist_x1y3 +
						q0[3 * x2 + 2] * dist_x2y0 +
						q1[3 * x2 + 2] * dist_x2y1 +
						q2[3 * x2 + 2] * dist_x2y2 +
						q3[3 * x2 + 2] * dist_x2y3 +
						q0[3 * x3 + 2] * dist_x3y0 +
						q1[3 * x3 + 2] * dist_x3y1 +
						q2[3 * x3 + 2] * dist_x3y2 +
						q3[3 * x3 + 2] * dist_x3y3);

					float thre = 198.0f;
					if ((abs(p[3 * j] - q1[3 * x1]) > thre) || (abs(p[3 * j + 1] - q1[3 * x1 + 1]) > thre) ||
						(abs(p[3 * j + 2] - q1[3 * x1 + 2]) > thre))
					{
						p[3 * j] = q1[3 * x1];
						p[3 * j + 1] = q1[3 * x1 + 1];
						p[3 * j + 2] = q1[3 * x1 + 2];
					}
					break;
				}
				}
			}
		}
	}
}


//zh���ӵ�ֱ��ͼ��ǿ
//int gray[256] = { 0 };  //��¼ÿ���Ҷȼ����µ����ظ���
//double gray_prob[256] = { 0 };  //��¼�Ҷȷֲ��ܶ�
//double gray_distribution[256] = { 0 };  //��¼�ۼ��ܶ�
//int gray_equal[256] = { 0 };  //���⻯��ĻҶ�ֵ
//int gray_sum = 0;  //��������
//Mat equalize_hist(Mat& input)
//{
//	Mat output = input.clone();
//	gray_sum = input.cols * input.rows;
//
//	//ͳ��ÿ���Ҷ��µ����ظ���
//	for (int i = 0; i < input.rows; i++)
//	{
//		uchar* p = input.ptr<uchar>(i);
//		for (int j = 0; j < input.cols; j++)
//		{
//			int vaule = p[j];
//			gray[vaule]++;
//		}
//	}
//
//	//ͳ�ƻҶ�Ƶ��
//	for (int i = 0; i < 256; i++)
//	{
//		gray_prob[i] = ((double)gray[i] / gray_sum);
//	}
//
//	//�����ۼ��ܶ�
//	gray_distribution[0] = gray_prob[0];
//	for (int i = 1; i < 256; i++)
//	{
//		gray_distribution[i] = gray_distribution[i - 1] + gray_prob[i];
//	}
//
//	//���¼�����⻯��ĻҶ�ֵ���������롣�ο���ʽ��(N-1)*T+0.5
//	for (int i = 0; i < 256; i++)
//	{
//		gray_equal[i] = (uchar)(255 * gray_distribution[i] + 0.5);
//	}
//
//	//ֱ��ͼ���⻯,����ԭͼÿ���������ֵ
//	for (int i = 0; i < output.rows; i++)
//	{
//		uchar* p = output.ptr<uchar>(i);
//		for (int j = 0; j < output.cols; j++)
//		{
//			p[j] = gray_equal[p[j]];
//		}
//	}
//
//	return output;
//}




//***Detect��OnImage_zrk main_Rect main_allRects case11 case12�� showϵ�У����ã�
void CNEllipseDetector::PrePeocessing(Mat1b& I,
	Mat1b& DP,
	Mat1b& DN
)
{

	Tic(0); //edge detection
	// Mid smooth
	//medianBlur(I,I,3);   //��ֵ�˲�
	//namedWindow("medianBlur", CV_WINDOW_AUTOSIZE);
	//imshow("medianBlur", I);
	// Smooth image
	//�δ������Ǹ�˹�˲������ﻻΪ����ֵ�˲�
	GaussianBlur(I, I, _szPreProcessingGaussKernelSize, _dPreProcessingGaussSigma);    //��˹�˲�
	namedWindow("GaussianBlur", CV_WINDOW_AUTOSIZE);
	imshow("GaussianBlur", I);
	// Temp variables
	Mat1b E;				//edge mask
	Mat1s DX, DY;			//sobel derivatives


	//�Ŵ�
	//imageResize(I, 2);


	//˫������ֵ���ź���
	Mat dst;
	float transMat[3][3] = { { 2.0, 0, 0 },{ 0, 2.0, 0 },{ 0, 0, 1 } };

	MyScaleBiCubicInter(I, dst, transMat);
	namedWindow("out_image", CV_WINDOW_AUTOSIZE);
	imshow("out_image", I);

	//��
	//sharpen(I, I);

	//ֱ��ͼ��ǿ
	//I = equalize_hist(I);


	// Detect edges
	Canny2(I, E, DX, DY, 0, 10, 3, false);
	namedWindow("Canny2", CV_WINDOW_AUTOSIZE);
	imshow("Canny2", I);

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
	
	//namedWindow("preprocessing", CV_WINDOW_AUTOSIZE);
	//imshow("preprocessing", I);
};

//�޵���
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

//zh
void CNEllipseDetector::DetectAfterPreProcessing_zh(vector<Ellipse>& ellipses, Mat1b& E, Mat1f& PHI, Rectangle rect)
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
	DetectEdges13_zrk(DP, points_1, points_3);
	DetectEdges24_zrk(DN, points_2, points_4);


	// Find triplets
	Triplets124(points_1, points_2, points_4, centers, ellipses, rect);
	Triplets231(points_2, points_3, points_1, centers, ellipses, rect);
	Triplets342(points_3, points_4, points_2, centers, ellipses, rect);
	Triplets413(points_4, points_1, points_3, centers, ellipses, rect);


	// Sort detected ellipses with respect to score
	sort(ellipses.begin(), ellipses.end());

	//free accumulator memory
	delete[] accN;
	delete[] accR;
	delete[] accA;

	//cluster detections
	//ClusterEllipses(ellipses);
};

//1.18  Xt���� validation
void CNEllipseDetector::ValidationForOneArcEllipse(Ellipse ell, VP& edge, vector<Ellipse>& ellipses)
{
	countsOfFindEllipse++;
	Tac(4); //validation
	size_t sz_e = edge.size();
			// Get the score. See Sect [3.3.1] in the paper

			// Find the number of edge pixel lying on the ellipse
	float _cos = cos(-ell._rad);
	float _sin = sin(-ell._rad);

	float invA2 = 1.f / (ell._a * ell._a);
	float invB2 = 1.f / (ell._b * ell._b);

	float invNofPoints = 1.f / float(sz_e);
	int counter_on_perimeter = 0;

	for (ushort l = 0; l < sz_e; ++l)
	{
		float tx = float(edge[l].x) - ell._xc;
		float ty = float(edge[l].y) - ell._yc;
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
		return;
	}
	cout << "***********ͨ��validation1.1" << endl;

	// Compute score
	float score = float(counter_on_perimeter) * invNofPoints;
	if (score < 0.5f)
	{
		//validation
		Toc(4);
		cout << "û��ͨ���ĵ÷��ǣ�" << score << endl;
		return;
	}
	cout << "***********ͨ��validation1.2" << endl;

	//validation�ĵڶ�����
	float d;
	{
		Point2f p1(float(edge[0].x), float(edge[0].y));
		Point2f p2(float(edge[sz_e - 1].x), float(edge[sz_e - 1].y));
		p1.x -= ell._xc;
		p1.y -= ell._yc;
		p2.x -= ell._xc;
		p2.y -= ell._yc;
		Point2f r1((p1.x*_cos - p1.y*_sin), (p1.x*_sin + p1.y*_cos));
		Point2f r2((p2.x*_cos - p2.y*_sin), (p2.x*_sin + p2.y*_cos));
		d = abs(r2.x - r1.x) + abs(r2.y - r1.y);
	}

	// This allows to get rid of thick edges
	float rel = min(1.f, (d / (ell._a + ell._b)));

	if (rel < _fMinReliability)
	{
		Toc(4); //validation
		cout << "δͨ����relֵΪ��" <<  rel <<endl;
		return;
	}
	cout << "***********ͨ��validation2.1" << endl;

	// Assign the new score!
	ell._score = (score + rel) * 0.5f;//need to change

									  // The tentative detection has been confirmed. Save it!


	ellipses.push_back(ell);
	cout << "***********************�ѳɹ�ͨ��validation���һ����Բ��������****************" << endl;

};


//zrk
//**OnImage_zrk main_Rect main_allRects case11 case12
//Ԥ����+��⻡��+���λ��ж��Ҽ�����Բ����
void CNEllipseDetector::Detect(Mat1b& I, vector<Ellipse>& ellipses, Rectangle rect, double* EllPara)
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

	Tac(1); //preprocessing
			// Detect edges and find convexities
	DetectEdges13_zrk(DP, points_1, points_3);
	DetectEdges24_zrk(DN, points_2, points_4);

	Toc(1); //preprocessing

	////////����Ϊ�����ĵ�/˫������
	//Xt�������������������޵Ļ�ʱ
	
	if (points_1.size() == 0 && (points_2.size() == 0 || points_3.size() == 0 || points_4.size() == 0) ||
		points_2.size() == 0 && (points_3.size() == 0 || points_4.size() == 0) ||
		points_3.size() == 0 && points_4.size() == 0) {
		VP arc, arc1, arc2, arc3, arc4;
		
		int arc_counter = 0;
		if (points_1.size() != 0)
		{
			drawEllipseFromOneArc(points_1, rect, 1, arc1);
			if(arc1.size() != 0)
				arc_counter++;
		}
		if (points_2.size() != 0)
		{
			drawEllipseFromOneArc(points_2, rect, 2, arc2);
			if (arc2.size() != 0)
				arc_counter++;
		}
		if (points_3.size() != 0)
		{
			drawEllipseFromOneArc(points_3, rect, 3, arc3);
			if (arc3.size() != 0)
				arc_counter++;
		}
		if (points_4.size() != 0)
		{
			drawEllipseFromOneArc(points_4, rect, 4, arc4);
			if (arc4.size() != 0)
				arc_counter++;
		}

		if (arc_counter < 2)
		{
			cout << "���󣡲����������mean-shift��ϵ�����" << endl;
			getchar();
		}
			
		else cout << "#########����" << arc_counter << "�������л�" << endl;

		arc.insert(arc.end(), arc1.begin(), arc1.end());
		arc.insert(arc.end(), arc2.begin(), arc2.end());
		arc.insert(arc.end(), arc3.begin(), arc3.end());
		arc.insert(arc.end(), arc4.begin(), arc4.end());
		cout << "���֮�󻡶��Ϲ���" << arc.size() << "����" << endl;

		VP arcA, arcB;
		if (arc1.size() != 0 && arc2.size() != 0)
		{
			arcA = arc1;
			arcB = arc2;
		}
		if (arc1.size() != 0 && arc3.size() != 0)
		{
			arcA = arc1;
			arcB = arc3;
		}
		if (arc1.size() != 0 && arc4.size() != 0)
		{
			arcA = arc1;
			arcB = arc4;
		}
		if (arc2.size() != 0 && arc3.size() != 0)
		{
			arcA = arc2;
			arcB = arc3;
		}
		if (arc2.size() != 0 && arc4.size() != 0)
		{
			arcA = arc2;
			arcB = arc4;
		}
		if (arc3.size() != 0 && arc4.size() != 0)
		{
			arcA = arc3;
			arcB = arc4;
		}
		
		double N, K, rho;
		VP PointCloud = FindPointCloud(arcA, arcB,N, K, rho);
		cout << "���ҵ�7��������!!!!!!" << N << "	" << K << endl;
		Point center;
 		center = MeanShift(PointCloud);
		cout << "ͨ��mean-shift�ҵ������ĵ���(" << center.x << "," << center.y << ")************"  << endl;
		Point p0 = arcA[arcA.size() / 2];
		double x_0 = ((p0.x - center.x) + (p0.y - center.y) * K) / sqrt(K * K + 1),
			   y_0 = ((p0.y - center.y) + (p0.x - center.x) * K) / sqrt(K * K + 1);
		double Ax = sqrt((x_0 * x_0 * N * N + y_0 * y_0) / (N * N * (1 + K * K)));
		cout << Ax << endl;
		double a = Ax * acos(rho), b = a * N;
		cout << "����mean-shift������������в����ֱ���" << a << "	" << b << "	" << rho << endl;

		//��˹�ֽ�
		//LeastSquareFittingEllipse_LU(arc, FittingResult);
		

		//opencv�Դ��ĺ���
		RotatedRect box = fitEllipse(arc);
		EllipsePara EP_t;
		getEllipsePara(box, EP_t);
		EllPara[0] = EP_t.c.x;
		EllPara[1] = EP_t.c.y;
		EllPara[2] = EP_t.a;
		EllPara[3] = EP_t.b;
		EllPara[4] = EP_t.theta;


		//QR�ֽ�
		//Point2f result = LeastSquareFittingEllipse(arc);
		/*
		FittingResult[0] = result.x;
		FittingResult[1] = result.y;
		*/
		////////1.18  Xt   validation
		Ellipse ell(EllPara[0], EllPara[1], EllPara[2], EllPara[3], fmod(EllPara[4] + float(CV_PI)*2.f, float(CV_PI)));
		ValidationForOneArcEllipse(ell, arc, ellipses);
	}


	Tic(2); //grouping
			//find triplets
	Triplets124(points_1, points_2, points_4, centers, ellipses, rect);
	Triplets231(points_2, points_3, points_1, centers, ellipses, rect);
	Triplets342(points_3, points_4, points_2, centers, ellipses, rect);
	Triplets413(points_4, points_1, points_3, centers, ellipses, rect);

	
	//namedWindow("Triplets", CV_WINDOW_AUTOSIZE);
	//imshow("Triplets", I);
	
	Toc(2); //grouping	
			// time estimation, validation inside
	_times[2] -= (_times[3] + _times[4]);//����׼ȷ

	Tac(4); //validation
			// Sort detected ellipses with respect to score
	sort(ellipses.begin(), ellipses.end());//���շ������� ������ͬ Խ�ӽ�ԲԽС
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



void CNEllipseDetector::Detect(Mat1b& I, vector<Ellipse>& ellipses)
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


			// time estimation, validation  inside

	Tic(2); //grouping
			//find triplets
	namedWindow("1", CV_WINDOW_AUTOSIZE);
	//imshow("1", I);

	Triplets124(I,points_1, points_2, points_4, centers, ellipses);
	Triplets231(I,points_2, points_3, points_1, centers, ellipses);
	Triplets342(I,points_3, points_4, points_2, centers, ellipses);
	Triplets413(I,points_4, points_1, points_3, centers, ellipses);


	CNEllipseDetector cned;
	Mat1b image1=showT(I, cned, ellipses);
	namedWindow("Triplets", CV_WINDOW_AUTOSIZE);
	imshow("Triplets", image1);
	////image1 = Mat(ellipses, true);
	////fromVectorToMat(ellipses, image1);
	//cned.DrawDetectedEllipses(image1, ellipses);//////////////////////////////////////////////////////////////////////////////////////////////
	//namedWindow("finding", CV_WINDOW_AUTOSIZE);
	////imshow("finding", image1);
	//cvSaveImage("C:/Users/zhangruike/Desktop/ell/CNC-Ellipse/image/result/result.jpg",&IplImage(resultImage));
	//imwrite("C:/Users/zhangruike/Desktop/ell/CNC-Ellipse/image/result/033_0022_normal.jpg", resultImage);
	//}


	Toc(2); //grouping	
			// time estimation, validation inside                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
	_times[2] -= (_times[3] + _times[4]);//����׼ȷ

	Tac(4); //validation
			// Sort detected ellipses with respect to score
	sort(ellipses.begin(), ellipses.end());//���շ������� ������ͬ Խ�ӽ�ԲԽС
	Toc(4); //validation


			// Free accumulator memory
	delete[] accN;
	delete[] accR;
	delete[] accA;

	Tic(5);
	// Cluster detections
	ClusterEllipses(ellipses);

	Mat1b image2 = showT(I, cned, ellipses);
	namedWindow("ClusterEllipses", CV_WINDOW_AUTOSIZE);
	//imshow("ClusterEllipses", image2);

	Toc(5);

};



//��Բ������̡��μ����ĵ�[3.3.2]��
// Ellipse clustering procedure. See Sect [3.3.2] in the paper.
//δ���� ����
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


//�����������⵽����Բ
//Draw at most iTopN detected ellipses.
//***showT showT_zrk��OnImage_zrk main_Rect main_allRects case11 case12�� OnVideo(����) OnImage(����)
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

//�����������⵽����Բ zh�ӵĻ��Ҷ�ͼ
//Draw at most iTopN detected ellipses.
void CNEllipseDetector::DrawDetectedEllipses(Mat1b& output, vector<Ellipse>& ellipses, int iTopN, int thickness)
{
	int sz_ell = int(ellipses.size());
	int n = (iTopN == 0) ? sz_ell : min(iTopN, sz_ell);
	for (int i = 0; i < n; ++i)
	{
		Ellipse& e = ellipses[n - i - 1];
		int g = cvRound(e._score * 255.f);//��������
		Scalar color(0, g, 0);

		ellipse(output, Point(cvRound(e._xc), cvRound(e._yc)), Size(cvRound(e._a), cvRound(e._b)), e._rad*180.0 / CV_PI, 0.0, 360.0, color, thickness);// thickness���
	}
}

// add 2015��12��24��  zrk �����޵��á�����
void CNEllipseDetector::showEdgeInPic(Mat1b& I) {
	_szImg = I.size();
	// initialize temporary data structures
	Mat1b DP = Mat1b::zeros(_szImg);		// arcs along positive diagonal
	Mat1b DN = Mat1b::zeros(_szImg);		// arcs along negative diagonal

											//other temporary 
	vector<vector<Point>> points_1, points_2, points_3, points_4;		//vector of points, one for each convexity class

																		//Preprocessing
																		//From input image I, find edge point with coarse convexity along positive (DP) or negative (DN) diagonal


	PrePeocessing(I, DP, DN);

	// show DP
	Mat1b DTMP = Mat1b::ones(_szImg) * 255;
	Mat1b tmp = DTMP - DP;
	Mat1b tmp2 = DTMP - DP;
	/*cvShowImage("DP",&(IplImage)tmp );
	cvSaveImage("DP.jpg",&IplImage(tmp));
	tmp=DTMP-DN;
	cvShowImage("DN",&(IplImage)tmp );
	cvSaveImage("DN.jpg",&IplImage(tmp));*/
	tmp = DP + DN;
	tmp2 = DTMP - tmp;
	cvShowImage("DNP", &(IplImage)tmp2);
	cvSaveImage("DPN.jpg", &IplImage(tmp2));
	//end show DpDn
	float t_fMinOrientedRectSide = _fMinOrientedRectSide;
	_fMinOrientedRectSide = 4;
	//detect edges and find convexities// 4������
	DetectEdges13(DP, points_1, points_3);
	DetectEdges24(DN, points_2, points_4);
	_fMinOrientedRectSide = t_fMinOrientedRectSide;
	//��ʾѡ���ı�
	Mat picture(_szImg, CV_8UC3, Scalar(255, 255, 255));
	IplImage *image = cvCreateImage(_szImg, IPL_DEPTH_8U, 3);
	// Mat picture(_szImg,CV_8UC3);
	/*for(int ih=0;ih<_szImg.height;ih++){
	for(int iw=0;iw<_szImg.width;iw++){
	picture.at<Vec3b>(ih,iw)[0]=I.at<uchar>(ih,iw);
	picture.at<Vec3b>(ih,iw)[1]=I.at<uchar>(ih,iw);
	picture.at<Vec3b>(ih,iw)[2]=I.at<uchar>(ih,iw);

	}
	}*/
	vector<Mat> imgs(4);
	Mat picture1 = picture.clone();
	Mat picture2 = picture.clone();
	Mat picture3 = picture.clone();
	Mat picture4 = picture.clone();
	showEdge(points_1, picture1);
	showEdge(points_2, picture2);
	showEdge(points_3, picture3);
	showEdge(points_4, picture4);
	/*
	imshow("��1����", picture1);
	imshow("��2����", picture2);
	imshow("��3����", picture3);
	imshow("��4����", picture4);*/

	Mat picture5 = picture.clone();
	showEdge(points_1, picture5);
	showEdge(points_2, picture5);
	showEdge(points_3, picture5);
	showEdge(points_4, picture5);

	imgs[0] = picture2;
	imgs[1] = picture1;
	imgs[2] = picture3;
	imgs[3] = picture4;

	cvNamedWindow("all arcs", CV_WINDOW_NORMAL);
	cvShowImage("all arcs", &(IplImage)picture5);
	cvSaveImage("arcs.jpg", &IplImage(picture5));

	//imshow("all arcs", picture5);

	MultiImage_OneWin("�����еĻ�", imgs, cvSize(2, 2), cvSize(400, 800));
	//cvWaitKey(0);  
	//cvDestroyWindow("all arcs");
}
//zrk �޵���
void CNEllipseDetector::showEdgeInPic_zrk(Mat1b& I) {
	_szImg = I.size();
	// initialize temporary data structures
	Mat1b DP = Mat1b::zeros(_szImg);		// arcs along positive diagonal
	Mat1b DN = Mat1b::zeros(_szImg);		// arcs along negative diagonal

											//other temporary 
	vector<vector<Point>> points_1, points_2, points_3, points_4;		//vector of points, one for each convexity class

																		//Preprocessing
																		//From input image I, find edge point with coarse convexity along positive (DP) or negative (DN) diagonal


	PrePeocessing(I, DP, DN);

	// show DP
	Mat1b DTMP = Mat1b::ones(_szImg) * 255;
	Mat1b tmp = DTMP - DP;
	Mat1b tmp2 = DTMP - DP;
	/*cvShowImage("DP",&(IplImage)tmp );
	cvSaveImage("DP.jpg",&IplImage(tmp));
	tmp=DTMP-DN;
	cvShowImage("DN",&(IplImage)tmp );
	cvSaveImage("DN.jpg",&IplImage(tmp));*/
	tmp = DP + DN;
	tmp2 = DTMP - tmp;
	cvShowImage("DNP", &(IplImage)tmp2);
	cvSaveImage("DPN.jpg", &IplImage(tmp2));
	//end show DpDn
	float t_fMinOrientedRectSide = _fMinOrientedRectSide;
	_fMinOrientedRectSide = 4;
	//detect edges and find convexities// 4������
	DetectEdges13_zrk(DP, points_1, points_3);
	DetectEdges24_zrk(DN, points_2, points_4);
	_fMinOrientedRectSide = t_fMinOrientedRectSide;
	//��ʾѡ���ı�
	Mat picture(_szImg, CV_8UC3, Scalar(255, 255, 255));
	IplImage *image = cvCreateImage(_szImg, IPL_DEPTH_8U, 3);
	// Mat picture(_szImg,CV_8UC3);
	/*for(int ih=0;ih<_szImg.height;ih++){
	for(int iw=0;iw<_szImg.width;iw++){
	picture.at<Vec3b>(ih,iw)[0]=I.at<uchar>(ih,iw);
	picture.at<Vec3b>(ih,iw)[1]=I.at<uchar>(ih,iw);
	picture.at<Vec3b>(ih,iw)[2]=I.at<uchar>(ih,iw);

	}
	}*/
	vector<Mat> imgs(4);
	Mat picture1 = picture.clone();
	Mat picture2 = picture.clone();
	Mat picture3 = picture.clone();
	Mat picture4 = picture.clone();
	showEdge(points_1, picture1);
	showEdge(points_2, picture2);
	showEdge(points_3, picture3);
	showEdge(points_4, picture4);
	/*
	imshow("��1����", picture1);
	imshow("��2����", picture2);
	imshow("��3����", picture3);
	imshow("��4����", picture4);*/

	Mat picture5 = picture.clone();
	showEdge(points_1, picture5);
	showEdge(points_2, picture5);
	showEdge(points_3, picture5);
	showEdge(points_4, picture5);

	imgs[0] = picture2;
	imgs[1] = picture1;
	imgs[2] = picture3;
	imgs[3] = picture4;

	cvNamedWindow("all arcs", CV_WINDOW_NORMAL);
	cvShowImage("all arcs", &(IplImage)picture5);
	cvSaveImage("arcs.jpg", &IplImage(picture5));

	imshow("all arcs", picture5);

	MultiImage_OneWin("�����еĻ�", imgs, cvSize(2, 2), cvSize(400, 800));
	//cvWaitKey(0);  
	//cvDestroyWindow("all arcs");
}
//Ϊ��showAllEdgeInPic������//��Ϊ�޵���
void CNEllipseDetector::showAllEdgeInPic(Mat1b& I) {
	_szImg = I.size();
	// initialize temporary data structures
	Mat1b DP = Mat1b::zeros(_szImg);		// arcs along positive diagonal
	Mat1b DN = Mat1b::zeros(_szImg);		// arcs along negative diagonal

											//other temporary 
	vector<vector<Point>> points_1, points_2, points_3, points_4;		//vector of points, one for each convexity class

																		//Preprocessing
																		//From input image I, find edge point with coarse convexity along positive (DP) or negative (DN) diagonal

																		//this->_fMaxRectAxesRatio=0;
																		//this->_iMinEdgeLength=0;
																		//this->_uNs=0;

	PrePeocessing(I, DP, DN);
	cvNamedWindow("aaas", CV_WINDOW_NORMAL);
	Mat1b tmp = Mat1b::zeros(_szImg);


	cvShowImage("aaas", &(IplImage)DP);
	//detect edges and find convexities// 4������
	//DetectEdges13(DP, points_1, points_3);
	//DetectEdges24(DN, points_2, points_4);
}
//***detect4s database_4s  case 5 case8  ����
int CNEllipseDetector::showEdgeInPic(Mat1b& I, bool showedge) {
	_szImg = I.size();
	// initialize temporary data structures
	Mat1b DP = Mat1b::zeros(_szImg);		// arcs along positive diagonal
	Mat1b DN = Mat1b::zeros(_szImg);		// arcs along negative diagonal

											//other temporary 
	vector<vector<Point>> points_1, points_2, points_3, points_4;		//vector of points, one for each convexity class

																		//Preprocessing
																		//From input image I, find edge point with coarse convexity along positive (DP) or negative (DN) diagonal


	PrePeocessing(I, DP, DN);
	//detect edges and find convexities// 4������
	DetectEdges13(DP, points_1, points_3);
	DetectEdges24(DN, points_2, points_4);
	int EdgesNumber = points_1.size() + points_2.size() + points_3.size() + points_4.size();
	//��ʾѡ���ı�
	return EdgesNumber;
	//cvWaitKey(0);  
	//cvDestroyWindow("all arcs");
}
//zrk �޵���
int CNEllipseDetector::showEdgeInPic_zrk(Mat1b& I, bool showedge) {
	_szImg = I.size();
	// initialize temporary data structures
	Mat1b DP = Mat1b::zeros(_szImg);		// arcs along positive diagonal
	Mat1b DN = Mat1b::zeros(_szImg);		// arcs along negative diagonal

											//other temporary 
	vector<vector<Point>> points_1, points_2, points_3, points_4;		//vector of points, one for each convexity class

																		//Preprocessing
																		//From input image I, find edge point with coarse convexity along positive (DP) or negative (DN) diagonal


	PrePeocessing(I, DP, DN);
	//detect edges and find convexities// 4������
	DetectEdges13_zrk(DP, points_1, points_3);
	DetectEdges24_zrk(DN, points_2, points_4);
	int EdgesNumber = points_1.size() + points_2.size() + points_3.size() + points_4.size();

	//��ʾѡ���ı�
	return EdgesNumber;
	//cvWaitKey(0);  
	//cvDestroyWindow("all arcs");
}