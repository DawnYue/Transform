#include <opencv2/opencv.hpp>  
#include <iostream>
#include<cassert>
#include<opencv2/core/matx.hpp>
#include<vector>

//练习4
using namespace cv;
using namespace std;

int main()
{
	Mat srcImage, dstImage;
	srcImage = imread("E:\\lena.jpg");
	if (!srcImage.data) {  return -1; }	
	float angle = -10.0, scale = 1;	
	cv::Point2f pts1[] = {//变换前的四点坐标
		cv::Point2f(150,150),
		cv::Point2f(150,300),
		cv::Point2f(350,300),
		cv::Point2f(350,150)

	};
	cv::Point2f pts2[] = {//变换后的四点坐标
		cv::Point2f(200,150),
		cv::Point2f(200,300),
		cv::Point2f(340,270),
		cv::Point2f(340,180)

	};
	// 仿射变换矩阵
	const cv::Mat perspective_matrix = cv::getPerspectiveTransform(pts1, pts2);
	// 仿射变换函数
	cv::warpPerspective(srcImage, dstImage, perspective_matrix, srcImage.size());

	imshow("srcImage", srcImage);
	imshow("dstImage", dstImage);
	waitKey(0);

}
