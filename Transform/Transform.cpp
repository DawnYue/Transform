#include <opencv2/opencv.hpp>  
#include <iostream>
#include<cassert>
#include<opencv2/core/matx.hpp>
#include<vector>

//练习3
using namespace cv;
using namespace std;

int main()
{
	Mat srcImage, dstImage;
	srcImage = imread("E:\\lena.jpg");
	if (!srcImage.data) {  return -1; }	
	float angle = -10.0, scale = 1;	
	cv::Point2f src_pt[] = {//变换前的三点坐标
		cv::Point2f(200,200),
		cv::Point2f(250,200),
		cv::Point2f(200,100)

	};
	cv::Point2f dst_pt[] = {//变换后的三点坐标
		cv::Point2f(300,100),
		cv::Point2f(300,50),
		cv::Point2f(200,100)

	};
	// 仿射变换矩阵
	const cv::Mat affine_matrix = cv::getAffineTransform(src_pt, dst_pt);
	// 仿射变换函数
	cv::warpAffine(srcImage, dstImage, affine_matrix, srcImage.size());

	imshow("srcImage", srcImage);
	imshow("dstImage", dstImage);
	waitKey(0);

}
