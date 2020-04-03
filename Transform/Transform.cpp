#include <opencv2/opencv.hpp>  
#include <iostream>
#include<cassert>
#include<opencv2/core/matx.hpp>
#include<vector>

//练习2
using namespace cv;
using namespace std;

int main()
{
	Mat srcImage, dstImage;
	srcImage = imread("E:\\lena.jpg");
	if (!srcImage.data) { printf("读取图片错误 \n"); return -1; }	
	float angle = -10.0, scale = 1;	
	cv::Point2f center(srcImage.cols / 2, srcImage.rows / 2);
	// 变换矩阵
	const cv::Mat affine_matrix = cv::getRotationMatrix2D(center, angle, scale);
	// 仿射变换函数
	cv::warpAffine(srcImage, dstImage, affine_matrix, srcImage.size());

	imshow("srcImage", srcImage);
	imshow("dstImage", dstImage);
	waitKey(0);

}
