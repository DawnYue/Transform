#include<opencv2/opencv.hpp>
#include<iostream>
using namespace cv;
using namespace std;

//练习1
int main(int argc, char*argv)
{
	//实例化一个videocapture类，名称为cap
	VideoCapture cap;
	//cap(0)表示打开本机的第一个摄像头
	cap.open(0);
	//isOpened()检查视频是否开启，正常开启返回1	
	if (!cap.isOpened())
	{
		std::cout << "不能打开视频文件" << std::endl;
		return -1;
	}
	//获得视频的fps
	double fps = cap.get(CAP_PROP_FPS);
	std::cout << "fps" << fps << std::endl;

	while (1)
	{
		cv::Mat frame;
		cv::Mat dx, dy,canny2, canny1;
		bool rSucess = cap.read(frame);

		if (!rSucess)
		{
			std::cout << "不能从视频中读取帧" << std::endl;
			break;
		}
		else
		{
			cvtColor(frame, frame, COLOR_BGR2GRAY);// 将原图像转换为灰度图像
			threshold(frame, frame, 0, 255, CV_THRESH_OTSU);//二值化
			convertScaleAbs(frame, frame);
			Sobel(frame, dx, CV_16SC1, 2, 0, 9);//X方向
			Sobel(frame, dy, CV_16SC1, 0, 2, 9);//Y方向
			Canny(frame, canny2, 20, 60);
			Canny(dx, dy, canny1, 20, 60);

			imshow("canny2边缘检测", canny2);
			imshow("canny1边缘检测", canny1);
		}
		waitKey(30); //延时30ms

	}
	return 0;
}
