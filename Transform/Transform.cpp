#include <opencv2/opencv.hpp>  
#include <iostream>
#include<cassert>
#include<opencv2/core/matx.hpp>
#include<vector>


using namespace cv;
using namespace std;

#define WINDOW_NAME1 "【原始图窗口】"					//为窗口标题定义的宏 
#define WINDOW_NAME2 "【经过Warp后的图像】"        //为窗口标题定义的宏 
#define WINDOW_NAME3 "【经过Warp和Rotate后的图像】"        //为窗口标题定义的宏 
//Transform
int main()
{


	//【1】参数准备
	//定义两组点，代表两个三角形
	Point2f srcTriangle[3];
	Point2f dstTriangle[3];
	//定义一些Mat变量
	Mat rotMat(2, 3, CV_32FC1);
	Mat warpMat(2, 3, CV_32FC1);
	Mat srcImage, dstImage_warp, dstImage_warp_rotate;

	//【2】加载源图像并作一些初始化
	srcImage = imread("E:\\1.png");
	if (!srcImage.data) { printf("读取图片错误，请确定目录下是否有imread函数指定的图片存在~！ \n"); return false; }
	// 设置目标图像的大小和类型与源图像一致
	dstImage_warp = Mat::zeros(srcImage.rows, srcImage.cols, srcImage.type());

	//【3】设置源图像和目标图像上的三组点以计算仿射变换
	srcTriangle[0] = Point2f(0, 0);
	srcTriangle[1] = Point2f(static_cast<float>(srcImage.cols - 1), 0);
	srcTriangle[2] = Point2f(0, static_cast<float>(srcImage.rows - 1));

	dstTriangle[0] = Point2f(static_cast<float>(srcImage.cols*0.0), static_cast<float>(srcImage.rows*0.33));
	dstTriangle[1] = Point2f(static_cast<float>(srcImage.cols*0.65), static_cast<float>(srcImage.rows*0.35));
	dstTriangle[2] = Point2f(static_cast<float>(srcImage.cols*0.15), static_cast<float>(srcImage.rows*0.6));

	//【4】求得仿射变换
	warpMat = getAffineTransform(srcTriangle, dstTriangle);

	//【5】对源图像应用刚刚求得的仿射变换
	warpAffine(srcImage, dstImage_warp, warpMat, dstImage_warp.size());

	//【6】对图像进行缩放后再旋转
	// 计算绕图像中点顺时针旋转50度缩放因子为0.6的旋转矩阵
	Point center = Point(dstImage_warp.cols / 2, dstImage_warp.rows / 2);
	double angle = -50.0;
	double scale = 0.6;
	// 通过上面的旋转细节信息求得旋转矩阵
	rotMat = getRotationMatrix2D(center, angle, scale);
	// 旋转已缩放后的图像
	warpAffine(dstImage_warp, dstImage_warp_rotate, rotMat, dstImage_warp.size());

	//【7】显示结果
	imshow(WINDOW_NAME1, srcImage);
	imshow(WINDOW_NAME2, dstImage_warp);
	imshow(WINDOW_NAME3, dstImage_warp_rotate);

	// 等待用户按任意按键退出程序
	waitKey(0);

	return 0;
}

//最近邻域插值
//根据目标图像的像素点（浮点坐标）找到原始图像中的4个像素点，取距离该像素点最近的一个原始像素值作为该点的值。


namespace mycv {
	void nearestIntertoplation(cv::Mat& src, cv::Mat& dst, const int rows, const int cols);
}//mycv

double distance(const double x1, const double y1, const double x2, const double y2);//两点之间距离，这里用欧式距离

void mycv::nearestIntertoplation(cv::Mat& src, cv::Mat& dst, const int rows, const int cols)
{
	//比例尺
	const double scale_row = static_cast<double>(src.rows) / rows;
	const double scale_col = static_cast<double>(src.rows) / cols;

	//扩展src到dst
	dst = cv::Mat(rows, cols, src.type());
	assert(src.channels() == 1 && dst.channels() == 1);

	for (int i = 0; i < rows; ++i)//dst的行
		for (int j = 0; j < cols; ++j)//dst的列
		{
			//求插值的四个点
			double y = (i + 0.5) * scale_row + 0.5;
			double x = (j + 0.5) * scale_col + 0.5;
			int x1 = static_cast<int>(x);//col对应x
			if (x1 >= (src.cols - 2)) x1 = src.cols - 2;//防止越界
			int x2 = x1 + 1;
			int y1 = static_cast<int>(y);//row对应y
			if (y1 >= (src.rows - 2))  y1 = src.rows - 2;
			int y2 = y1 + 1;
			//根据目标图像的像素点（浮点坐标）找到原始图像中的4个像素点，取距离该像素点最近的一个原始像素值作为该点的值。
			assert(0 < x2 && x2 < src.cols && 0 < y2 &&  y2 < src.rows);
			std::vector<double> dist(4);
			dist[0] = distance(x, y, x1, y1);
			dist[1] = distance(x, y, x2, y1);
			dist[2] = distance(x, y, x1, y2);
			dist[3] = distance(x, y, x2, y2);

			int min_val = dist[0];
			int min_index = 0;
			for (int i = 1; i < dist.size(); ++i)
				if (min_val > dist[i])
				{
					min_val = dist[i];
					min_index = i;
				}

			switch (min_index)
			{
			case 0:
				dst.at<uchar>(i, j) = src.at<uchar>(y1, x1);
				break;
			case 1:
				dst.at<uchar>(i, j) = src.at<uchar>(y1, x2);
				break;
			case 2:
				dst.at<uchar>(i, j) = src.at<uchar>(y2, x1);
				break;
			case 3:
				dst.at<uchar>(i, j) = src.at<uchar>(y2, x2);
				break;
			default:
				assert(false);
			}
		}
}

double distance(const double x1, const double y1, const double x2, const double y2)//两点之间距离，这里用欧式距离
{
	return (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2);//只需比较大小，返回距离平方即可
}


//线性内插（双线性插值)


namespace mycv {
	void bilinearIntertpolatioin(cv::Mat& src, cv::Mat& dst, const int rows, const int cols);
}//mycv


void mycv::bilinearIntertpolatioin(cv::Mat& src, cv::Mat& dst, const int rows, const int cols)
{
	//比例尺
	const double scale_row = static_cast<double>(src.rows) / rows;
	const double scale_col = static_cast<double>(src.rows) / cols;

	//扩展src到dst
	dst = cv::Mat(rows, cols, src.type());
	assert(src.channels() == 1 && dst.channels() == 1);

	for (int i = 0; i < rows; ++i)//dst的行
		for (int j = 0; j < cols; ++j)//dst的列
		{
			//求插值的四个点
			double y = (i + 0.5) * scale_row + 0.5;
			double x = (j + 0.5) * scale_col + 0.5;
			int x1 = static_cast<int>(x);//col对应x
			if (x1 >= (src.cols - 2)) x1 = src.cols - 2;//防止越界
			int x2 = x1 + 1;
			int y1 = static_cast<int>(y);//row对应y
			if (y1 >= (src.rows - 2))  y1 = src.rows - 2;
			int y2 = y1 + 1;

			assert(0 < x2 && x2 < src.cols && 0 < y2 &&  y2 < src.rows);
			//插值公式,参考维基百科矩阵相乘的公式https://zh.wikipedia.org/wiki/%E5%8F%8C%E7%BA%BF%E6%80%A7%E6%8F%92%E5%80%BC

			cv::Matx12d matx = { x2 - x, x - x1 };
			cv::Matx22d matf = { static_cast<double>(src.at<uchar>(y1, x1)), static_cast<double>(src.at<uchar>(y2, x1)),
								 static_cast<double>(src.at<uchar>(y1, x2)), static_cast<double>(src.at<uchar>(y2, x2)) };
			cv::Matx21d maty = {
				y2 - y,
				y - y1
			};

			auto  val = (matx * matf * maty);
			dst.at<uchar>(i, j) = val(0, 0);
		}

}


//双三次插值
namespace mycv {
	void bicubicInsterpolation(cv::Mat& src, cv::Mat& dst, const int rows, const int cols);
	std::vector<double> getW(double coor, double a = -0.5);//a默认-0.5
}//mycv



void mycv::bicubicInsterpolation(cv::Mat& src, cv::Mat& dst, const int rows, const int cols)
{
	dst = cv::Mat(rows, cols, src.type(), cv::Scalar::all(0));//初始化dst
	//比例尺
	double row_scale = static_cast<double>(src.rows) / rows;
	double col_scale = static_cast<double>(src.cols) / cols;

	switch (src.channels())
	{

	case 1://灰度
		for (int i = 2; i < dst.rows - 2; ++i)
			for (int j = 2; j < dst.cols - 2; ++j)
			{
				//计算系数w
				double r = static_cast<double>(i *  row_scale);
				double c = static_cast<double>(j *  col_scale);

				//防止越界
				if (r < 1.0) r += 1.0;
				if (c < 1.0) c += 1.0;


				std::vector<double> wr = mycv::getW(r);
				std::vector<double> wc = mycv::getW(c);

				//最后计算插值得到的灰度值
				double val = 0;
				int cc = static_cast<int>(c);
				int rr = static_cast<int>(r);
				//防止越界
				if (cc > src.cols - 3)
				{
					cc = src.cols - 3;

				}
				if (rr > src.rows - 3) rr = src.rows - 3;

				assert(0 <= (rr - 1) && 0 <= (cc - 1) && (rr + 2) < src.rows && (cc + 2) < src.cols);

				//4x4数量的点(rr, cc) -> (y, x)
				std::vector<std::vector<int> > src_arr = {
					{ src.at<uchar>(rr - 1, cc - 1), src.at<uchar>(rr, cc - 1), src.at<uchar>(rr + 1, cc - 1), src.at<uchar>(rr + 2, cc - 1)},
					{ src.at<uchar>(rr - 1, cc), src.at<uchar>(rr, cc), src.at<uchar>(rr + 1, cc), src.at<uchar>(rr + 2, cc)},
					{ src.at<uchar>(rr - 1, cc + 1), src.at<uchar>(rr, cc + 1), src.at<uchar>(rr + 1, cc + 1), src.at<uchar>(rr + 2, cc + 1)},
					{ src.at<uchar>(rr - 1, cc + 2), src.at<uchar>(rr, cc + 2), src.at<uchar>(rr + 1, cc + 2), src.at<uchar>(rr + 2, cc + 2)}
				};
				for (int p = 0; p < 3; ++p)
					for (int q = 0; q < 3; ++q)
					{
						//val(p, q) = w(p,q) * src(p, q)
						val += wr[p] * wc[q] * static_cast<double>(src_arr[p][q]);
					}
				assert(i < dst.rows && j < dst.cols);
				dst.at<uchar>(i, j) = static_cast<int>(val);

			}
		break;
	case 3://彩色(原理一样多了两个通道而已)
		break;
	default:
		break;
	}
}

std::vector<double> mycv::getW(double coor, double a)
{
	std::vector<double> w(4);
	int base = static_cast<int>(coor);//取整作为基准
	double e = coor - static_cast<double>(base);//多出基准的小数部分
	std::vector<double> tmp(4);//存放公式中 |x| <= 1,  1 < |x| < 2四个值
	// 4 x 4的16个点，所以tmp[0]和tmp[4]距离较远，值在[1, 2]区间
	tmp[0] = 1.0 + e;// 1 < x < 2
	tmp[1] = e;//x <= 1
	tmp[2] = 1.0 - e;// x <= 1
	tmp[3] = 2.0 - e;// 1 < x < 2

	//按照bicubic公式计算系数w
	// x <= 1
	w[1] = (a + 2.0) * std::abs(std::pow(tmp[1], 3)) - (a + 3.0) * std::abs(std::pow(tmp[1], 2)) + 1;
	w[2] = (a + 2.0) * std::abs(std::pow(tmp[2], 3)) - (a + 3.0) * std::abs(std::pow(tmp[2], 2)) + 1;
	// 1 < x < 2
	w[0] = a * std::abs(std::pow(tmp[0], 3)) - 5.0 * a * std::abs(std::pow(tmp[0], 2)) + 8.0*a*std::abs(tmp[0]) - 4.0*a;
	w[3] = a * std::abs(std::pow(tmp[3], 3)) - 5.0 * a * std::abs(std::pow(tmp[3], 2)) + 8.0*a*std::abs(tmp[3]) - 4.0*a;

	return w;
}