opencv简单实践6：canny算子、 旋转及缩放、仿射和投影变换


练习1 canny算子
读取摄像头图像，并对摄像头图像进行中值滤波
0. x方向梯度信息，(CV_16SC1 or CV_16SC3)。1. y方向梯度信息，(CV_16SC1 or CV_16SC3)。2. 输出，边缘检测结果3. double类型的threshold1，第一个滞后性阈值。4. double类型的threshold2，第二个滞后性阈值。5. bool类型的L2gradient，一个计算图像梯度幅值的标识，有默认值false。
0. 输入图像，即源图像，需为单通道8位图像。1. 输出，边缘检测结果2. double类型的threshold1，第一个滞后性阈值。3. double类型的threshold2，第二个滞后性阈值。4. 表示应用Sobel算子的核大小，默认值3。5. bool类型的L2gradient，一个计算图像梯度幅值的标识，有默认值false。
重载函数是函数的一种特殊情况，为方便使用，C++允许在同一范围中声明几个功能类似的同名函数，但是这些同名函数的形式参数（指参数的个数、类型或者顺序）必须不同，也就是说用同一个函数完成不同的功能。这就是重载函数。重载函数常用来实现功能类似而所处理的数据类型不同的问题。不能只有函数返回值类型不同。


练习2 旋转及缩放
根据角度及缩放比例生成仿射变换矩阵
仿射变换函数


练习3 仿射变换
变换前的三点坐标 ， 变换后的三点坐标，  根据角度及缩放比例，生成仿射变换矩阵，仿射变换函数


练习4 投影变换
变换前的四点坐标，变换后的四点坐标，根据角度及缩放比例，生成仿射变换矩阵，仿射变换函数
仿射矩阵T1：比例、旋转、对称、错切T2：平移T3：投影T4：整体缩放