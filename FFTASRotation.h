/**************************************************************************************
* version 1.0 2019/2/26
* Designer 金铭
* E-mail jinmingaps@163.com
* Fuction：
* 坐标旋转：
* 输入 Eu0, Ev0，得到坐标旋转后的Eut,Evt;
* 需要输入的参数，频率，目标距离，计算点数N*M;
目标面倾斜角（theta和phi）、ds（默认波长/3.5）
***************************************************************************************/
#pragma once
#ifndef FFTASRotation_H
#define FFTASRotation_H

#include <cmath>
#include <complex>
#include <vector>


#include "FFT.h"
class GraphTrans;
class _declspec(dllexport) FFTASRotation
{
public:
	FFTASRotation(double _f = 10.65e9, int _Nu = 361, int _Nv = 361);
	~FFTASRotation();
	void SetParas(double _f, double _Nu, double _Nv, double _ds);
	void SetParas(double _f, double _Nu, double _Nv, double _du,double dv);
	void SetRotationParas(GraphTrans GT0, GraphTrans GTt);
	void Allocate();
	void FreeCal();
	//输出Ex Ey Ez并调用FreeCal
	void output(std::complex<double> ** EuOut, std::complex<double> ** EvOut);
	void output(std::complex<double> ** EuOut, std::complex<double> ** EvOut, std::complex<double> ** EnOut);
	//设置输入并补0
	void SetInput(std::complex<double> ** EuIn, std::complex<double> ** EvIn);
	void PerformRotate();
	//输出Ex Ey Ez Hx Hy Hz 并调用FreeCal

private:

private:
	double freq; // 频率 默认 10.65e9
	double lambda; // 波长
	double k;
	double du,dv;
	double theta_rou, phi_rou;//旋转
	int Nu, Nv; //计算的点数 一般设为奇数 默认361
	int Nu2, Nv2; // N2 = 2 * N -1

	std::complex<double> ** Eu0, ** Ev0, ** En0;
	std::complex<double> ** Eut, ** Evt, ** Ent;

};



#endif // !CALCUlATION_H
