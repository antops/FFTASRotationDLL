/**************************************************************************************
* version 1.0 2019/2/26
* Designer ����
* E-mail jinmingaps@163.com
* Fuction��
* ������ת��
* ���� Eu0, Ev0���õ�������ת���Eut,Evt;
* ��Ҫ����Ĳ�����Ƶ�ʣ�Ŀ����룬�������N*M;
Ŀ������б�ǣ�theta��phi����ds��Ĭ�ϲ���/3.5��
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
	//���Ex Ey Ez������FreeCal
	void output(std::complex<double> ** EuOut, std::complex<double> ** EvOut);
	void output(std::complex<double> ** EuOut, std::complex<double> ** EvOut, std::complex<double> ** EnOut);
	//�������벢��0
	void SetInput(std::complex<double> ** EuIn, std::complex<double> ** EvIn);
	void PerformRotate();
	//���Ex Ey Ez Hx Hy Hz ������FreeCal

private:

private:
	double freq; // Ƶ�� Ĭ�� 10.65e9
	double lambda; // ����
	double k;
	double du,dv;
	double theta_rou, phi_rou;//��ת
	int Nu, Nv; //����ĵ��� һ����Ϊ���� Ĭ��361
	int Nu2, Nv2; // N2 = 2 * N -1

	std::complex<double> ** Eu0, ** Ev0, ** En0;
	std::complex<double> ** Eut, ** Evt, ** Ent;

};



#endif // !CALCUlATION_H
