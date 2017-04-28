#include "expcode.h"
#include <iostream>]
//#include<math.h>
//#include<stdlib.h>
//#include"math.h"

//#define pi 3.14159
//int MAX_MV = int(256 * sqrt(2.0) + 2);
//using namespace std;

#define N 256
#define M 8
#define PI 3.1415926
using namespace std;

void reverse(float* rdata, float* idata,int RI){
	int i,j;
	int k;
	float temp;
	j=N/2;
	for(i=1;i<N-2;i++){
		if(i<j){
			temp=rdata[i];
			rdata[i]=rdata[j];
			rdata[j]=temp;
			temp=idata[i];
			if(RI){
				temp=idata[i];
				idata[i]=idata[j];
				idata[j]=temp;
			}

		}
		k=N/2;
		while(1){
			if(j<k){
				j=j+k;
				break;
			}
			else{
				j=j-k;
				k=k/2;

			}
		}
	}
}

void fft(float* rdata,float* idata){
	float Tr,Ti;
	float temp;
	int p,J,L,B;
	int i,j,k;
	for(L=1;L<=M;L++){
		B=1;
		i=L-1;
		while(i){
			B*=2;
			i--;
		}
		for(J=0;J<=B-1;J++){
			p=1;
			i=M-L;
			while(i){
				p=p*2;
				i--;
			}
			p=J*p;
			for(k=J;k<=N-1;k=k+2*B){
				Tr=rdata[k];
				Ti=idata[k];
				temp=rdata[k+B];
				rdata[k]=rdata[k]+rdata[k+B]*cos(2.0*PI*p/N)+idata[k+B]*sin(2.0*PI*p/N);
				idata[k]=idata[k]+idata[k+B]*cos(2.0*PI*p/N)-rdata[k+B]*sin(2.0*PI*p/N);
				rdata[k+B]=Tr-rdata[k+B]*cos(2.0*PI*p/N)-idata[k+B]*sin(2.0*PI*p/N);
				idata[k+B]=Ti-idata[k+B]*cos(2.0*PI*p/N)+temp*sin(2.0*PI*p/N);
			}
		}
	}
}

void ifft(float* rdata,float* idata){
	float Tr,Ti;
	float temp;
	int p,J,L,B;
	int i,j,k;
	for(L=1;L<=M;L++){
		B=1;
		i=L-1;
		while(i){
			B*=2;
			i--;
		}
		for(J=0;J<=B-1;J++){
			p=-1;
			i=M-L;
			while(i){
				p=p*2;
				i--;
			}
			p=J*p;
			for(k=J;k<=N-1;k=k+2*B){
				Tr=rdata[k];
				Ti=idata[k];
				temp=rdata[k+B];
				rdata[k]=rdata[k]+rdata[k+B]*cos(2.0*PI*p/N)+idata[k+B]*sin(2.0*PI*p/N);
				idata[k]=idata[k]+idata[k+B]*cos(2.0*PI*p/N)-rdata[k+B]*sin(2.0*PI*p/N);
				rdata[k+B]=Tr-rdata[k+B]*cos(2.0*PI*p/N)-idata[k+B]*sin(2.0*PI*p/N);
				idata[k+B]=Ti-idata[k+B]*cos(2.0*PI*p/N)+temp*sin(2.0*PI*p/N);
			}
		}
	}
}


void idft(float** rdct,float** idct){
	float* rcolumn=new float[N];
	float* icolumn=new float[N];
	int i,j;

	for(i=0;i<N;i++){
		reverse(rdct[i],idct[i],1);
		ifft(rdct[i],idct[i]);
	}
	for(j=0;j<N;j++){
		for(i=0;i<N;i++){
			rcolumn[i]=rdct[i][j]*1.0/N;
			icolumn[i]=idct[i][j]*1.0/N;
		}
		reverse(rcolumn,icolumn,1);
		ifft(rcolumn,icolumn);
		for(i=0;i<N;i++){
			rdct[i][j]=rcolumn[i];
			idct[i][j]=icolumn[i];
		}
	}

}



/**********************************************
**********自定义函数****************************
***********************************************/
float gaussrand()
{
	float v1, v2, s;
	int phase = 0;
	float x;

	if(phase == 0)
	{
		do
		{
			float u1 = (float)rand() / RAND_MAX;
			float u2 = (float)rand() / RAND_MAX;

			v1 = 2 * u1 - 1;
			v2 = 2 * u2 - 1;
			s = v1 * v1 + v2 * v2;
		}while(s >= 1 || s == 0);

		x = v1 * sqrt(-2 * log(s) / s);
	}
	else
		x= v2 * sqrt(-2 *log(s) / s);
	phase = 1 - phase;
	return x;
}

int** img_padding(int** pixelmat, int mheight = 256, int mwidth = 256, int ker = 3)
{
	int halfk = ker/2;
	int h = mheight + 2*halfk, w = mwidth + 2*halfk;	//the new m,h
	int** img_pad = new int*[h];
	int i, j;

	for(i = 0; i < h; i++)
		img_pad[i] = new int[w];
	for(j = halfk; j < mwidth + halfk; j++)
	{
		for(i = 0; i < halfk; i++)
		{
			img_pad[halfk - i -1][j] = pixelmat[i][j - halfk];
			img_pad[mheight + halfk + i][j] = pixelmat[mheight - 1 - i][j - halfk];
		}
		for(i = halfk; i < mheight + halfk; i++)
			img_pad[i][j] = pixelmat[i - halfk][j - halfk];
	}
	for(i = 0; i < h; i++)
	{
		for(j = 0; j < halfk; j++)
		{
			img_pad[i][halfk - 1 - j] = img_pad[i][halfk + j];
			img_pad[i][mwidth + halfk + j] = img_pad[i][mwidth + halfk - j - 1];
		}
	}
	return img_pad;
}

int getmedinum(int arr[], int ker)
{
	int i, j, tmp, medi = 0;

	//选择排序
	for(i = 0; i < ker * ker; i++)
		for(j = i + 1; j < 9; j++)
			if(arr[i] > arr[j])
			{
				tmp = arr[i];
				arr[i] = arr[j];
				arr[j] = tmp;
			}
			return arr[ker * ker / 2];
}

int** kernal_res(int** img_tmp, int** pixelmat, int mheight, int mwidth, int kernal[], int ker, int flag)
{
	int i, j, k, l, m;
	int res;
	int halfk = ker/2;
	for(i = halfk; i < mheight + halfk; i++)
		for(j = halfk; j < mwidth + halfk; j++)
		{
			res = 0;
			m = 0;
			for(k = 0 - halfk; k <= halfk; k++)
				for(l = 0 - halfk; l <= halfk; l++)
					res  += img_tmp[i + k][j + l] * kernal[m++];
			if(flag)
			{
				res = res > 255 ? 255 : res;
				res = res < 0 ? 0 : res;
			}
			pixelmat[i - 1][j -1] = res;   //update the pixel value
			//img_tmp[i][j] = res;
		}
		return pixelmat;
}



//示例: 求图像中心点像素的灰度值
int midptvalue(int** pixelmat, int mheight, int mwidth)
{
	//pixelmat为指向图像像素值二维数组指针, 指向图像左上角的像素, 像素值类型为int;
	//mheight为图像的高度(行), mwidth为图像的宽度(列);
	int middlerow = mheight / 2 - 1;
	int middlecol = mwidth / 2 - 1;
	return pixelmat[middlerow][middlecol];
}

//求最大值
int maxvalue(int** pixelmat, int mheight, int mwidth)
{
	int max = 0, i, j;
	for(i = 0; i < mheight; i++)
		for(j = 0; j < mwidth; j++)
			max = max < pixelmat[i][j] ? pixelmat[i][j] : max;
	return max;
}

//求最小值
int minvalue(int** pixelmat, int mheight, int mwidth)
{
	int min = 255, i, j;
	for(i = 0; i < mheight; i++)
		for(j = 0; j < mwidth; j++)
			min = min > pixelmat[i][j] ? pixelmat[i][j] : min;
	return min;
}

//求平均值
float avgvalue(int** pixelmat, int mheight, int mwidth)
{
	int sum = 0, i, j;
	for(i = 0; i < mheight; i++)
		for(j = 0; j < mwidth; j++)
			sum += pixelmat[i][j];
	return (float)sum / (mheight * mwidth);
}

//求方差
float varvalue(int** pixelmat, int mheight, int mwidth)
{
	int i, j;
	float avg, var;
	avg = avgvalue(pixelmat, mheight, mwidth);
	for(i = 0; i < mheight; i++)
		for(j = 0; j < mwidth; j++)
			var += pow(pixelmat[i][j] - avg, 2);
	return var/(mheight*mwidth);
}

//统计直方图, 返回长度为256的1维数组
int* histogram(int** pixelmat, int mheight, int mwidth)
{
	int* hg = new int[256];
	int i , j;
	for(i = 0; i < 256; i++)
		hg[i] = 0;
	for(i = 0; i < mheight; i++)
		for(j = 0; j < mwidth; j++)
			hg[pixelmat[i][j]] ++;
	return hg;
}

//示例,将灰度图转化为二值图像,返回处理后的图像
int** binaryimg(int** pixelmat, int mheight, int mwidth)
{
	for(int i = 0; i < mheight; i++)
	{
		for(int j = 0; j < mwidth; j++)
		{
			//从左上角开始遍历整幅图像, 实现二值化;
			pixelmat[i][j] = pixelmat[i][j] > 128 ? 255 : 0;
		}
	}
	//根据实验要求返回对应的值;
	return pixelmat;
}

//直方图均衡, 返回处理后的图像
int** histogramequ(int** pixelmat, int mheight, int mwidth)
{
	int* hg = new int[256];
	int i, j, pixel_sum = mheight * mwidth;
	hg = histogram(pixelmat, mheight, mwidth);
	for(i = 1; i < 256; i++)
		hg[i] += hg[i - 1];
	for(i = 0; i < 256; i++)
		hg[i] = floor(255 * ((float)hg[i] / pixel_sum) + 0.5);
	for(i = 0; i < mheight; i++)
		for(j = 0; j < mwidth; j++)
			pixelmat[i][j] = hg[pixelmat[i][j]];
	return pixelmat;
}

//灰度拉伸, 返回处理后的图像
int** graystretch(int** pixelmat, int mheight, int mwidth)
{
	int min = minvalue(pixelmat, mheight, mwidth);
	int max = maxvalue(pixelmat, mheight, mwidth);

	int i, j;

	int s1 = min, r1 = 0;
	int s2 = max, r2 = 255;
	for(i = 0; i < mheight; i++)
		for(j = 0; j < mwidth; j++)
			pixelmat[i][j] =int((float)((pixelmat[i][j] - s1) * (r2 - r1)) / (s2 - s1));
	return pixelmat;
}


//中值滤波, 返回处理后的图像
int** medianfit(int** pixelmat, int mheight, int mwidth)
{
	int i, j, k, l, m;		//循环计数用
	int ker = 3, halfk = ker/2;	//模板核size = ker X ker
	int medi;	//中值
	int* arr = new int[ker * ker];	//	模板计算后得到的数组
	int** img_tmp = new int*[mheight + halfk + 1];
	for(i = 0; i < mheight; i++)
		img_tmp[i] = new int[mwidth + halfk + 1];
	img_tmp = img_padding(pixelmat, mheight, mwidth, ker);   //镜像补充

	for(i = halfk; i < mheight + halfk; i++)
		for(j = halfk; j < mwidth + halfk; j++)
		{
			m = 0;
			for(k = 0 - halfk; k <= halfk; k++)
				for(l = 0 - halfk; l <= halfk; l++)
					arr[m++] = img_tmp[i + k][j + l];
			medi = getmedinum(arr, ker);
			pixelmat[i - 1][j -1] = medi;   //update the pixel value
			img_tmp[i][j] = medi;
		}
		return pixelmat;
}


//均值滤波, 返回处理后的图像
int** averagefit(int** pixelmat, int mheight, int mwidth)
{
	int kernal[9] = {1,1,1,
		1,1,1,
		1,1,1};
	int i, j;		//循环计数用
	int ker = 3, halfk = ker/2;	//模板核size = ker X ker
	int** img_tmp = new int*[mheight + halfk + 1];

	for(i = 0; i < mheight; i++)
		img_tmp[i] = new int[mwidth + halfk + 1];

	img_tmp = img_padding(pixelmat, mheight, mwidth, ker);   //镜像补充

	pixelmat = kernal_res(img_tmp, pixelmat, mheight, mwidth, kernal, ker, 0);
	for(i = 0; i <mheight; i++)
		for(j = 0; j < mwidth; j++)
			pixelmat[i][j] = pixelmat[i][j] / 9;
	return pixelmat;
}



//理想低通滤波, 返回处理后的图像

int** lowpassfit(int** pixelmat, int mheight, int mwidth)
{
	float** rdct=new float*[mheight];
	float** idct=new float*[mheight];
	float* rcolumn=new float[N];
	float* icolumn=new float[N];
	int i,j;
	for(i=0;i<N;i++){
		idct[i]=new float[N];
		rdct[i]=new float[N];
		for(j=0;j<N;j++){
			idct[i][j]=0;
			rdct[i][j]=1.0*pixelmat[i][j];
		}
	}
	for(i=0;i<N;i++){
		reverse(rdct[i],idct[i],0);
		fft(rdct[i],idct[i]);
	}
	for(j=0;j<N;j++){
		for(i=0;i<N;i++){
			rcolumn[i]=rdct[i][j]*1.0/N;
			icolumn[i]=idct[i][j]*1.0/N;
		}
		reverse(rcolumn,icolumn,1);
		fft(rcolumn,icolumn);
		for(i=0;i<N;i++){
			rdct[i][j]=rcolumn[i];
			idct[i][j]=icolumn[i];
		}
	}
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			if(pow(i*i+j*j,0.5)>50){
				if(pow((i-N)*(i-N)+j*j,0.5)>50){
					if(pow((j-N)*(j-N)+i*i,0.5)>50){
						if(pow((i-N)*(i-N)+(j-N)*(j-N),0.5)>50){
							rdct[i][j]=0;
							idct[i][j]=0;
						}
					}
				}
			}
		}
	}

	idft(rdct,idct);

	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			pixelmat[i][j]=rdct[i][j];
		}
	}

	return pixelmat;
}

//sobel算子, 返回处理后的图像
int** sobel(int** pixelmat, int mheight, int mwidth)
{
	int sobel[9] = {1,2,1,
		0,0,0,
		-1,-2,-1};
	int i;
	int ker = 3, halfk = ker/2;	//模板核size = ker X ker
	int** img_tmp = new int*[mheight + halfk + 1];
	for(i = 0; i < mheight; i++)
		img_tmp[i] = new int[mwidth + halfk + 1];

	img_tmp = img_padding(pixelmat, mheight, mwidth, ker);   //镜像补充

	pixelmat = kernal_res(img_tmp, pixelmat, mheight, mwidth, sobel, ker, 1);
	return pixelmat;
}

//laplace算子, 返回处理后的图像
int** laplace(int** pixelmat, int mheight, int mwidth)
{
	int laplace[9] = {0,-1,0,
		-1,4,-1,
		0,-1,0};
	int i;
	int ker = 3, halfk = ker/2;	//模板核size = ker X ker
	int** img_tmp = new int*[mheight + halfk + 1];
	for(i = 0; i < mheight; i++)
		img_tmp[i] = new int[mwidth + halfk + 1];

	img_tmp = img_padding(pixelmat, mheight, mwidth, ker);   //镜像补充

	pixelmat = kernal_res(img_tmp, pixelmat, mheight, mwidth, laplace, ker, 1);
	return pixelmat;
}

//理想高通滤波, 返回处理后的图像
int** highpassfit(int** pixelmat, int mheight, int mwidth)
{
	return NULL;
}

//示例, 将图像平移到显示区域的中心
int** centralize(int** framemat, int** pixelmat, int mheight, int mwidth)
{
	//framemat为指向显示区域(画板)的二维数组指针, 大小为FRAME_HEIGHT x FRAMEWIDTH = 800 x 800
	int xpt = (FRAME_HEIGHT - mheight) / 2;
	int ypt = (FRAME_WIDTH - mwidth) / 2;
	for(int i = 0; i < mheight; i++)
	{
		for(int j = 0; j < mwidth; j++)
		{
			framemat[i + xpt][j + ypt] = pixelmat[i][j];
		}
	}
	return framemat;
}


//旋转图像, 返回显示区域(画板)指针
int** rotation(int** framemat, int** pixelmat, int mheight, int mwidth)
{
	int i;
	return NULL;
}

//平移图像, 返回显示区域(画板)指针
int** moveimage(int** framemat, int** pixelmat, int mheight, int mwidth)
{
	int xm = -100;
	int ym = 10;
	int i, j;
	for(i = 0; i < mheight; i++)
		for(j = 0; j < mwidth; j++)
		{
			if((i + xm) >= FRAME_HEIGHT || (j + ym) >= FRAME_WIDTH || (i + xm) < 0 || (j + ym) < 0)
				continue;
			else
				framemat[i + xm][j + ym] = pixelmat[i][j];
		}
		return framemat;
}

//缩放图像, 返回显示区域(画板)指针
int** scaling(int** framemat, int** pixelmat, int mheight, int mwidth)
{
	int i, j, is, js, flag = 0;
	float hs = 1.5;
	float ws = 1.5;
	int h = mheight*hs;
	int w = mwidth*ws;

	for(i = 0; i < h; i++)
		for(j = 0; j < w; j++)
			framemat[i][j] = -1;
	for(i = 0; i < mheight; i++)
		for(j = 0; j < mwidth; j++)
		{
			is=i*hs;
			js=j*ws;
			if(framemat[is][js] == -1)
				framemat[is][js] = pixelmat[i][j];
			else
				continue;
		}
		for(i = 0; i < h; i++)
		{
			if(framemat[i][0] != -1)
				flag = 1;
			else
				flag = 0;
			for(j = 0; j < w; j++)
				if(flag)
					framemat[i][j] = (framemat[i][j] == -1) ? framemat[i][j - 1] : framemat[i][j];
				else
					framemat[i][j] = (framemat[i][j] == -1) ? framemat[i - 1][j] : framemat[i][j];
		}
		return framemat;
}

//DFT变换, 返回处理后的图像, 注意缩放到0~255的整型
int** DFT(int** pixelmat, int mheight, int mwidth)
{
	return NULL;
}

//DCT变换, 返回处理后的图像
int** DCT(int** pixelmat, int mheight, int mwidth)
{
	int i, j, u, v;
	float sum, tmp = pi/(2*mheight);
	float au, av, a[2];// = {sqrt(1/mheight), sqrt(2/mwidth)};
	a[0] = 1;
	a[1] = sqrt(2.0);
	int** img_tmp = new int*[mheight];
	for(i = 0; i < mheight; i++)
		img_tmp[i] = new int[mwidth];

	for(u = 0; u < mheight; u++)
		for(v = 0; v < mwidth; v++)
		{
			sum = 0;
			au = (u == 0) ? a[0] : a[1];
			av = (v == 0) ? a[0] : a[1];
			for(i = 0; i < mheight; i++)
				for(j = 0; j < mwidth; j++)
					sum += pixelmat[i][j] * cos((2*i + 1)* u * tmp) * cos((2*j + 1) * v * tmp);
			img_tmp[u][v] = (int)(sum * au * av / mheight);
		}
	return img_tmp;
}

//walsh变换, 返回处理后的图像
int** walsh(int** pixelmat, int mheight, int mwidth)
{
	return NULL;
}

//haar变换, 返回处理后的图像
int** haar(int** pixelmat, int mheight, int mwidth)
{
	return NULL;
}

//生成随机噪声, 返回处理后的图像;
int** randomnoise(int** pixelmat, int mheight, int mwidth)
{
	float e = 4, var = 12;
	int i, j;
	for(i = 0; i < mheight; i++)
		for(j = 0; j < mwidth; j++)
		{
			pixelmat[i][j] += int(var * gaussrand() + e);
			pixelmat[i][j] = pixelmat[i][j] > 255 ? 255 : pixelmat[i][j];
			pixelmat[i][j] = pixelmat[i][j] < 0 ? 0 : pixelmat[i][j];
		}
		return pixelmat;
}


//生成椒盐噪声, 返回处理后的图像
int** impulsenoise(int** pixelmat, int mheight, int mwidth)
{
	int num_sp = mheight * mwidth / 50;
	int i, j ;
	for(int k = 0; k < num_sp; k++)
	{
		i = rand() % mheight;
		j = rand() % mwidth;
		pixelmat[i][j] = 255;
	}
	for(int k = 0; k < num_sp; k++)
	{
		i = rand() % mheight;
		j = rand() % mwidth;
		pixelmat[i][j] = 0;
	}
	return pixelmat;
}

//逆滤波复原
int** inversefit(int** pixelmat, int mheight, int mwidth)
{
	return NULL;
}

//维纳滤波
int** wienerfit(int** pixelmat, int mheight, int mwidth)
{
	return NULL;
}


//示例: JPEG压缩及解压缩
int** jpeg(int** pixelmat, int mheight, int mwidth)
{
	return NULL;
}
