#include "gauss.h"
#include <iostream>
#include <cmath>
#include <vector>

void gaussSolver(std::vector<std::vector<double>>& a, std::vector<double>& b,std::vector<Cell>& cells) //a 为方阵 ，b 为列向量
//求线性方程组的解(ax=b ,求 x)，矩阵 a 为方阵并且方程组有唯一解时返回 true
{
	/*构造增广矩阵*/
	int raws=a.size();
	int cols=a[0].size();
	std::vector<std::vector<double>> augmat(raws, std::vector<double>(cols+1, 0.0));
	for (int i = 0; i < raws; i++) {
		for (int j = 0; j < cols; j++) {
			augmat[i][j] = a[i][j];
		}
	}
	for (int i = 0; i < raws; i++) {
		augmat[i][cols] = b[i];
	}



	//下面代码将增广矩阵化为上三角矩阵，并判断增广矩阵秩是否为 n


	if (cols != raws)//只求解是方阵时的情形
	{
		std::cout << "系数矩阵不是方阵" << std::endl;
		return; 
	}

	//下面代码将增广矩阵化为上三角矩阵，并判断增广矩阵秩是否为 n
	for (int i = 0; i < cols; i++)
	{
		std::cout << "上三角矩阵化简第："<<i<<"列"<<std::endl; 
		//寻找第 i 列不为零的元素
		int k;
		for (k = i; k < cols; k++)
		{
			if (fabs(augmat[k][i]) > 1e-10) //满足这个条件时，认为这个元素不为0
				break;
		}
		if (k < cols)//说明第 i 列有不为0的元素
		{
			//交换第 i 行和第 k 行所有元素
		for (int j = i; j < cols-1 ; j++) {
    		std::swap(augmat[i][j], augmat[k][j]);
		}
			for (int j = i + 1; j < cols; j++)
			{
				double c = -augmat[j][i] / augmat[i][i];
				std::vector<double>& row_i = augmat[i]; // 指向第 i 行
    			std::vector<double>& row_j = augmat[j]; // 指向第 j 行
    			for (k = i; k < cols + 1; k++) 
				{
        			row_j[k] += c * row_i[k];
    			}
			}
		}
		else //没有找到则说明系数矩阵秩不为 n ，说明方程组中有效方程的个数小于 n
		{
			std::cout << "系数矩阵奇异，线性方程组无解或有无数解" << std::endl;
			return;
		}
	}
	//自下而上求解
	for (int i = raws-1; i >= 0; i--)
	{
		std::cout << "求解第："<<i<<"列"<<std::endl; 
		double result = augmat[i][cols];
		for (int j = raws-1; j > i; j--)
		{

			result =result - cells[j].physics * augmat[i][j];
		}
		if (fabs(augmat[i][i]) < 1e-10) 
		{
            throw std::runtime_error("奇异矩阵，线性方程组无解或有无数解");
        }
		result /= augmat[i][i];
		cells[i].physics=result;
		std::cout << cells[i].physics << ",";
	}
	return;
}
 