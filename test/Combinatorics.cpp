#include <iostream>
#include <cmath>
// 递归计算阶乘
int fibonacci(int n) {
    // 基本情况：0的阶乘为1
    if (n == 0) {
        // std::cout << "n=0" << std::endl;
        return 0;
    }
    else if (n==1){
        // std::cout << "n=1" << std::endl;
        return 1;
    }

    // 递归情况：n的阶乘等于n乘以(n-1)的阶乘
    else {
        int result = fibonacci(n - 1) + fibonacci(n - 2);
        // std::cout << result << " "; // 输出每级的数
        return result;
    }
}

// 计算从 n 到 1 的斐波那契数列的函数
void fibonacciReverse(int n , int m) {
    double y1 = 1/sqrt(5);
    double y2 = -1/sqrt(5);
    double q1 = (1+sqrt(5))/2;
    double q2 = (1-sqrt(5))/2;
    if (n <= 0) {
        std::cout << "Invalid input. Please enter a positive integer." << std::endl;
    } 
    else if (m == 0){//1题目的奇数项和
        int sum = 0;
        for (int i = n; i >= 0; --i) {
            int result = fibonacci(2*i+1);
            sum += result;
            std::cout << result << " ";
        }
        std::cout << std::endl;
        std::cout << "Sum: " << sum << std::endl;
        y1 = y1*q1;
        y2 = y2*q2;
        q1 = q1*q1;
        q2 = q2*q2;
        double mysum = y1*(1-pow(q1,(n+1)))/(1-q1)+y2*(1-pow(q2,(n+1)))/(1-q2);
        std::cout << "特征数列Sum: " << mysum << std::endl;
    }
    else if (m == 1){
        int sum = 0;
        for (int i = n; i >= 0; --i) {
            int result = fibonacci(2*i);
            sum += result;
            std::cout << result << " ";
        }
        std::cout << std::endl;
        std::cout << "Sum: " << sum << std::endl;
        q1 = q1*q1;
        q2 = q2*q2;
        double mysum = y1*(1-pow(q1,(n+1)))/(1-q1)+y2*(1-pow(q2,(n+1)))/(1-q2);
        std::cout << "特征数列Sum: " << mysum << std::endl;
    }
    else if (m == 2){
    int sum = 0;
    for (int i = n; i >= 0; --i) {
        int result = pow(-1,i)*fibonacci(i);
        sum += result;
        std::cout << result << " ";
    }
    std::cout << std::endl;
    std::cout << "Sum: " << sum << std::endl;
    q1 = -q1;
    q2 = -q2;
    double mysum = y1*(1-pow(q1,(n+1)))/(1-q1)+y2*(1-pow(q2,(n+1)))/(1-q2);
    std::cout << "特征数列Sum: " << mysum << std::endl;
    }
    else if (m == 3){
    int sum = 0;
    for (int i = n; i >= 0; --i) {
        int result = pow(fibonacci(i),2);
        sum += result;
        std::cout << result << " ";
    }
    std::cout << std::endl;
    std::cout << "Sum: " << sum << std::endl;
    double y3 = 2*y1*y2;
    double q3 = q1*q2;
    y1 = pow(y1,2);
    y2 = pow(y2,2);
    q1 = pow(q1,2);
    q2 = pow(q2,2);
    double mysum = y1*(1-pow(q1,(n+1)))/(1-q1)+y2*(1-pow(q2,(n+1)))/(1-q2)+y3*(1-pow(q3,(n+1)))/(1-q3);
    std::cout << "特征数列Sum: " << mysum << std::endl;
    }
    else {
        std::cout << "Invalid input. Please enter 0 or 1." << std::endl;
    }
}

int main() {
    int num;
    int mode=3 ;
    std::cout << "斐波那契数列: ";
    std::cin >> num;

    if (num < 0) {
        std::cout << "Please enter a non-negative integer." << std::endl;
    } else {
        fibonacciReverse(num,mode);
    }
    

    double y= 0.5/sqrt(5)*(sqrt(5)+3)*pow(0.5*(1+sqrt(5)),num)+0.5/sqrt(5)*(sqrt(5)-3)*pow(0.5*(1-sqrt(5)),num);
    std::cout << "y: " << y << std::endl; 
    return 0;
}
