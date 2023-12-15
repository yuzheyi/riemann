#include <cmath>
#include <iostream> 
int fibonacci(int n) {
    // 基本情况：0的阶乘为1
    if (n == 0) {
        // std::cout << "n=0" << std::endl;
        return 2;
    }
    // else if (n==1){
    //     // std::cout << "n=1" << std::endl;
    //     return 1;
    // }

    // 递归情况：n的阶乘等于n乘以(n-1)的阶乘
    else {
        // std::cout << result << " "; // 输出每级的数
        return  fibonacci(n - 1) - n +3;
    }
}
void fibonacciReverse(int n , int m) {
    if (n <= 0) {
        std::cout << "Invalid input. Please enter a positive integer." << std::endl;
    } 
    else if (m == 0){//1题目的奇数项和
        int sum = 0;
        for (int i = n; i >= 0; --i) {
            int result = fibonacci(i);
            // sum += result;
            std::cout <<"数列"<<i<<"结果："<< result << " ";
        }
        std::cout << std::endl;
        // std::cout << "Sum: " << sum << std::endl;
        // y1 = y1*q1;
        // y2 = y2*q2;
        // q1 = q1*q1;
        // q2 = q2*q2;
        // double mysum = y1*(1-pow(q1,(n+1)))/(1-q1)+y2*(1-pow(q2,(n+1)))/(1-q2);
        // std::cout << "特征数列Sum: " << mysum << std::endl;
    }
    else {
        std::cout << "Invalid input. Please enter 0 or 1." << std::endl;
    }
}

int main() {
    int num;
    int mode=0 ;
    std::cout << "数列: ";
    std::cin >> num;

    if (num < 0) {
        std::cout << "Please enter a non-negative integer." << std::endl;
    } else {
        fibonacciReverse(num,mode);
    }
}