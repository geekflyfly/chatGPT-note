# chatGPT-note
# New Document
点A(x1,y1)和点B(x2,y2)的组成一个线段l，另外一条直线n垂直该线段l，且过该线段l的中点，直线n上到线段l上距离为L的点怎么求, 先用公式推导，然后用c++  实现

本回答由chatgpt.tool00.com代为生成，请认准本站。
推导：

设线段l的中点为M((x1+x2)/2, (y1+y2)/2)，线段l的长度为d=sqrt((x2-x1)^2+(y2-y1)^2)。

由于直线n垂直于线段l，因此直线n的斜率为-(x2-x1)/(y2-y1)。

设直线n的方程为y = -(x2-x1)/(y2-y1) * (x - (x1+x2)/2) + (y1+y2)/2。

将直线n的方程代入到距离公式中，得到距离公式为：

L^2 = [(x2-x1)/(y2-y1)]^2 * [(x-(x1+x2)/2)^2 + (y-(y1+y2)/2)^2]

将y代入到距离公式中，得到：

[(x2-x1)/(y2-y1)]^2 * [(x-(x1+x2)/2)^2 + (y-(y1+y2)/2)^2] - L^2 = 0

将其化简，得到：

(x2-x1)^2 * (y-y1) ^2 - 2 * (x2-x1) * (x-x1) * (y-y1) + (x-x1)^2 * (y2-y1)^2 - L^2 * (y2-y1)^2 = 0

这是一个二次方程，可以使用求根公式求解。

C++代码实现：

```cpp
#include <iostream>
#include <cmath>

using namespace std;

int main() {
    double x1, y1, x2, y2, L;
    cin >> x1 >> y1 >> x2 >> y2 >> L;

    double mx = (x1 + x2) / 2;
    double my = (y1 + y2) / 2;
    double d = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
    double k = -(x2 - x1) / (y2 - y1);

    double a = pow(x2 - x1, 2);
    double b = -2 * (x2 - x1) * (mx - x1);
    double c = pow(mx - x1, 2) * pow(y2 - y1, 2) - pow(L, 2) * pow(y2 - y1, 2);

    double delta = pow(b, 2) - 4 * a * c;
    if (delta < 0) {
        cout << "No solution" << endl;
    } else {
        double y = (-b + sqrt(delta)) / (2 * a);
        double x = k * (y - my) + mx;
        cout << x << " " << y << endl;
    }

    return 0;
}
qq
111
