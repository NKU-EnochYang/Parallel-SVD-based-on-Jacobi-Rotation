# JacobiSVD程序说明

```bash
.
├── BlockParallel # 分块并行算法
│   ├── BlockParallel.cpp # 分块并行算法，实现了矩阵的旋转及奇异值的求解
│   └── BlockParallel.h # 分块并行算法头文件，包括rotate()和SVD()两个函数
├── CMakeLists.txt # CMake文件，内含编译信息
├── common.h # 公有的一些头文件及全局变量，在此更改矩阵规模、线程数目等参数
├── main.cpp # 主函数，内含对矩阵的初始化及三个算法的测试
├── Matrix #矩阵类
│   ├── Matrix.cpp # 一些常用的矩阵操作
│   └── Matrix.h # 矩阵类头文件
├── Parallel # 直接并行算法
│   ├── Parallel.cpp # 直接并行算法，实现了矩阵的旋转及奇异值的求解
│   └── Parallel.h # 直接并行算法头文件，包括rotate()和SVD()两个函数
└── Serial # 串行算法
    ├── Serial.cpp # 串行算法，实现了矩阵的旋转及奇异值的求解
    └── Serial.h # 串行算法头文件，包括rotate()和SVD()两个函数
```