# ConvexityMeasure 

### 基本介绍
本仓库使用C++实现了一些不同论文中涉及的2D凸性度量，提供了网格布局中相对统一的凸性接口调用:

代码是以下工作的一部分 [https://github.com/thu-vis/Cluster-Aware-Grid-Layout](https://github.com/thu-vis/Cluster-Aware-Grid-Layout) ，即将发表在TVCG2023

![图片](https://github.com/ChenJash/ConvexityMeasure/assets/70832896/51ba3085-14b6-4ec2-a9d7-9c36b691f797)
![图片](https://github.com/ChenJash/ConvexityMeasure/assets/70832896/7bbe5f67-f390-474b-b5c0-329c47a878b8)



### 仓库细节
#### 代码接口

每个实现中保留一个类似如下对外接口：

```C
std::vector<double> checkConvexForT(
const std::vector<int> &_grid_asses, // 网格位置数据
const std::vector<int> &_cluster_labels); // 标准数据
```

#### 文件管理

- convexity 文件夹放置measure核心代码
- utils 文件夹放置通用函数代码
- test 文件夹放置使用的测试代码
- other 文件夹及setup,py、mainop_alter2.cpp等文件
  - 目前是 pybind11 及 python 中的一些实际使用，代码编写中可以不使用，也可以用于测试


#### 记录完成情况
【腾讯文档】Convexity Measure 完成情况
https://docs.qq.com/sheet/DRlVPYkhNd1ZQVkNF?tab=BB08J2
![图片](https://github.com/ChenJash/ConvexityMeasure/assets/70832896/8686d30f-b35b-40b4-94b2-8d1c1dddde85)
