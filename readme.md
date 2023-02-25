##### 代码接口

每个实现中保留一个类似如下对外接口：

```C
std::vector<double> checkConvexForT(
const std::vector<int> &_grid_asses, // 网格位置数据
const std::vector<int> &_cluster_labels); // 标准数据
```

其余简洁、准确即可



##### 文件管理

- convexity 文件夹放置measure核心代码
- utils 文件夹放置通用函数代码
- test 文件夹放置使用的测试代码
- other 文件夹及setup,py、mainop_alter2.cpp等文件
  - 目前是 pybind11 及 python 中的一些实际使用，代码编写中可以不使用，也可以用于测试
  - 如果用该方式测试，请一定**不要把生成的库文件及build文件夹上传git**



##### 代码管理

每次新项目前**开一个分支**，可以用自己名字，也可以用measure名字

可以参考 ./git.png 使用git管理代码，但可以不用那么严格



##### 记录完成情况

做之前/做完一个measure在下面的文档更新一下状态，保证不重复工作即可：

【腾讯文档】Convexity Measure 完成情况
https://docs.qq.com/sheet/DRlVPYkhNd1ZQVkNF?tab=BB08J2