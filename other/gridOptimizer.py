import numpy as np
import os
import random
from scipy.spatial.distance import cdist
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import lapjv
import gridlayoutOpt
import time


class gridOptimizer(object):
    def __init__(self):
        super().__init__()

    def solveKM(self, cost_matrix):
        # print('start KM')
        # row_asses, _, _ = lapjv.lapjv(cost_matrix)
        row_asses = np.array(gridlayoutOpt.solveKM(cost_matrix))
        N = row_asses.shape[0]
        col_asses = np.zeros(shape=N, dtype='int')
        for i in range(N):
            col_asses[row_asses[i]] = i
        # print('end KM')
        return row_asses, col_asses

    def solveJV(self, cost_matrix):
        row_asses = np.array(gridlayoutOpt.solveLap(cost_matrix, True, 50))
        N = row_asses.shape[0]
        col_asses = np.zeros(shape=N, dtype='int')
        for i in range(N):
            col_asses[row_asses[i]] = i
        # print('end KM')
        return row_asses, col_asses

    # DFS检查联通性
    def checkConnectDFS(self, row_asses, labels, x, y, lb, checked, num, square_len, connect):
        gid = x * square_len + y
        if checked[gid] == 1:
            return
        checked[gid] = 1
        f_x = [-1, 0, 0, 1, 1, -1, 1, -1]
        f_y = [0, -1, 1, 0, 1, 1, -1, -1]
        for i in range(connect):
            x2 = x + f_x[i]
            y2 = y + f_y[i]
            if x2 < 0 or x2 >= square_len or y2 < 0 or y2 >= square_len:
                continue
            gid2 = x2 * square_len + y2
            if row_asses[gid2] < num and labels[row_asses[gid2]] == lb:
                self.checkConnectDFS(row_asses, labels, x2, y2, lb, checked, num, square_len, connect)

    # 检查联通性
    def checkConnect(self, row_asses, labels, gid, show=False, connect=4):
        num = labels.shape[0]
        if row_asses[gid] >= num:
            return 0
        id = row_asses[gid]
        N = row_asses.shape[0]
        square_len = math.ceil(np.sqrt(N))
        checked = np.zeros(N, dtype='int')
        lb = labels[id]
        x = int(gid / square_len)
        y = gid - square_len * x
        self.checkConnectDFS(row_asses, labels, x, y, lb, checked, num, square_len, connect)
        per = checked.sum() / (labels == lb).sum()
        # if show:
        #     print('show per', checked.sum(), (labels == lb).sum())
        if per > 0.55:
            return 0
        return 1

    # 优化gridlayout
    def grid_op(self, ori_embedded, row_asses, labels, useGlobal=True, useLocal=True, convex_type="E", maxit=10, maxit2=5, only_compact=False, swap_cnt=2147483647):

        start = time.time()
        N = row_asses.shape[0]

        def solve_op(ori_embedded, row_asses, type, alpha, beta, alter=False, alter_best=None, maxit=5, maxit2=2, swap_cnt=2147483647):
            if alter_best is None:
                alter_best = [1, 1, 1, 1]

            # 迭代二分图匹配优化
            new_row_asses = row_asses.copy()
            new_row_asses2 = row_asses.copy()
            ans_row_asses = row_asses.copy()
            new_cost = np.array([2147483647, 2147483647, 2147483647])
            best_cost = new_cost.copy()
            print("````````````````````")
            print('alpha', alpha, beta)
            print("````````````````````")
            change = np.ones(shape=N, dtype='bool')

            if maxit > 0:
                # tmp_row_asses = np.array(gridlayoutOpt.optimizeBA(ori_row_asses, row_asses, labels, change, type, alpha, beta, maxit))
                tmp_row_asses = np.array(gridlayoutOpt.optimizeBA(ori_embedded, row_asses, labels, change, type, alpha, beta, alter, alter_best, maxit))
                for i in range(N):
                    new_row_asses[i] = round(tmp_row_asses[i])
                new_cost = np.array([tmp_row_asses[N], tmp_row_asses[N+1], tmp_row_asses[N+2]])
                ans_row_asses = new_row_asses.copy()
                best_cost = new_cost.copy()
                # print("cost1", best_cost)

            # 枚举交换优化
            if maxit2 >= 0:
                # alpha = 1
                # beta = 0
                change = np.ones(shape=N, dtype='bool')
                seed_list = [10]
                # if (type == "E") and alter:
                # if type == "E":
                #     seed_list = [0, 10, 20]
                flag = 0
                for seed in seed_list:
                    tmp_row_asses = np.array(
                        gridlayoutOpt.optimizeSwap(ori_embedded, new_row_asses, labels, change, type, alpha, beta,
                                                maxit2, seed, True, swap_cnt))
                    for i in range(N):
                        new_row_asses2[i] = round(tmp_row_asses[i])
                    if tmp_row_asses[N] != -1:
                        new_cost = np.array([tmp_row_asses[N], tmp_row_asses[N+1], tmp_row_asses[N+2]])
                        alpha_full = np.array([1-alpha-beta, beta, alpha])
                        print("seed", seed, new_cost, (alpha_full*new_cost).sum(), (alpha_full*best_cost).sum())
                        if flag == 0 or (alpha_full*new_cost).sum() <= (alpha_full*best_cost).sum():
                            best_cost = new_cost.copy()
                            ans_row_asses = new_row_asses2.copy()
                            flag = 1

            return ans_row_asses, best_cost

        compact_it = 3
        global_it = 5
        if convex_type == "O":
            compact_it = 0
            global_it = 0

        ans = row_asses
        new_cost = np.array([2147483647, 2147483647, 2147483647])
        if useGlobal:
            row_asses1, new_cost1 = solve_op(ori_embedded, row_asses, "Global", 0, 0, False, None, 0, 0)
            # self.show_grid(row_asses1, labels, square_len, "1.png", False, False)
            print("new_cost1", new_cost1)
            print("time1", time.time()-start)
            row_asses2, new_cost2 = solve_op(ori_embedded, row_asses, "Global", 0, 1, False, None, global_it, 0)
            # self.show_grid(row_asses2, labels, square_len, "2.png", False, False)
            print("new_cost2", new_cost2)
            print("time2", time.time()-start)

            alter_best = [0, 0, 1, 1]
            alter_best[0] = min(new_cost1[0], new_cost2[0])
            alter_best[2] = max(new_cost1[0], new_cost2[0])-min(new_cost1[0], new_cost2[0])
            alter_best[1] = min(new_cost1[1], new_cost2[1])
            alter_best[3] = max(new_cost1[1], new_cost2[1])-min(new_cost1[1], new_cost2[1])

            print(alter_best)

            if only_compact:
                ans, new_cost = row_asses2, new_cost2
            elif alter_best[3] == 0:
                ans, new_cost = row_asses1, new_cost1
            else:
                ans, new_cost = solve_op(ori_embedded, row_asses, "Global", 0.33, 0.33, True, alter_best, global_it, 0)
                # self.show_grid(ans, labels, square_len, "4.png", False, False)
            print("new_cost4", new_cost)
            print("time4", time.time()-start)

        end2 = time.time()
        t2 = end2 - start

        if useLocal:

            # data = np.load('T-global.npz')
            # ori_row_asses = data['row_asses']
            # ans = data['row_asses']

            # if convex_type == "E":
            #     ans, new_cost = solve_op(ori_embedded, ans, "C", 1, 0, False, None, maxit, maxit2, 2147483647)
            # else:
            #     swap_cnt = 2147483647

            # if convex_type == "C":
            #     ans, new_cost = solve_op(ori_embedded, ans, "E", 1, 0, False, None, maxit, maxit2, 2147483647)
            #     maxit2 = 1
            #     swap_cnt = 1

            if convex_type == "CE":
                ans, new_cost = solve_op(ori_embedded, ans, "C", 1, 0, False, None, 0, 3, 2147483647)
                convex_type = "E"

            if convex_type == "EC":
                ans, new_cost = solve_op(ori_embedded, ans, "E", 1, 0, False, None, 0, 3, 2147483647)
                convex_type = "C"

            if convex_type == "ST":
                ans, new_cost = solve_op(ori_embedded, ans, "S", 1, 0, False, None, 0, 3, 2147483647)
                convex_type = "T"

            if convex_type == "TS":
                ans, new_cost = solve_op(ori_embedded, ans, "T", 1, 0, False, None, 5, 0, 2147483647)
                convex_type = "S"

            ans, new_cost = solve_op(ori_embedded, ans, convex_type, 1, 0, False, None, maxit, maxit2, swap_cnt)

            print("new_cost5", new_cost)
            print("time5", time.time()-start)

        end = time.time()
        t1 = end - start

        # 计算其他指标
        # if convex_type != "T":
        #     _, new_cost2 = solve_op(ori_embedded, ans, "T", 1, 0, False, None, 0, 0)
        #     new_cost = np.append(new_cost, [new_cost2[2]], None)
        #     new_cost[2], new_cost[3] = new_cost[3], new_cost[2]
        # else:
        #     # _, new_cost2 = solve_op(ori_embedded, ans, "E", 1, 0, False, None, 0, 0)
        #     _, new_cost2 = solve_op(ori_embedded, ans, "2020", 1, 0, False, None, 0, 0)
        #     new_cost = np.append(new_cost, [new_cost2[2]], None)

        _, new_cost2 = solve_op(ori_embedded, ans, "T", 1, 0, False, None, 0, 0)
        new_cost[2] = new_cost2[2]
        _, new_cost2 = solve_op(ori_embedded, ans, "S", 1, 0, False, None, 0, 0)
        new_cost = np.append(new_cost, [new_cost2[2]], None)
        _, new_cost2 = solve_op(ori_embedded, ans, "C", 1, 0, False, None, 0, 0)
        new_cost = np.append(new_cost, [new_cost2[2]], None)
        _, new_cost2 = solve_op(ori_embedded, ans, "E", 1, 0, False, None, 0, 0)
        new_cost = np.append(new_cost, [new_cost2[2]], None)

        return ans, t1, t2, new_cost

    # 生成gridlayout
    def grid(self, X_embedded: np.ndarray, labels: np.ndarray = None, type='E', maxit=10, maxit2=5, use_global=True, only_compact=False, swap_cnt=2147483647):
        # 初始化信息
        start = time.time()
        ans = None
        best = 2147483647
        X_embedded -= X_embedded.min(axis=0)
        X_embedded /= X_embedded.max(axis=0)
        num = X_embedded.shape[0]
        square_len = math.ceil(np.sqrt(num))
        print('num', num)
        N = square_len * square_len
        maxLabel = 0
        for id in range(num):
            label = labels[id]
            maxLabel = max(maxLabel, label + 1)

        labelCheckList = np.zeros(shape=(maxLabel + 1, maxLabel + 1), dtype=bool)
        for i in range(maxLabel):
            for j in range(maxLabel):
                if i != j:
                    labelCheckList[i][j] = True

        def getDist(x1, y1, x2, y2):
            return np.sqrt(np.square(x1 - x2) + np.square(y1 - y2))

        # 检查近邻关系变化情况
        def checkNeighbor(old_row_asses, col_asses, show=False):
            dist_inc_list = []
            dist_inc_list_same = []
            for old_gid1 in range(N):
                for old_gid2 in range(N):
                    if old_gid1 >= old_gid2:
                        continue
                    old_x1 = int(old_gid1 / square_len)
                    old_y1 = old_gid1 - square_len * old_x1
                    old_x2 = int(old_gid2 / square_len)
                    old_y2 = old_gid2 - square_len * old_x2
                    old_dist = getDist(old_x1, old_y1, old_x2, old_y2)
                    if old_dist > 2:
                        continue
                    id1 = old_row_asses[old_gid1]
                    id2 = old_row_asses[old_gid2]
                    if id1 >= num or id2 >= num:
                        continue
                    gid1 = col_asses[id1]
                    x1 = int(gid1 / square_len)
                    y1 = gid1 - square_len * x1
                    gid2 = col_asses[id2]
                    x2 = int(gid2 / square_len)
                    y2 = gid2 - square_len * x2
                    dist_inc_list.append(getDist(x1, y1, x2, y2) - old_dist)
                    if labels[id1] == labels[id2]:
                        dist_inc_list_same.append(getDist(x1, y1, x2, y2) - old_dist)
            dist_inc_list = np.array(dist_inc_list)
            dist_inc_list_same = np.array(dist_inc_list_same)
            tot = dist_inc_list.shape[0]
            tot_same = dist_inc_list_same.shape[0]
            if show:
                X = []
                Y = []
                X_same = []
                Y_same = []
                for x in range(1000):
                    x /= 100
                    X.append(x)
                    Y.append((dist_inc_list <= x).sum() / tot)
                    X_same.append(x)
                    Y_same.append((dist_inc_list_same <= x).sum() / tot_same)
                plt.plot(X, Y)
                plt.plot(X_same, Y_same)
                plt.show()
                print('all', dist_inc_list.sum() / tot)
                print('same', dist_inc_list_same.sum() / tot_same)
            return dist_inc_list, dist_inc_list_same


        # 生成初始layout
        grids = np.dstack(np.meshgrid(np.linspace(0, 1 - 1.0 / square_len, square_len),
                                      np.linspace(0, 1 - 1.0 / square_len, square_len))) \
            .reshape(-1, 2)

        tmp = grids[:,0].copy()
        grids[:,0] = grids[:,1]
        grids[:,1] = tmp

        # print(grids)

        original_cost_matrix = cdist(grids, X_embedded, "euclidean")
        # knn process
        dummy_points = np.ones((N - original_cost_matrix.shape[1], 2)) * 0.5
        # dummy at [0.5, 0.5]
        dummy_vertices = (1 - cdist(grids, dummy_points, "euclidean")) * 100
        cost_matrix = np.concatenate((original_cost_matrix, dummy_vertices), axis=1)

        cost_matrix = np.power(cost_matrix, 2)

        # row_asses, col_asses, info = fastlapjv(cost_matrix, k_value=50 if len(cost_matrix)>50 else len(cost_matrix))
        # row_asses, col_asses = self.solveKM(cost_matrix)
        row_asses, col_asses = self.solveJV(cost_matrix)
        # col_asses = col_asses[:num]
        # self.show_grid(row_asses, labels, square_len, 'new.png')

        # 补全labels
        # labels = np.concatenate((labels, np.full(shape=(N-num), dtype='int', fill_value=maxLabel)), axis=0)
        # num = N
        # maxLabel = maxLabel+1

        old_X_embedded = X_embedded
        ori_row_asses = row_asses.copy()
        ori_embedded = grids[col_asses]
        ori_labels = labels.copy()

        t0 = time.time()-start
        # 简单聚类
        print("start cluster")
        labels = np.array(gridlayoutOpt.getClusters(ori_row_asses, labels))
        maxLabel = labels.max()+1
        print("end cluster")
        # self.show_grid(row_asses, ori_labels, square_len, 'new0.png')
        # self.show_grid(row_asses, labels, square_len, 'new1.png')
        # return row_asses, 0, 0

        # datas = np.load("T-base.npz")
        # labels = datas['labels']
        # row_asses = datas['row_asses']
        # ori_row_asses = row_asses

        # data = np.load('T-global.npz')
        # labels = data['labels']

        # 开始优化
        print("start optimize")
        print("--------------------------------------------------")
        start = time.time()

        ans, t1, t2, new_cost = self.grid_op(ori_embedded, row_asses, labels, use_global, True, type, maxit, maxit2, only_compact, swap_cnt)

        end = time.time()
        print("end optimize")
        print("--------------------------------------------------")
        print('time:', end - start)

        # self.show_grid(ans, labels, square_len, 'new1.png')
        # change_list = []
        # for i in range(N):
        #     if (ans[i]<num)and(labels[ans[i]]==0):
        #         change_list.append(i)
        # change_list = np.array(change_list)
        # ans = self.translateAdjust(ori_embedded, ans, labels, use_global, True, type, int(maxit*1.5), int(maxit2*1.5), change_list, np.array((0, 0.3)))
        #
        # self.show_grid(ans, labels, square_len, 'new3.png')

        # self.show_grid(row_asses, ori_labels, square_len, 'new3.png')
        # np.savez("tmp.npz", row_asses=row_asses, labels=labels)

        if (maxit+maxit2) == 0:
            t1 = t2
        return ans, t1, t0, labels, new_cost

    def translateAdjust(self, ori_embedded, row_asses, labels, useGlobal=True, useLocal=True, convex_type="E", maxit=5, maxit2=2, change_list=np.array([]), translate=np.array([0, 0])):
        N = row_asses.shape[0]
        square_len = round(np.sqrt(N))
        num = labels.shape[0]

        # X_embedded = np.zeros((N, 2))
        # for x in range(square_len):
        #     bias = x*square_len
        #     for y in range(square_len):
        #         gid = bias+y
        #         X_embedded[ori_row_asses[gid]][0] = x*1.0/square_len
        #         X_embedded[ori_row_asses[gid]][1] = y*1.0/square_len
        X_embedded = ori_embedded.copy()

        X_embedded[row_asses[change_list]] += translate

        # 生成初始layout
        grids = np.dstack(np.meshgrid(np.linspace(0, 1 - 1.0 / square_len, square_len),
                                      np.linspace(0, 1 - 1.0 / square_len, square_len))) \
            .reshape(-1, 2)

        tmp = grids[:, 0].copy()
        grids[:, 0] = grids[:, 1]
        grids[:, 1] = tmp

        cost_matrix = cdist(grids, X_embedded, "euclidean")

        cost_matrix = np.power(cost_matrix, 2)

        # row_asses, col_asses, info = fastlapjv(cost_matrix, k_value=50 if len(cost_matrix)>50 else len(cost_matrix))
        new_ori_row_asses, col_asses = self.solveJV(cost_matrix)

        # self.show_grid(new_ori_row_asses, labels, square_len, 'new2.png')

        new_row_asses, _, _, _ = self.grid_op(X_embedded, new_ori_row_asses, labels, useGlobal, useLocal, convex_type, maxit, maxit2)
        # print("done")
        change = np.ones(shape=N, dtype='bool')
        new_row_asses = np.array(gridlayoutOpt.optimizeInnerCluster(ori_embedded, new_row_asses, labels, change))
        return new_row_asses

    # 检查凸性
    def checkConvex(self, row_asses, labels, type='T', labelCheckList=None, edgePer=1):
        if type == 'T':
            T_list = np.array(gridlayoutOpt.checkConvexForT(row_asses, labels))
            score1 = 0
            if T_list[1] > 0:
                score1 = T_list[0] / T_list[1]
            return score1
        elif type == 'E':
            E_list = np.array(gridlayoutOpt.checkConvexForE(row_asses, labels))
            score1 = 0
            if E_list[1] > 0:
                score1 = E_list[0] / E_list[1]
            return score1

    # 绘图
    def show_grid(self, row_asses, grid_labels, square_len, path='new.png', showNum=True, just_save=False):
        def highlight_cell(x, y, ax=None, **kwargs):
            rect = plt.Rectangle((x - .5, y - .5), 1, 1, fill=False, **kwargs)
            ax = ax or plt.gca()
            ax.add_patch(rect)
            return rect

        from colors import MyColorMap
        cm = MyColorMap()

        data = []
        num = math.ceil(square_len)
        for i in range(num - 1, -1, -1):
            row = []
            for j in range(num):
                row.append(plt.cm.tab20(grid_labels[row_asses[num * i + j]]))
                # row.append(cm.color(grid_labels[row_asses[num * i + j]]))
            data.append(row)
        plt.cla()
        plt.imshow(data)
        for i in range(num - 1, -1, -1):
            for j in range(num):
                highlight_cell(i, j, color="white", linewidth=1)
        if showNum:
            for i in range(num):
                for j in range(num):
                    text = plt.text(j, num - i - 1, row_asses[num * i + j], fontsize=7, ha="center", va="center",
                                    color="w")
        plt.axis('off')
        plt.savefig(path)
        if not just_save:
            plt.show()