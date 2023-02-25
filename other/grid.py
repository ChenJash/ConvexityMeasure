dict = {}
dict2 = {}
dict0 = {}
save_file = True

def body():
    import numpy as np
    np.set_printoptions(suppress=True)
    np.set_printoptions(precision=4)
    import os
    processed_data = "./processed_data"
    tsne_result = "./tsne_result"

    # feature_path = os.path.join(processed_data, "{}_full.npz".format(dataset))
    tsne_path = os.path.join(tsne_result, dataset, "position.npz")
    # features = np.load(feature_path, allow_pickle=True)['features']
    tsne_file = np.load(tsne_path, allow_pickle=True)
    embeddings = tsne_file['positions']
    labels = tsne_file['labels']
    label_names = tsne_file['label_names']

    if dataset == "Cats-vs-dogs":
        label_names = label_names[:2]

    print(label_names)
    print(labels.shape[0])
    if labels.shape[0] < grid_width*grid_width:
        return

    size = labels.shape[0]
    # grid_width = 45
    # flag = 7
    path = "./grid_result"
    if not os.path.exists(path):
        os.mkdir(path)

    path = os.path.join(path, dataset)
    if not os.path.exists(path):
        os.mkdir(path)

    path = os.path.join(path, str(grid_width))
    if not os.path.exists(path):
        os.mkdir(path)

    path = os.path.join(path, str(flag) + "-" + str(grid_width) + "-" + str(tt))
    if save_file and not os.path.exists(path):
        os.mkdir(path)

    samples_N = grid_width * grid_width

    if flag == 'all':
        if (grid_width-1)*(grid_width-1) >= size:
            return
        if samples_N > size:
            samples_N = size
        samples = np.random.choice(size, samples_N, replace=False)
        np.save("samples_{}0.npy".format(dataset), samples)
    else:
        if flag >= label_names.shape[0]:
            return
        if label_names.shape[0] <= 4:
            return

        import collections
        print(collections.Counter(labels))

        chosen_labels = np.random.choice(label_names.shape[0], flag, replace=False)
        # chosen_labels = np.array([3, 0, 9])
        print(chosen_labels)
        idx = np.zeros(size, dtype='bool')
        for lb in chosen_labels:
            idx = (idx | (labels == lb))
        print(idx)
        idx = (np.arange(size))[idx]
        chosen_size = idx.shape[0]
        print('chosen_size', chosen_size)

        if (grid_width-1)*(grid_width-1) >= chosen_size:
            return

        if samples_N > chosen_size:
            samples_N = chosen_size

        samples = np.random.choice(chosen_size, samples_N, replace=False)
        samples = idx[samples]
        np.save("samples_{}0.npy".format(dataset), samples)

    samples = np.load("samples_{}0.npy".format(dataset))
    # s_features = features[samples]
    s_embeddings = embeddings[samples]
    print(s_embeddings.shape)
    s_labels = labels[samples]
    print(s_labels.shape)
    # print(samples.shape[0], s_features.shape)

    import numpy as np

    new_labels = s_labels.copy()

    from gridOptimizer import gridOptimizer
    import os

    # type = "T"
    # # showT = ""
    # showT = "-NoneText"

    # for type in ["O", "S", "T", "E", "C"]:
    # for type in ["T", "S", "ST", "TS"]:
    for type in ["T", "S", "ST", "TS", "C", "E", "EC", "CE"]:
    # for type in ["E"]:
        # for showT in ["", "-NoneText"]:
        for showT in ["-NoneText"]:

            # for op_type in ["base", "compact", "global", "full"]:
            for op_type in ["global", "full"]:
                if (type != "T") and (op_type != "full"):
                    continue
                m1 = 0
                m2 = 3
                # if alpha+beta == 0:
                #     m1 = 1
                #     m2 = 0
                # if type == "S" or type == "C" or type == "TB":
                #     m1 = 0
                #     m2 = 5
                if type == "E":
                    m1 = 0
                    m2 = 5
                if type == "C":
                    m1 = 0
                    m2 = 5
                if type == "T":
                    m1 = 10
                    m2 = 0
                if type == "ST":
                    m1 = 5
                    m2 = 0
                if type == "S":
                    m1 = 0
                    m2 = 5
                if type == "TS":
                    m1 = 0
                    m2 = 3
                if type == "O":
                    m1 = 0
                    m2 = 0

                if (op_type == "base") or (op_type == "global") or (op_type == "compact"):
                    m1 = 0
                    m2 = 0

                use_global = True
                if (op_type == "base") or (op_type == "local"):
                    use_global = False

                only_compact = False
                if op_type == "compact":
                    only_compact = True

                file_path = os.path.join(path, type + "-" + op_type + showT + ".png")
                save_path = os.path.join(path, type + "-" + op_type + ".npz")

                Optimizer = gridOptimizer()
                # print("check done", BASolver.checkConvex(np.array(row_asses_c), np.array(s_labels)))
                # row_asses_m, heat = BASolver.grid3(s_embeddings, s_labels, 'E')
                row_asses_m, t1, t2, new_labels, new_cost = Optimizer.grid(s_embeddings, s_labels, type, m1, m2,
                                                                           use_global, only_compact, swap_cnt=2147483647)

                show_labels = new_labels
                show_labels = np.array(show_labels)
                tmp = np.full(row_asses_m.shape[0] - show_labels.shape[0], dtype='int', fill_value=-1)
                show_labels = np.concatenate((show_labels, tmp), axis=0)
                print(row_asses_m.shape[0], show_labels.shape[0])

                print(t1, t2)
                showText = True
                if showT == "-NoneText":
                    showText = False
                    s_samples = samples
                    if dataset=="OoDAnimals":
                        s_samples = tsne_file['true_id'][samples]
                    if dataset=="OoDAnimals3":
                        s_samples = tsne_file['true_id'][samples]
                    np.savez(save_path, row_asses=row_asses_m, labels=new_labels, samples=s_samples)
                print("new_cost", new_cost)
                name = "\'" + dataset + "\'-" + str(grid_width) + "-" + str(flag) + "-" + type + "-" + op_type
                print(name, tt)
                new_cost = np.append(new_cost, [t1 + t2, t2], None)

                new_cost[0] = np.exp(-new_cost[0] / grid_width / grid_width)
                new_cost[1] = np.exp(-new_cost[1] / grid_width / grid_width)
                new_cost[2] = 1 - new_cost[2] / grid_width / grid_width
                new_cost[3] = 1 - new_cost[3] / grid_width / grid_width
                new_cost[4] = 1 - new_cost[4] / grid_width / grid_width
                new_cost[5] = 1 - new_cost[5] / grid_width / grid_width

                if name not in dict:
                    dict.update({name: new_cost})
                    dict2.update({name: 1})
                else:
                    dict[name] += new_cost
                    dict2[name] += 1

                name0 = name+"-"+str(tt)
                dict0.update({name0: new_cost.copy()})

                print(show_labels.max())
                if show_labels.max() < 30:
                    Optimizer.show_grid(row_asses_m, show_labels, grid_width, file_path, showText, just_save=True)
                # Optimizer.show_grid(row_asses_m, show_labels, grid_width, "E-full-NoneText.svg", showText)
                # Optimizer.show_grid(row_asses_m, show_labels, grid_width, "test"+name+".png", showText)
                # Optimizer.show_grid(row_asses_m, s_labels, grid_width)


# for dataset in ["MNIST", "STL-10", "CIFAR10", "USPS"]:
for dataset in ["MNIST"]:
# for dataset in ["MNIST", "USPS", "CIFAR10", "STL-10"]:
# for dataset in ["MNIST", "STL-10", "CIFAR10", "USPS", "Cats-vs-dogs", "Weather", "Wifi", "Indian food"]:
#     for grid_width in [20, 30, 40]:
    for grid_width in [30, 40]:
        for flag in [3, 5, 7, 'all']:
        # for flag in ['all']:
            max_tt = 20
            # if grid_width == 20:
            #     max_tt = 50
            for tt in range(max_tt):
                body()
    import numpy as np
    import pickle
    f = open("grid_result/"+dataset+"/data.pkl", 'wb+')
    pickle.dump({"dict": dict, "dict2": dict2}, f, 0)
    f.close()
    # np.savez("grid_result/"+dataset+"/data.npz", dict=dict, dict2=dict2)

# for dataset in ["MNIST"]:
#     for grid_width in [30, 40]:
#         for flag in ['all']:
#             max_tt = 50
#             for tt in range(max_tt):
#                 body()
#
# for dataset in ["MNIST"]:
#     for grid_width in [20]:
#         for flag in [3, 5, 7]:
#             max_tt = 50
#             for tt in range(max_tt):
#                 body()

# for key in dict:
#     print(key,"----", dict[key])

# f = open("grid_result/" + "MNIST" + "/data.pkl", 'rb+')
# data = pickle.load(f)
# dict = data['dict']
# dict2 =data['dict2']

for key in dict:
    t = -2
    if "base" in key:
        t = -1
    dict[key] /= dict2[key]
    print(key, "----", "%.4lf"%dict[key][0], "&", "%.4lf"%dict[key][1], "&", "%.4lf"%dict[key][2], "&", "%.4lf"%dict[key][3], "&", "%.4lf"%dict[key][4], "&", "%.4lf"%dict[key][5], "&", "%.3lf"%dict[key][t])

# for key in dict0:
#     t = -2
#     if "base" in key:
#         t = -1
#     if "-C-" in key:
#         continue
#     print(key, "----", "%.4lf"%dict0[key][3], "&", "%.4lf"%dict0[key][4])