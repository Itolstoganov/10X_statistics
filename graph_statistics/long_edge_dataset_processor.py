import pandas as pd
import matplotlib as mpl
import numpy as np
from sklearn import svm
from sklearn.model_selection import train_test_split
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter
import os

mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def prepare_dataset(data):
    print(len(data))
    data["MinSize"] = np.minimum(data["LeftSize"], data["RightSize"])
    data["ContIndex"] = data["Intersection"] / data["MinSize"]
    data["AvCov"] = (data["LeftCov"] + data["RightCov"]) / 2
    data["MinCov"] = np.minimum(data["LeftCov"], data["RightCov"])
    data["MaxCov"] = np.maximum(data["LeftCov"], data["RightCov"])
    data["SCovContIndex"] = data["ContIndex"] / np.sqrt(data["MinSize"]) * 100
    data["InterToCov"] = data["Intersection"] / data["AvCov"]
    data["InterToMinSizeSquare"] = data["Intersection"] / (data["MinSize"] * data["MinSize"])
    data["DistanceBin"] = data["Distance"] // 500
    correct_data = data[data["Correct"] == 1]
    print("{} correct pairs".format(len(correct_data)))
    print("{} close pairs(10000)".format(len([x for x in correct_data["Distance"] if x < 10000])))
    print("{} close pairs(5000)".format(len([x for x in correct_data["Distance"] if x < 5000])))
    # size_counter = Counter(data["Genome"])
    # data = data[data.apply(lambda x: size_counter[x["Genome"]] > 20, axis=1)]
    # genome_dict = {0: "Distant", 1: "1", 4: "2", 5: "3", 6: "4"}
    # genome_list = list(data["Genome"])
    # print(genome_list[:10])
    # data["Reference"] = [genome_dict[x] for x in genome_list]
    # print(size_counter)
    return data


def get_stats(data):
    correct_data = data[data["Correct"] == 1]
    zero_intersection = correct_data[correct_data["Intersection"] == 0]
    print("Zero intersection: {}".format(len(zero_intersection)))
    small_containment = correct_data[correct_data["ContIndex"] < 0.03]
    print("Small contaiment: {}".format(len(small_containment)))
    small_containment_close = small_containment[small_containment["Distance"] < 1000]
    print("Small containment close: {}".format(len(small_containment_close)))
    # distance_counter = Counter(zero_intersection["DistanceBin"])
    # print(distance_counter)
    print(small_containment_close["Genome"])
    random_data = data[data["Correct"] == 0]
    large_random = random_data[random_data["ContIndex"] > 0.05]
    print(len(large_random))
    print(np.percentile(data["LeftSize"], 50))

def draw_2d_plots(data, output_path):
    x_names = ["AvCov", "MinSize", "Distance"]
    y_names = ["ContIndex", "Intersection", "InterToCov"]
    for x_name in x_names:
        for y_name in y_names:
            draw_score_to_param(data=data, param_name=x_name, score_name=y_name, output_path=output_path)

    draw_score_to_param(data=data, param_name="AvCov", score_name="ContIndex", output_path=output_path, hue_name="Genome")


def draw_cont_index_histogram(data, output_path):
    correct_data = data[data["Correct"] == 1]
    random_data = data[data["Correct"] == 0]
    plt.hist(correct_data["ContIndex"], bins=50, normed=1, cumulative=True)
    plt.savefig(os.path.join(output_path, "cont_index_hist"))
    plt.clf()


def draw_violin_plot(data, output_path):
    correct_data = data[data["Correct"] == 1]
    # sns.lvplot(x="DistanceBin", y="ContIndex", data=correct_data)
    sns.lvplot(x="Distance", data=correct_data)
    plt.ylim(0, np.percentile(correct_data["ContIndex"], 95))
    plt.xlim(0, np.percentile(correct_data["Distance"], 95))
    # print(np.percentile(correct_data["ContIndex"], 20))
    plt.savefig(os.path.join(output_path, "distance_lvplot"))
    plt.clf()


def draw_score_to_param(data, param_name, score_name, output_path, hue_name="Genome"):
    sns.lmplot(x=param_name, y=score_name, hue=hue_name, data=data, fit_reg=False)
    plt.xlabel(param_name)
    plt.ylabel(score_name)
    plt.xlim(0, np.percentile(data[param_name], 99))
    if param_name == "Distance":
        plt.xlim(0, 20000)
    plt.ylim(0, np.percentile(data[score_name], 99))
    plt.savefig(plt.savefig(os.path.join(output_path, score_name + "_to_" + param_name)))
    plt.clf()