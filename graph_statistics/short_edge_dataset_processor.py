import pandas as pd
import matplotlib as mpl
import numpy as np
from sklearn import svm
from sklearn.model_selection import train_test_split
from mpl_toolkits.mplot3d import Axes3D
import os
import shutil

mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def clean_and_add_columns(short_data):
    MIN_LENGTH = 5000
    short_data = short_data[short_data["Length"] < MIN_LENGTH]
    short_data["MinInter"] = np.minimum(short_data["LeftInter"], short_data["RightInter"])
    short_data["AvCov"] = (short_data["LeftCov"] + short_data["RightCov"]) / 2
    short_data["CovRatio"] = short_data["Coverage"] / short_data["AvCov"]
    short_data["CovSquare"] = short_data["Coverage"] * (short_data["LeftCov"] + short_data["RightCov"])
    short_data["CovCovLen"] = (short_data["CovSquare"] * short_data["Length"]) / 1000000
    short_data["Score"] = np.minimum(short_data["LeftInter"], short_data["RightInter"]) / short_data["CovCovLen"]
    short_data["ContIndex"] = short_data["MinInter"] / short_data["Barcodes"]
    print("Small covcovlen: {}".format(len([x for x in short_data["CovCovLen"] if x < 0.01])))
    print("Overall: {}".format(len(short_data)))
    short_data = short_data[short_data["CovRatio"] < 3]
    random_short_data = short_data[short_data["Correct"] == 0]
    correct_short_data = short_data[short_data["Correct"] == 1]
    random_sample = random_short_data.sample(n=50000)
    short_data = pd.concat([correct_short_data, random_sample])
    print(len(short_data))
    return short_data


def prepare_dataset(short_data):
    MIN_LENGTH = 5000
    print(len(short_data))
    short_data = short_data[short_data["Barcodes"] == short_data["LeftInter"]]
    short_data = short_data[short_data["Barcodes"] == short_data["RightInter"]]
    print(len(short_data))
    short_data = short_data[short_data["Length"] < MIN_LENGTH]

    short_data["MinInter"] = min(short_data["LeftInter"], short_data["RightInter"])

    correct_data = short_data[short_data["Correct"] == 1]
    random_data = short_data[short_data["Correct"] == 0]
    random_sample = np.random.choice(len(random_data), len(random_data) / 10)
    print(len(correct_data))
    print(len(random_sample))


def draw_violin_plots(short_data, output_path):
    print("Drawing violin plots")
    short_data["NormBySquare"] = short_data["MinInter"] / short_data["CovSquare"]
    param_names = ["Score", "MinInter", "NormBySquare"]
    violin_path = os.path.join(output_path, "violin_plots")
    if not os.path.exists(violin_path):
        os.mkdir(violin_path)
    else:
        shutil.rmtree(violin_path)
    for name in param_names:
        draw_param_violin_plot(short_data, name, violin_path)


def draw_param_violin_plot(short_data, param_name, output_path):
    sns.violinplot(x="Correct", y=param_name, data=short_data)
    plt.ylim(0, np.percentile(short_data[param_name], 95))
    plt.savefig(os.path.join(output_path, param_name))
    plt.clf()


def draw_score_to_params(short_data, output_path):
    print("Drawing 2d score plots")
    correct = short_data[short_data["Correct"] == 1]
    random = short_data[short_data["Correct"] == 0]
    param_names = ["CovRatio", "Length", "Coverage", "CovSquare", "CovCovLen"]
    for name in param_names:
        draw_score_to_param(data=short_data, param_name=name, output_path=output_path)


def containment_index_criteria(element):
    series = element[1]
    return series["Length"] < 50 or series["ContIndex"] > 0.006


def covcovlen_criteria(element):
    series = element[1]
    score = series["MinInter"] / series["CovCovLen"]
    return series["CovCovLen"] < 0.5 or score >= 10.0


def cov_and_len_criteria(element):
    series = element[1]
    score_cov = series["MinInter"] / series["CovSquare"]
    score_len = series["MinInter"] / series["Length"]
    return series["Length"] < 50 or (score_cov >= 0.00005 and score_len > 0.01)


def test_criteria(data, criteria, name):
    print("Testing " + name)
    correct = data[data["Correct"] == 1]
    random = data[data["Correct"] == 0]
    print("Correct: {}".format(len(correct)))
    print("Random: {}".format(len(random)))
    tp = len([x for x in correct.iterrows() if criteria(x)])
    fp = len([x for x in random.iterrows() if criteria(x)])
    sensitivity = tp / (len(correct))
    specificity = 1 - fp / len(random)
    print("Sensitivity: {}".format(sensitivity))
    print("Specificity: {}".format(specificity))


def test_criterias(data):
    criterias = [(covcovlen_criteria, "covcovlen"),
                 (cov_and_len_criteria, "cov_and_len"),
                 (containment_index_criteria, "containment index")]
    for func, name in criterias:
        test_criteria(data=data, criteria=func, name=name)


def draw_score_to_param(data, param_name, output_path):
    sns.lmplot(x=param_name, y="MinInter", hue="Correct", data=data, fit_reg=False)
    plt.xlabel(param_name)
    plt.ylabel("Min intersection")
    plt.xlim(0, np.percentile(data[param_name], 99))
    # plt.ylim(0, np.percentile(data["Score"], 99))
    if param_name == "CovCovLen":
        plt.xlim(0, 2)
    plt.ylim(0, 10)
    plt.savefig(plt.savefig(os.path.join(output_path, param_name)))
    plt.clf()


def launch_svm(data, output_path):
    y = data["Correct"]
    param_name = "CovSquare"
    X = data[[param_name, "Score"]]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    coef, intercept = train_svm(X_train=X_train, y_train=y_train)
    draw_svm_plot(data=data, coef=coef, intercept=intercept, param_name=param_name, output_path=output_path)


def train_svm(X_train, y_train):
    clf = svm.LinearSVC(class_weight={0: 5})
    clf.fit(X_train, y_train)
    print(clf.coef_)
    return clf.coef_, clf.intercept_


def get_plot_function(coef, intercept):
    print(coef)
    a, b = coef[0][0], coef[0][1]
    return lambda x: -(a * x + intercept) / b


def draw_svm_plot(data, coef, intercept, param_name, output_path):
    sns.lmplot(x=param_name, y="Score", hue="Correct", data=data, fit_reg=False)
    plt.xlabel(param_name)
    plt.ylabel("Min intersection")
    xlimits = (0, np.percentile(data[param_name], 99))
    plt.xlim(xlimits)
    plt.ylim(0, np.percentile(data["Score"], 99))
    plot_function = get_plot_function(coef, intercept)
    plt.plot([xlimits[0], xlimits[1]], [plot_function(x) for x in [xlimits[0], xlimits[1]]], 'k-', lw=3, color='r')
    plt.savefig(plt.savefig(os.path.join(output_path, param_name + "_svm")))
    plt.clf()


def draw_3d_plot(short_data, output_path):
    correct_data = short_data[short_data["Correct"] == 1]
    random_data = short_data[short_data["Correct"] == 0]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    cor_x, cor_y, cor_z = get_3d_data(correct_data)
    ax.scatter(cor_x, cor_y, cor_z, c='g')
    rand_x, rand_y, rand_z = get_3d_data(random_data)
    kwargs = {'xs': rand_x, 'ys': rand_y, 'zs': rand_z, 'c': 'r', 'alpha': 0.15}
    ax.scatter(**kwargs)

    ax.set_xlabel("Length")
    ax.set_ylabel("Coverage")
    ax.set_zlabel("Score")
    ax.set_zlim(0.0, 0.1)

    plt.savefig(os.path.join(output_path, "cov_and_len_to_score"))


def get_3d_data(data):
    scores = data["Score"]
    lengths = data["Length"]
    coverages = data["Coverage"]
    ratios = data["CovRatio"]
    return lengths, coverages, scores
