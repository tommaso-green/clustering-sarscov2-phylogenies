import math
import os
import argparse
import numpy as np
from Bio import SeqIO
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib import rc

rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)


def get_dict_of_nations(filename):
    """Using the file of all sequences, get for each sequence its nation"""
    sequences = list(SeqIO.parse(filename, "fasta"))
    offset = 10
    sequences_by_nation = {}
    nations = []
    for sequence in sequences:
        code_start = sequence.description.find("ISL")  # find index where ISL is
        sequence_code = sequence.description[code_start:code_start + offset]  # get substring ISL_XXXXXX
        sequence_nation = sequence.description.split("/")[1]
        if sequence_nation in ['Wuhan', 'Wuhan-Hu-1', 'Shanghai', 'NanChang', 'Fujian', 'Guangdong']:
            sequence_nation = "China"
        if sequence_nation in ['England', 'Wales', 'Scotland']:
            sequence_nation = "UK"
        if sequence_nation not in nations:
            nations.append(sequence_nation)
        sequences_by_nation[sequence_code] = sequence_nation
    print("Found the following %d nations:" % len(nations))
    print(nations)
    return nations, sequences_by_nation


def compute_entropy(clusters):
    entropy_vector = np.empty(
        len(clusters))  # for cluster i, entropy_vector[i] has its entropy (cluster -1 is in entropy_vector[0])
    no_of_points = sum([len(clusters[x]) for x in clusters])  # getting size of cluster
    k = 0
    for cluster_id in clusters.keys():  # for each cluster
        if cluster_id != -1:
            size_of_cluster = len(clusters[cluster_id])
            nation_list = [element[1] for element in
                           clusters[cluster_id]]  # getting labels of nations present in cluster
            nations_in_cluster = {}
            for item in nation_list:
                nations_in_cluster[item] = nation_list.count(
                    item)  # getting a dict that stores nation:no_of_occurences_of_nation
            entropy = 0
            for nation in nations_in_cluster.keys():
                occurence_of_nation = nations_in_cluster[nation]
                p = occurence_of_nation / size_of_cluster
                entropy += -p * math.log(p, 2)  # entropy formula
            entropy_vector[k] = entropy  # storing value
            k += 1
        else:  # clusters with label -1 are actually singleton clusters with 0 entropy
            entropy_vector[k] = 0  # storing value
            k += 1
    total_entropy = 0
    i = 0
    for cluster_id in clusters.keys():
        total_entropy += (len(clusters[cluster_id]) / no_of_points) * entropy_vector[
            i]  # formula for total entropy of clustering
        i += 1
    return entropy_vector, total_entropy


def compute_entropy_by_class(clusters, nations):
    entropy_by_nation = dict.fromkeys(nations, 0)
    total_occurrences_by_nation = dict.fromkeys(nations, 0)
    for key in list(clusters):
        clusters[key] = [x[1] for x in clusters[key]]  # removing sequence code (useless for entropy)
    for nation in nations:
        sum = 0
        for key in list(clusters):
            sum += clusters[key].count(nation)
        total_occurrences_by_nation[nation] = sum  # dictionary nation:occurrence_of_nation_in_cluster
    for nation in nations:
        entropy = 0
        for key in list(clusters):
            if key == -1:  # formula for singleton clusters
                no_of_singletons = clusters[-1].count(nation)  # singletons with nation
                p = 1 / (total_occurrences_by_nation[nation])
                entropy_of_singleton = -p * math.log(p, 2)
                entropy += no_of_singletons * entropy_of_singleton  # all singletons have same entropy
            else:
                occurence_in_cluster = clusters[key].count(nation)
                if occurence_in_cluster > 0:  # avoid Math Exception with log(0)
                    p = occurence_in_cluster / total_occurrences_by_nation[nation]
                    entropy += -p * math.log(p, 2)
        entropy_by_nation[nation] = entropy
    return entropy_by_nation


def print_nations_in_clusters(clusters):
    keys = list(clusters.keys())
    if -1 in keys:
        keys.remove(-1)
    for cluster_id in keys:
        print("CLUSTER NO. " + str(cluster_id) + ": " + str(sorted([x[1] for x in clusters[cluster_id]])))
    if -1 in clusters.keys():
        print("SINGLETON CLUSTERS")
        max_cluster_id = max(clusters.keys())
        for element in clusters[-1]:
            max_cluster_id += 1
            print("CLUSTER NO. " + str(max_cluster_id) + ": " + element[1])


def main():
    """ Parsing Arguments """
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--talign', required=True, type=float)
    parser.add_argument('-td2', '--td2', required=True, type=float)
    parser.add_argument('-td2*', '--td2star', required=True, type=float)
    parser.add_argument('-m', '--method', required=True, type=str)
    parser.add_argument('-v', '--verbose', action='store_true', required=False)
    parser.add_argument('-tf', '--threshold_free', required=False, type=str)
    args = parser.parse_args()

    # Getting list of nations (labels) and for each sequence its nation
    nations, sequences_by_nation = get_dict_of_nations("allsequences.fasta")
    nations = sorted(nations)
    directory = "tree_files/converted/"
    file_dict = {"Pairwise": directory + "pairwise_converted.nwk",
                 "ClustalW": directory + "multiplealignment_converted.nwk", "D2": directory + "outtreeD2_converted.nwk",
                 "D2Star": directory + "outtreeD2Star_converted.nwk"}

    # PARAMETERS
    alignment_threshold = args.talign
    d2_threshold = args.td2
    d2star_threshold = args.td2star
    method = args.method
    verbose = args.verbose
    threshold_free = args.threshold_free

    print("~" * 30)
    for key in list(file_dict):
        print("Starting clustering for " + key)
        command = "python TreeCluster.py -i "
        command += file_dict[key]
        if key in ["Pairwise", "ClustalW"]:
            command += " -t " + str(alignment_threshold)
        elif key == "D2":
            command += " -t " + str(d2_threshold)
        else:
            command += " -t " + str(d2star_threshold)
        command += " -m " + method
        outfile = "outputs/" + key + ".txt"
        command += " -o " + outfile
        if threshold_free:
            command += " -tf argmax_clusters"
        if verbose:
            command += " -v"
        print(command)
        os.system(command)
        clustering = open(outfile, "r+").readlines()[1:]  # not first line because it's SequenceName, Clustering
        clustering = [x.rstrip() for x in clustering]  # remove \n at the end of each line
        clusters = {}  # this dictionary will have as key the cluster_id, as value the list of all elements of such
        # cluster
        for line in clustering:  # reading the file
            sequence, cluster_id = line.split("	")  # getting ISL_XXXXXX cluster_id
            cluster_id = int(cluster_id)  # cluster_id is an integer
            nation = sequences_by_nation[sequence]  # get nation of considered sequence
            # print(sequence, cluster_id, nation)
            if cluster_id not in clusters.keys():  # if a new cluster_id is considered
                clusters[cluster_id] = [(sequence, nation)]  # initialise its value
            else:
                clusters[cluster_id].append((sequence, nation))  # otherwise append new tuple
        if -1 not in clusters.keys():  # if there are no singleton clusters
            no_of_clusters = len(clusters.keys())
        else:  # otherwise this is the formula
            no_of_clusters = len(clusters.keys()) - 1 + (len(clusters[-1]))
        print(">>>>NUMBER OF CLUSTERS FOUND: " + str(no_of_clusters))
        print_nations_in_clusters(clusters)
        cluster_names = np.array(list(clusters.keys()))  # getting cluster entropy for table
        entropy_vector, total_entropy = compute_entropy(clusters)  # computing entropy for table
        cluster_sizes = np.array([len(clusters[key]) for key in clusters.keys()])  # list of sizes of each cluster
        data = Table([cluster_names, cluster_sizes, entropy_vector], names=['Cluster', 'Size of Cluster', 'Entropy'])
        print(data)
        print(">>>>TOTAL ENTROPY: " + str(total_entropy))

        ### PLOTTING ENTROPY FOR CLUSTERING
        x_pos = [i for i, _ in enumerate(list(clusters))]
        plt.bar(x_pos, entropy_vector, color='red')
        plt.xlabel("Clusters")
        plt.ylabel("Entropy")
        plt.axhline(y=math.log(len(nations), 2), color='gray', linestyle='--')  # max_entropy
        plt.annotate("Max entropy", xy=(1, math.log(len(nations), 2) - 0.1))
        if key == "D2":
            print("hey")
            plt.title("Entropy of Clusters for " + "$EP_{2}$" + " Tree Clustering")
        elif key == "D2Star":
            plt.title("Entropy of Clusters for " + "$EP_{2}^{*}$" + " Tree Clustering")
        else:
            plt.title("Entropy of Clusters for " + key + " Tree Clustering")
        plt.xticks(x_pos, list(clusters))
        plt.show()
        ###

        entropy_by_nation = compute_entropy_by_class(clusters, nations)
        nation_array = np.array(list(entropy_by_nation))
        entropy_nation_vector = np.array(list(entropy_by_nation.values()))
        data2 = Table([nation_array, entropy_nation_vector], names=['Nation', 'Entropy of nation'])
        print(data2)

        # PLOTTING ENTROPY FOR NATIONS
        x_pos = [i for i, _ in enumerate(nation_array)]
        plt.bar(x_pos, entropy_nation_vector, color='blue')
        plt.xlabel("Nations")
        plt.ylabel("Entropy of Nation")
        plt.axhline(y=math.log(no_of_clusters, 2), color='gray', linestyle='--')
        plt.annotate("Max entropy", xy=(1, math.log(no_of_clusters, 2) + 0.03))
        if key == "D2":
            plt.title("Entropy of Nations for " + "$EP_{2}$" + " Tree Clustering")
        elif key == "D2Star":
            plt.title("Entropy of Nations for " + "$EP_{2}^{*}$" + " Tree Clustering")
        else:
            plt.title(r"Entropy of Nations for " + key + " Tree Clustering")
        plt.xticks(x_pos, nation_array)
        plt.xticks(rotation=90)
        plt.show()
        ###

        print("~" * 30)


main()
