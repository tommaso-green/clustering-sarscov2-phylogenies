import os


def tree_convert():
    """ This program converts Newick Tree in tree_files into a single line file for TreeCluster.py"""
    your_directory = "tree_files"
    tree_files = [f for f in os.listdir(your_directory) if os.path.isfile(os.path.join(your_directory, f))]
    for filename in tree_files:
        f = open("/".join([your_directory, filename]), "r+")
        lines = [x.rstrip() for x in f.readlines()]
        tree_string = "".join(lines)
        f.close()
        filename_no_ext = filename.split(".")[0]
        out = open("/".join([your_directory, "converted", filename_no_ext]) + "_converted.nwk", "w+")
        out.write(tree_string)
        out.close()


tree_convert()
