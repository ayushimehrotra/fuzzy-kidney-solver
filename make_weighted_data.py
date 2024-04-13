import random

def start():
    for i in range(200):
        output = open("choice_data/choice_{}.txt".format(i), "w")
        input = open("data/data{}.txt".format(i), "r")
        lines = input.readlines()

        input_lines = [line for line in lines if len(line.strip()) > 0]
        n_digraph_edges = int(input_lines[0].split()[1])
        digraph_lines = input_lines[1:n_digraph_edges + 2]
        output.write(str(input_lines[0].split()[0]) + " " + str(input_lines[0].split()[1]) +'\n')
        for j in digraph_lines:
            tokens = [x for x in j.split()]
            output.write(str(tokens[0]) + " " + str(tokens[1]) + " " + str(tokens[2]) + " 1.0 1.0\n")

        if len(input_lines) > len(digraph_lines):
            ndd_lines = input_lines[n_digraph_edges + 2:]
            ndd_count, edge_count = [int(x) for x in ndd_lines[0].split()]
            ndd_lines = ndd_lines[1:]
            output.write(str(ndd_count) + " " + str(edge_count) +'\n')
            for j in ndd_lines:
                tokens = [x for x in j.split()]
                output.write(str(tokens[0]) + " " + str(tokens[1]) + " " + str(tokens[2]) + " 1.0 1.0\n")


if __name__ == "__main__":
    start()