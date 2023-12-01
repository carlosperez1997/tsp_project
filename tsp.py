import numpy as np
import re
import math

from helper import *

# Define data load class
class TspDataLoader:
    """ Class to load and process data for TSP instances. """

    def __init__(self, dir='data'):
        self.dir = dir

    def load_problem(self, problem):
        """ Load a TSP problem and return its vertices count and cost matrix. """
        file_name = "data/"+problem

        # Read file content
        with open(file_name, 'r') as file:
            content = file.read()

        problem_type = problem.split('.')[-1]

        n = int(re.search(r'\d+', problem).group())
        # 64 is actually 65 
        if n == 64:
            n = 65

        # Process the file content
        data_lines = content.split('\n')

        try:
            start_index = data_lines.index('EDGE_WEIGHT_SECTION')
        except:
            start_index = data_lines.index('NODE_COORD_SECTION')
            
        # Extract the data lines after 'EDGE_WEIGHT_SECTION'
        edge_weight_lines = data_lines[start_index + 1:][:-2]

        # Parse the distance matrix (lines are not rows)
        if problem_type == 'atsp':
            distance_matrix = [list(map(int, line.split())) for line in edge_weight_lines if line.strip()]

            def flatten(l):
                return [item for sublist in l for item in sublist]

            distance_list = flatten(distance_matrix)
            my_matrix = np.array(distance_list).reshape((n, n))
            
        else:
            # Avoid index
            try:
                ys = [ float(line.split(' ')[3])  for line in edge_weight_lines if line.strip()]
                xs = [ float(line.split(' ')[2])  for line in edge_weight_lines if line.strip()]
            except:
                ys = [ float(line.split(' ')[2])  for line in edge_weight_lines if line.strip()]
                xs = [ float(line.split(' ')[1])  for line in edge_weight_lines if line.strip()]

        for line in data_lines:
            if 'EDGE_WEIGHT_TYPE' in line:
                weight_type = line.split(':')[1].strip()
        print(weight_type)

        V = range(n)
        if problem_type == 'atsp':
            cost = {
                (i, j): my_matrix[i,j]
                for i in V for j in V if i != j
            }
        else:
            if weight_type == 'GEO':
                cost = {
                    (i, j): haversine_distance((xs[i], ys[i]), (xs[j], ys[j]))
                    for i in V for j in V if i != j
                }
            else:
                cost = {
                    (i, j): math.sqrt((xs[i] - xs[j]) ** 2 + (ys[i] - ys[j]) ** 2)
                    for i in V for j in V if i != j
                }

        return n, cost