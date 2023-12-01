import numpy as np
import os
import re
import math

from typing import Dict, Tuple, List
from math import sqrt
from dataclasses import dataclass
from gurobipy import Model, GRB, tupledict
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

import numpy as np
import os
import re

def calculate_distance(node1, node2):
    x1, y1 = node1
    x2, y2 = node2
    return math.ceil(math.sqrt((x2 - x1)**2 + (y2 - y1)**2))


def haversine_distance(coord1, coord2):
    # Convert latitude and longitude from degrees to radians
    lat1, lon1 = math.radians(coord1[0]), math.radians(coord1[1])
    lat2, lon2 = math.radians(coord2[0]), math.radians(coord2[1])

    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat / 2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    # Radius of the Earth in kilometers (you can adjust this value)
    R = 6371.0

    # Calculate the distance
    distance = R * c

    return distance


def get_opt_sol(problem):
    problem_short = problem.split('.')[0]
    if 'atsp' in problem:
        file = 'data/opt_atsp.txt'
    else:
        file = 'data/opt_tsp.txt'
    
    with open(file, 'r') as file:
        content = file.read()
    
    lines = content.split('\n')
    for line in lines:
        if line.split(':')[0].strip() == problem_short:
            return int(line.split(':')[1].strip())
    return None


def lazy_constraints_size_count(ctrs, max_level=3):
    constraints_count = {}
    for ctr in ctrs:
        n = len(ctr)
        if n > max_level:
            n = 'more than ' + str(max_level)
        if n in constraints_count:
            constraints_count[n] += 1
        else:
            constraints_count[n] = 1

    return constraints_count


def get_solution_info(solution):
    return {'cost': solution.cost,
            'tour': solution.tour,
            'lazy_constraints': solution.lazy_constraints_added,
            'lazy_constraints_count': lazy_constraints_size_count(solution.lazy_constraints_added),
            'exec_time': solution.exec_time,
            'time_limit': solution.time_limit
            }