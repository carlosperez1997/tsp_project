from typing import Dict, Tuple, List
from dataclasses import dataclass
from gurobipy import Model, GRB, tupledict, LinExpr, GurobiError

class TspInstance:
    """ An instance of the Travelling Salesman Problem. """

    def __init__(self, n_vertices, cost, x=None, y=None):
        """
        Initialize a TSP instance.

        Parameters
        ----------
        n_vertices : int
            Number of vertices in the graph.
        cost : dict[tuple[int, int], float]
            Dictionary with arcs as key, and their corresponding costs as values.
        x : list[float], optional
            X coordinates of the vertices (for plotting).
        y : list[float], optional
            Y coordinates of the vertices (for plotting).
        """
        self.n_vertices = n_vertices
        self.cost = cost
        self.x = x if x is not None else []
        self.y = y if y is not None else []
        
@dataclass
class TspSolution:
    """ Solution to an instance of the Travelling Salesman Problem.

        Attributes
        ----------
        instance: TspInstance
            Reference to the TSP instance being solved.
        tour: List[int]
            Ordered list of vertices visited (without repeating the first one).
        cost: float
            Travel cost of the tour.
    """

    instance: TspInstance
    tour: List[int]
    cost: float
    lazy_constraints_added: List[dict]
    time_limit: float
    exec_time: float
    optimal: bool


class TspIntegerBCSolver:
    """ 
    Solver for the Travelling Salesman Problem using the branch-and-cut algorithm.

    This solver applies the Dantzig-Fulkerson-Johnson (DFJ) formulation of the TSP
    and focuses on separating subtour elimination constraints (SECs) only on integer solutions.

    Attributes:
    _instance : TspInstance
        The instance of the TSP to be solved.
    _V : List[int]
        List of vertex indices in the TSP instance.
    _A : List[Tuple[int, int]]
        List of arcs (edges between pairs of vertices) in the TSP instance.
    _model : Model
        The Gurobi model for the TSP instance.
    _x : tupledict
        Gurobi variables representing whether an arc is included in the tour.
    """

    def __init__(self, instance: TspInstance, cardinality=None, time_limit=None):
        """
        Initializes the solver with a TSP instance and builds the model.

        Parameters:
        instance (TspInstance): The TSP instance to be solved.
        cardinality=None
        time_limit: Maximum Exectime in minutes
        """
        # Store the TSP instance and create lists of vertices and arcs
        self._instance = instance
        self._V = list(range(self._instance.n_vertices))
        self._A = list(self._instance.cost.keys())
        self.time_limit = time_limit
        self.lazy_constraints_added = []  # Initialize an empty list to store lazy constraints
        self.cardinality = cardinality  # Store the cardinality parameter
        self.optimal = None

        # Build the optimization model
        self.__build_model()

    def __build_model(self):
        """ 
        Builds the initial DFJ model for the TSP without subtour elimination constraints (SECs).
        """
        # Initialize the Gurobi model
        self._model = Model()

        # Add variables representing whether each arc is included in the tour
        self._x = self._model.addVars(self._A, obj=self._instance.cost, vtype=GRB.BINARY, name='x')

        # Add constraints to ensure exactly one outgoing and one incoming arc for each vertex
        self._model.addConstrs((self._x.sum(i, '*') == 1 for i in self._V), name='outgoing')
        self._model.addConstrs((self._x.sum('*', i) == 1 for i in self._V), name='incoming')

        # Conditionally add cardinality constraints if the cardinality is set to 2
        if self.cardinality == 2 or self.cardinality == 3:
            self._model.addConstrs((self._x[i, j] + self._x[j, i] <= 1 for i in self._V for j in self._V if i != j), name='cardinality two')

        # Add constraints to ensure x(i,j) + x(j,k) +x(k,i) <= 2
        if self.cardinality == 3:
            self._model.addConstrs((self._x[i, j] + self._x[j, k] + self._x[k, i]<= 2 for i in self._V for j in self._V for k in self._V if i != j and j != k and k != i), name='cardinality three')


    def __next_vertex(self, i: int) -> int:
        """
        Finds the next vertex in the current solution, starting from vertex i.

        Parameters:
        i (int): The starting vertex index.

        Returns:
        int: The index of the next vertex in the tour.
        """
        # Ensure the vertex index is valid
        assert 0 <= i < self._instance.n_vertices

        # Loop through all vertices to find the next one in the tour
        for j in self._V:
            if i != j:
                try:
                    # Inside a callback, retrieve the solution using cbGetSolution
                    val = self._model.cbGetSolution(self._x[i, j])
                except:
                    # Otherwise, access the variable's value directly
                    val = self._x[i, j].X

                # If the arc is part of the solution, return the next vertex
                if val > 0.5:
                    return j

        # If no successor is found, raise an error
        assert False, f"Vertex {j} has no successor?!"

    def __tour_staring_at(self, i: int) -> List[int]:
        """
        Constructs a tour (or subtour) starting at vertex i based on the current solution.

        Parameters:
        i (int): The starting vertex index.

        Returns:
        List[int]: A list of vertex indices representing the tour.
        """
        # Ensure the vertex index is valid
        assert 0 <= i < self._instance.n_vertices

        # Initialize the tour with the starting vertex
        tour = [i]
        current = self.__next_vertex(i)

        # Continue adding vertices to the tour until it loops back to the start
        while current != i:
            tour.append(current)
            current = self.__next_vertex(current)

        return tour

    def __add_sec_for(self, S: List[int]) -> None:
        """
        Adds a subtour elimination constraint (SEC) for the vertices in set S.

        Parameters:
        S (List[int]): A list of vertex indices forming a subtour.
        """
        # Ensure the subtour has a valid number of vertices
        assert 0 < len(S) < self._instance.n_vertices
   
        # Record the constraint information
        constraint_info = {
            'subtour': S,
            'constraint': sum(self._x[i, j] for (i, j) in self._A if i in S and j not in S) >= 1
        }

        self.lazy_constraints_added.append(constraint_info['subtour'])
        
        # Add the actual lazy constraint to the model
        self._model.cbLazy(sum(self._x[i, j] for (i, j) in self._A if i in S and j not in S) >= 1)

    def __separate(self, where: int) -> None:
        """ Separates eventual violated SECs. """

        if where != GRB.Callback.MIPSOL:
            return
        
        # Set of vertices not yet placed in a tour
        remaining = set(self._V)

        while len(remaining) > 0:
            # Get any vertex in set `remaining`
            current = next(iter(remaining))

            # Get the subtour which includes the vertex
            subtour = self.__tour_staring_at(current)

            # If the tour visits all vertices, nothing to do
            if len(subtour) == self._instance.n_vertices:
                return

            # Otherwise, it's a subtour => Add a SEC
            self.__add_sec_for(subtour)

            # Update set `remaining`
            remaining -= set(subtour)

    def solve(self) -> TspSolution:
        """ Solves the TSP DFJ model and returns the solution.
        
            It throws a RuntimeError if Gurobi cannot solve the model to optimality.

            Returns
            -------
            A TspSolution object with details about the solution.
        """

        # We must set this parameter to use callbacks
        self._model.setParam(GRB.Param.LazyConstraints, 1)

        # Gurobi always passes two parameters to the callback: `model` and `where`
        self._model.optimize(lambda _, where: self.__separate(where))

        if self.time_limit is not None:
            self._model.setParam('TimeLimit', self.time_limit*60)

        if self._model.Status != GRB.OPTIMAL:
            #raise RuntimeError(f"Could not find the optimal solution. Gurobi status = {self._model.Status}")
            pass
        
        #if self._model.Status in [GRB.INFEASIBLE, GRB.INF_OR_UNBD]:
        #    print(f"Warning: The model is infeasible or unbounded. Gurobi status = {self._model.Status}")

        if self._model.Status not in [GRB.OPTIMAL, GRB.SUBOPTIMAL]:
            raise RuntimeError(f"Could not find the optimal solution. Gurobi status = {self._model.Status}")
        
        return TspSolution(
            instance=self._instance,
            tour=self.__tour_staring_at(0),
            optimal=self._model.Status != GRB.OPTIMAL,
            cost=self._model.ObjVal,
            lazy_constraints_added=self.lazy_constraints_added,
            time_limit=self.time_limit*60 if self.time_limit is not None else None,
            exec_time=self._model.Runtime
        )