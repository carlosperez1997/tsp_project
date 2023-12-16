from typing import Dict, Tuple, List
from math import sqrt
from dataclasses import dataclass
from gurobipy import Model, GRB, tupledict, LinExpr, GurobiError
import matplotlib.pyplot as plt

from graph_tool import Graph, EdgePropertyMap
from graph_tool.generation import complete_graph
from graph_tool.flow import boykov_kolmogorov_max_flow, min_st_cut

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

    def plot(self) -> Tuple[Figure, plt.Axes]:
        """ Plots the solution of the TSP.
        
            Vertices are black dots, the solution is represented by red lines.

            Returns
            -------
            A tuple with the figure and axes object from matplotlib.    
        """

        xy = [(self.instance.x[i], self.instance.y[i]) for i in self.tour + [self.tour[0]]]

        fig, ax = plt.subplots(figsize=(10, 10))
        ax.scatter(self.instance.x, self.instance.y, color='black')
        
        for pt1, pt2 in zip(xy[1:], xy[:-1]):
            x1, y1 = pt1
            x2, y2 = pt2
            ax.plot([x1, x2], [y1, y2], color='red')

        ax.set_title(f"TSP Solution. Num vertices: {self.instance.n_vertices}. Tour cost: {self.cost:.2f}.")
        ax.set_axis_off()

        return fig, ax

class TspFractionalBCSolver:
    """ Solver for the Travelling Salesman Problem.
    
        It uses the branch-and-cut algorithm applied to the "DFJ" formulation.
        Subtour elimination constraints are separated on integer and
        fractional solutions alike.
    """

    _EPS: float = 1e-6

    _instance: TspInstance
    _V: List[int]
    _A: List[Tuple[int, int]]
    _g: Graph
    _cap: EdgePropertyMap
    _model: Model
    _x: tupledict

    def __init__(self, instance: TspInstance, cardinality=None, time_limit=None):
        """ Initialises the solver and builds the model.
        
            Parameters
            ----------
            instance: TspInstance
                The TSP instance.
        """

        self._instance = instance
        self._V = list(range(self._instance.n_vertices))
        self._A = list(self._instance.cost.keys())
        self.time_limit = time_limit
        self.cardinality = cardinality  # Store the cardinality parameter
        self.__build_graph()
        self.__build_model()
        

    def __build_graph(self) -> None:
        """ Builds the graph for cut separation. """
        self._g = complete_graph(N=self._instance.n_vertices, self_loops=False, directed=True)
        self._cap = self._g.new_edge_property(value_type='double')
        self._g.edge_properties['cap'] = self._cap       

    def __build_model(self) -> None:
        """ Builds the DFJ model without SECs. """
         
        self._model = Model()
        self._x = self._model.addVars(self._A, obj=self._instance.cost, vtype=GRB.BINARY, name='x')
        self._model.addConstrs((self._x.sum(i, '*') == 1 for i in self._V), name='outgoing')
        self._model.addConstrs((self._x.sum('*', i) == 1 for i in self._V), name='incoming')
                # Conditionally add cardinality constraints if the cardinality is set to 2
        if self.cardinality == 2 or self.cardinality == 3:
            self._model.addConstrs((self._x[i, j] + self._x[j, i] <= 1 for i in self._V for j in self._V if i != j), name='cardinality two')

        # Add constraints to ensure x(i,j) + x(j,k) +x(k,i) <= 2
        if self.cardinality == 3:
            self._model.addConstrs((self._x[i, j] + self._x[j, k] + self._x[k, i]<= 2 for i in self._V for j in self._V for k in self._V if i != j and j != k and k != i), name='cardinality three')

    def __next_vertex(self, i: int) -> int:
        """ Gets the vertex visited immediately after i in the current TSP solution.

            This is the (only) vertex j such that x[i,j] == 1. This function works both
            when inside a callback (when the value of x must be retrieved via cbGetSolution)
            and after the optimisation is over (when the value of X can be accessed via X).
        """

        assert 0 <= i < self._instance.n_vertices

        for j in self._V:
            if i != j:
                try:
                    # If inside a call-back, we use cbGetSolution
                    val = self._model.cbGetSolution(self._x[i, j])
                except:
                    # Otherwise we use .X
                    val = self._x[i, j].X

                if val > 0.5:
                    return j

        assert False, f"Vertex {j} has no successor?!"

    def __tour_staring_at(self, i: int) -> List[int]:
        """ Gives the (sub)tour starting at a given vertex i in the current TSP solution. """

        assert 0 <= i < self._instance.n_vertices

        tour = [i]
        current = self.__next_vertex(i)

        while current != i:
            tour.append(current)
            current = self.__next_vertex(current)

        return tour

    def __add_sec_for(self, S: List[int]) -> None:
        """ Adds a Subtour Elimination Constraint for set S. """

        assert 0 < len(S) < self._instance.n_vertices

        self._model.cbLazy(sum(self._x[i, j] for (i, j) in self._A if i in S and j not in S) >= 1)

    def __set_capacity(self) -> None:
        """ Sets graph capacity based on current solution. """

        for e in self._g.edges():
            i, j = e.source(), e.target()

            try: # If we are in MIPSOL
                xval = self._model.cbGetSolution(self._x[i,j])
            except: # If we are in MIPNODE
                xval = self._model.cbGetNodeRel(self._x[i,j])

            self._cap[e] = xval

    def __separate(self, where: int) -> None:
        """ Separates eventual violated SECs. """

        if where not in [GRB.Callback.MIPSOL, GRB.Callback.MIPNODE]:
            return

        # In MIPNODE, we must ensure that we have the optimal solution to the
        # linear relaxation at that node. Otherwise we are not "authorised"
        # to call cbGetNodeRel. See:
        # https://www.gurobi.com/documentation/9.1/refman/py_model_cbgetnoderel.html
        if where == GRB.Callback.MIPNODE and self._model.cbGet(GRB.Callback.MIPNODE_STATUS) != GRB.OPTIMAL:
            return
        
        self.__set_capacity()

        source = self._g.vertex(0)
        already_added = set()

        print('Separation round!')

        for i in range(1, self._instance.n_vertices):
            if i in already_added:
                continue

            target = self._g.vertex(i)
            residual = boykov_kolmogorov_max_flow(
                g=self._g, source=source, target=target,
                capacity=self._cap
            )
            
            flow = residual.copy()
            flow.a = self._cap.a - residual.a

            total_flow = sum(flow[a] for a in target.in_edges())

            if total_flow < 1 - self._EPS: # Avoid numerical issues
                print(f"\tFound a subset which violate a SEC (flow = {total_flow:.3f} < 1)")
                
                cut = min_st_cut(g=self._g, source=source,capacity=self._cap, residual=residual)

                assert cut[source] == True
                assert cut[target] == False

                subtour = [j for j in self._V if cut[self._g.vertex(j)] == False]

                assert len(subtour) < self._instance.n_vertices

                print(f"\tSet size: {len(subtour)}")

                self.__add_sec_for(subtour)

                already_added = already_added.union(subtour)

    def solve(self) -> TspSolution:
        """ Solves the TSP DFJ model and returns the solution.
        
            It throws a RuntimeError if Gurobi cannot solve the model to optimality.

            Returns
            -------
            A TspSolution object with details about the solution.
        """

        # We must set this parameter to use callbacks
        self._model.setParam(GRB.Param.LazyConstraints, 1)
        
        if self.time_limit is not None:
            self._model.setParam('TimeLimit', self.time_limit*60)

        # Gurobi always passes two parameters to the callback: `model` and `where`
        self._model.optimize(lambda _, where: self.__separate(where))

        if self._model.Status != GRB.OPTIMAL:
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