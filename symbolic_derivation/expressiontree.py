#!/usr/bin/env python3

from sympy import *

from typing import List, Tuple, Dict


class Node:
    def __init__(self, expression: any) -> None:
        self.expression: any = expression
        self.init_components()
        self.init_equations()

    def init_components(self) -> None:
        """ TODO: extract the expression components. """
        self.components: List[any] = list() # components list
        pass

    def init_equations(self) -> None:
        """ TODO: solve the expression for each component and save the results. """
        self.equations: List[Tuple[any, any]] = list() # equations list created by solving for each component
        pass

class Tree:
    def __init__(self, root_nodes: List[Node], conditions: List[any]) -> None:
        self.conditions: List[any] = conditions # candidate conditions
        self.nodes: List[Node] = root_nodes # list of nodes in tree
        self.lineage: Dict[Node, List[Node]] # dict to track lineage of Nodes
        self.tree: Dict[Node, Dict[Tuple[any, any], Node]] = dict() # {src: {(component, substitution): dest} }
    
    def add_node(self, src: Node, equation: Tuple[any, any], dest: Node) -> None:
        """ TODO: add node '(src, (sub_src, sub_dest)) -> dest' to tree. """
        pass

    def search_expansion(self) -> Tuple[any, Tuple[any, any], any]:
        """ TODO: traverse self.nodes and find sets of (src, (sub_src, sub_dest), dest). """
        pass

    def get_candidates(self) -> List[Node]:
        """ TODO: traverse tree and find potential candidates based on set conditions. """
        pass
