#!/usr/bin/env python3

from argparse import ArgumentParser
import pandas as pd
from pandas import DataFrame

from typing import Dict, Tuple, List

# Column headers of results file entries 
headers: List[str] = [
    "Implementation", 
    "Input Size", 
    "Clusters", 
    "Dimension", 
    "Nodes", 
    "Tasks Per Node", 
    "Worst Calc. Time", 
    "Worst Exec. Time", 
    "Iterations"
]

# Implementation number as compared to Annes work
names: Dict[str, str] = {
    "own": "Implementation 1",
    "own_inc": "Implementation 2",
    "own_loc": "Implementation 3",
    "own_inc_loc": "Implementation 4"
}

def add_parameters(parser: ArgumentParser) -> None:
    # TODO
    pass

def process(df: DataFrame) -> Dict:
    # TODO
    pass

def process_reverse(df: DataFrame) -> Dict:
    # TODO
    pass

def serialize(data: Dict) -> Tuple[List, Dict]:
    # TODO
    pass

def calc_ci(values, z=1.96) -> Tuple[float, float]:
    # TODO
    pass

def ci_list(values) -> Tuple[List[float], List[float]]:
    # TODO
    pass

def create_confidence_interval(data: Dict, options) -> None:
    # TODO
    pass

def create_boxplots(data: Dict, options) -> None:
    # TODO
    pass
