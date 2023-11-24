#
# Author: Kristian Rietveld, Leiden University
#

from pathlib import Path
import itertools
import re


def create_job_name(params) -> str:
    return "{}_s{}_c{}_d{}_{}x{}".format(params["variant"],
                                         params["size"],
                                         params["clusters"],
                                         params["dimension"],
                                         params["n_nodes"],
                                         params["ntasks_per_node"])

def parse_job_name(jobname : str) -> dict:
    m = re.match(r"([\w_]+)_s(\d+)_c(\d+)_d(\d+)_(\d+)x(\d+)_id.*", jobname)
    if not m:
        raise ValueError("Cannot parse jobname '{}'.".format(jobname))

    labels = ["variant", "size", "clusters", "dimension", "n_nodes",
              "ntasks_per_node"]
    params = { k: v for k, v in zip(labels, m.groups()) }
    return params


def iter_param_combinations(params):
    keys, values = params.keys(), params.values()
    for t in itertools.product(*values):
        current = { key: value for key, value in zip(keys, t) }
        yield current


class DataSetRegistry:
    # Important: this not only indicates the required properties for each
    # dataset, but also specifies the order in which the directories with the
    # property values are nested in the datapath.
    required_properties = ["size", "clusters", "dimension", "seed"]


    def __init__(self, datapath : Path) -> None:
        self.datapath = datapath
        self.datasets = []   # type: list

        self._populate()

    @staticmethod
    def _get_properties_from_path(path : Path) -> dict:
        '''Turn path (e.g. size_28/clusters_4/dimension_4) into a dictionary
        of properties.'''
        properties = filter(lambda s: s != "", path.parts)
        return { key : value for key, value in map(lambda p: p.split("_"), properties) }

    @staticmethod
    def _valid_dataset(properties : dict) -> bool:
        '''A dataset is only valid if it contains all required properties.'''
        for p in DataSetRegistry.required_properties:
            if p not in properties:
                return False
        return True

    @staticmethod
    def _match_query(properties : dict, query : dict) -> bool:
        '''Returns True if all pairs of @query are present in @properties (so
        @properties matches the given @query), False otherwise.'''
        for key, value in query.items():
            if key not in properties or properties[key] not in value:
                return False

        return True

    def _populate(self):
        for datafile in self.datapath.glob("**/data.txt"):
            # remove self.datapath from components
            path = datafile.relative_to(self.datapath)

            properties = self._get_properties_from_path(path.parent)
            if not self._valid_dataset(properties):
                continue

            self.datasets.append((datafile.parent, properties))

    #
    # Public methods
    #

    def query(self, query : dict) -> list:
        # Make sure all specified values in the query are lists of strings
        real_query = dict()
        for k, v in query.items():
            if not isinstance(v, list):
                real_query[k] = [str(v)]
            else:
                real_query[k] = list(map(str, v))

        # Perform the actual query by visiting all known datasets
        result = list()
        for path, properties in self.datasets:
            if self._match_query(properties, real_query):
                result.append((path, properties))

        return result

    def find_datafiles(self, params : dict) -> list:
        # Remove non-required properties, these should not be included in
        # the query.
        query = { k: v for k, v in params.items() if k in DataSetRegistry.required_properties }
        result = self.query(query)
        if not result:
            return []
        return [item[0] for item in result]

    @staticmethod
    def summarize_query(query : dict) -> str:
        tmp = ", ".join([key + "=" + str(val) for key, val in query.items() \
                         if key in DataSetRegistry.required_properties])
        if not tmp:
            tmp = "*"
        return tmp

    @staticmethod
    def get_dir_name(params : dict) -> Path:
        if any(p not in params for p in DataSetRegistry.required_properties):
            return Path()

        components = []
        for p in DataSetRegistry.required_properties:
            components.append(p + "_" + str(params[p]))

        return Path(*components)

    @staticmethod
    def add_dataset_parameters(parser):
        parser.add_argument("--datapath", dest="datapath", type=str,
                            help="Location of the dataset")
        parser.add_argument("--size", dest="size", nargs="*", type=int,
                            help="Size of the dataset")
        parser.add_argument("--clusters", dest="clusters", nargs="*", type=int,
                            help="The number of clusters (k) to find")
        parser.add_argument("--dimension", dest="dimension", nargs="*", type=int,
                            help="The dimension of the data points")
        parser.add_argument("--seed", dest="seed", nargs="*", type=int,
                            help="Random seed the database was generated with")
