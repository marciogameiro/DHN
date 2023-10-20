### __init__.py
### MIT LICENSE 2023 Marcio Gameiro

from DHN._dhn import *
from DHN.Graphics import *
from DHN.DrawParameterGraph import *
from DHN.Query.Graph import *
from DHN.Query.Hexcodes import *
import sys
import os
import pickle

configuration().set_path(os.path.dirname(__file__) + '/Resources')
