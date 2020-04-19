import sys,os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from maptool.dispatcher.LocalContext import LocalSession
from maptool.dispatcher.LocalContext import LocalContext
from maptool.dispatcher.LazyLocalContext import LazyLocalContext
from maptool.dispatcher.SSHContext import SSHSession
from maptool.dispatcher.SSHContext import SSHContext
from maptool.dispatcher.Dispatcher import _split_tasks

from maptool.dispatcher.LocalContext import _identical_files

def setUpModule():
    os.chdir(os.path.abspath(os.path.dirname(__file__)))
