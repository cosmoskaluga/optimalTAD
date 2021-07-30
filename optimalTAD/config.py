import configparser
import logging
import glob
import os
import sys

log = logging.getLogger(__name__)

def get_configuration():
    config = configparser.ConfigParser()
    conf_list = config.read('config.ini')
    if not conf_list:
        log.error('Cannot opet configuration file!')
        sys.exit(1)
        
    return config

