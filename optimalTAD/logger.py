import logging
import sys
import os

def initialize_logger():
    filemode = "w"
    logger_name = "optimalTAD"
    format = "%(name)-5s: %(levelname)-10s %(message)s"
    filename = "optad.log"

    dirname = os.path.join(sys.path[0], "log")
    if not os.path.exists(dirname):
        os.makedirs(dirname, exist_ok=True)
    
    logging.basicConfig(level = logging.DEBUG,
                        filename = os.path.join(dirname, filename),
                        filemode = filemode)
        
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter(format)
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)
    logger = logging.getLogger(logger_name)
    return logger
