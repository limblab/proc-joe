#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 23:07:13 2019

@author: minyoungpark
"""

import os
import os.path
import toml


DEFAULT_CONFIG = {
    'filter': {
        'enabled': False,
        'medfilt': 13,
        'offset_threshold': 25,
        'score_threshold': 0.8,
        'spline': True
    },
    'triangulation': {
        'optim': False
    }
}

def full_path(path):
    path_user = os.path.expanduser(path)
    path_full = os.path.abspath(path_user)
    path_norm = os.path.normpath(path_full)
    return path_norm


def load_config(config_filename):
    if config_filename is None:
        config_filename = 'config.toml'

    if os.path.exists(config_filename):
        config = toml.load(config_filename)
    else:
        config = dict()
    
    for k, v in DEFAULT_CONFIG.items():
        if k not in config:
            config[k] = v
        elif isinstance(v, dict): # handle nested defaults
            for k2, v2 in v.items():
                if k2 not in config[k]:
                    config[k][k2] = v2

    return config