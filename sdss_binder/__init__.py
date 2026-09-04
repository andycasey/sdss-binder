#!/usr/bin/env python3

import os

__version__ = "0.0.5"

def setup_intro():
	return {
                'command': [],
		'launcher_entry': {
			'title': 'SDSS ❤️  BinderHub',
			'icon_path': os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sdss.svg'),
                        'path_info': 'lab/tree/notebooks/introduction.ipynb'
		},
	}

def setup_marimo():
	return {
		'command': [
			'marimo', 'edit',
			'--port', '{port}',
			'--base-url', os.environ['JUPYTERHUB_SERVICE_PREFIX'] + 'marimo',
			'--no-token',
			'--headless',
			'./'
		],
		'timeout': 60,
		'absolute_url': True,
		'launcher_entry': {
			'title': 'Marimo',
			'icon_path': os.path.join(os.path.dirname(os.path.abspath(__file__)), 'marimo.svg')
		},
	}
