#!/usr/bin/env python3

import os
import base64
import secrets

__version__ = "0.0.3"

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
	token = secrets.token_urlsafe(16)
	return {
		'command': [
			'marimo', 'edit',
			'--port', '{port}',
			'--base-url', os.environ['JUPYTERHUB_SERVICE_PREFIX'] + 'marimo',
			'--token', '--token-password', token,
			'--headless',
			'./'
		],
		'timeout': 60,
		'absolute_url': True,
		'request_headers_override': {
			'Authorization': 'Basic ' + base64.b64encode(b' :' + token.encode()).decode()
		},
		'launcher_entry': {
			'title': 'Marimo',
			'icon_path': os.path.join(os.path.dirname(os.path.abspath(__file__)), 'marimo.svg')
		},
	}
