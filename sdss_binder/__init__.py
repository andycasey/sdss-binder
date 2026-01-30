#!/usr/bin/env python3

import os

def setup():
	def mappath(path):
		return '/introduction.ipynb'

	return {
		'mappath': mappath,
		'launcher_entry': {
			'title': 'SDSS and BinderHub',
			'icon_path': os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sdss.svg')
		},
	}
