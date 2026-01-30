#!/usr/bin/env python3

import os
import setuptools

setuptools.setup(
	name='sdss-binder-instructions',
	version='0.0.1',
	url='https://github.com/andycasey/sdss-binder-instructions',
	author='Andy Casey',
	description='101',
	long_description=open(os.path.join(os.path.dirname(__file__), 'README.md'), 'r', encoding='utf-8').read(),
	long_description_content_type='text/markdown',
	packages=setuptools.find_packages(),
	package_data={
		'sdss-binder-instructions': ['icon.ico'],
	},
	entry_points={
		'jupyter_serverproxy_servers': [
			# name = packagename:function_name
			#'marimo = jupyter_marimo_proxy:setup_marimoserver',
                        'intro = sdss_binder_instructions:setup',
		]
	},
	install_requires=['jupyter-server-proxy'],
)
