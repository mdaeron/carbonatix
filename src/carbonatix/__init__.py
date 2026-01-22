import os
import typer
import subprocess
import pandas as pd
from . import co2irms

from .__version__ import __version__ as __version__

__author__ = 'Mathieu Daëron'
__contact__ = 'mathieu@daeron.fr'

app = typer.Typer(
	add_completion = True,
	context_settings={'help_option_names': ['-h', '--help']},
	rich_markup_mode = 'rich',
	)

@app.callback(invoke_without_command = True)
def main():
	"""
	GUI for reproducible processing of carbonate IRMS data
	"""
	here = os.path.abspath(os.path.dirname(__file__))
	subprocess.run(['streamlit', 'run',  f'{here}/gui.py'])

def run():
	app()

def process(
	rawdata_df,
	settings,
):

	data = rawdata_df.to_dict('records')
	anchors = {}
	for a in settings['anchors']:
		anchors[a] = settings['anchors'][a].copy()
		if 'aka' in anchors[a] and anchors[a]['aka']:
			for k in anchors[a]['aka']:
				anchors[k] = anchors[a]
			anchors[a].pop('aka')
	isoparams = {
		k: settings[k]
		for k in [
			'R13_VPDB',
			'R18_VSMOW',
			'R17_VSMOW',
			'LAMBDA_17',
			'd18O_VSMOW_of_VPDB',
			'D17O_VSMOW_of_VPDB',
			'alpha18_acid',
		]
	}
	constraints = settings['constraints']

	co2irms.config = co2irms.Config(**isoparams)

	S = co2irms.standardize(
		data = data,
		anchors = anchors,
		constraints = constraints,
	)

	return S
