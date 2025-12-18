import os
import typer
import subprocess
import pandas as pd

app = typer.Typer(
	add_completion = True,
	context_settings={'help_option_names': ['-h', '--help']},
	rich_markup_mode = 'rich',
	)

@app.callback(invoke_without_command = True)
def callback(ctx: typer.Context):
	"""
	GUI and CLI for reproducible processing of carbonate IRMS data
	"""
	if ctx.invoked_subcommand is None:
		typer.echo(ctx.get_help())
		raise typer.Exit()

@app.command(name = "gui")
def gui():
	"""
	Launch GUI
	"""
	here = os.path.abspath(os.path.dirname(__file__))
	subprocess.run(['streamlit', 'run',  f'{here}/gui.py'])

@app.command(name = "process")
def cli():
	"""
	Run CLI
	"""
	print("Running CLI now.")
	
def run():
	app()

def process(
	rawdata,
	anchors_df,
):
	try:
		rawdata_df = pd.read_excel(
			rawdata,
			sheet_name = 'Batch Report',
			header = 3,
		).dropna(subset = ['Id'])
		rawdata_df['UID'] = rawdata_df['Id']
		rawdata_df['Sample'] = rawdata_df['Name'].str.split().str[0]
		rawdata_df['d45'] = rawdata_df['45/44 δ⁴⁵CO₂ (raw)']
		rawdata_df['d46'] = rawdata_df['46/44 δ⁴⁶CO₂ (raw)']
		rawdata_df['Time'] = rawdata_df.index + 0.
		rawdata_df = rawdata_df[['UID', 'Time', 'Sample', 'd45', 'd46']]
	except ValueError:
		rawdata_df = pd.read_excel(
			rawdata,
			sheet_name = 'F1CO2-log',
			header = 1,
		)
		rawdata_df['UID'] = rawdata_df['RunIndex']
		rawdata_df['Sample'] = rawdata_df['FileText']
		rawdata_df['d45'] = rawdata_df['RawDelta1']
		rawdata_df['d46'] = rawdata_df['RawDelta2']
		rawdata_df['Time'] = rawdata_df.Index + 0.
		rawdata_df = rawdata_df[['UID', 'Time', 'Sample', 'd45', 'd46']]

	print(rawdata_df)