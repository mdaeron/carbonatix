import os
import io
import pathlib
import zipfile
import tomlkit
import pandas as pd
import datetime as dt
import streamlit as st
import carbonatix as cx

st.set_option("client.toolbarMode", "viewer")

st.set_page_config(
	page_title = f'Carbonatix (v{cx.__version__})',
	layout = 'wide',
)

st.write('''
	<style>
		h1 {
			margin-top: -4ex !important;
			}
	</style>
''', unsafe_allow_html = True)

with st.container(horizontal = True, horizontal_alignment = 'right'):
	st.markdown('# Carbonatix')
	downloadzip_button = st.empty()
	standardize_button = st.button("Standardize")

if standardize_button or 'S' in st.session_state:
	tab_rawdata, tab_settings, tab_constraints, tab_results, tab_plots = st.tabs(["Input Data", "Settings", 'Constraints', "Results", "Plots"])
else:
	tab_rawdata, tab_settings, tab_constraints = st.tabs(["Input Data", "Settings", 'Constraints'])

default_settings = dict(
	R13_VPDB           = 0.01118,     # Chang & Li (1990)
	R18_VSMOW          = 0.0020052,   # Baertschi (1976)
	R17_VSMOW          = 0.00038475,  # Assonov & Brenninkmeijer (2003) rescaled to match R13_VPDB = 0.01118
	LAMBDA_17          = 0.528,       # Barkan & Luz (2005)
	d18O_VSMOW_of_VPDB = 30.92,       # Coplen et al. (1983), corrected to match d18O_VPDB_of_NBS19 = -2.2 ‰
	D17O_VSMOW_of_VPDB = 0.0,         # arbitrary convention
	alpha18_acid       = 1.008129,    # Kim et al. (2007) for calcite at 90 °C

	anchors = {
		'IAEA603': dict(
			d13C_VPDB = 2.46,
			d18O_VPDB = -2.37,
			aka = ['IAEA-603', 'IAEA_603', 'IAEA 603'],
		),
		'NBS18': dict(
			d13C_VPDB = -5.01,
			d18O_VPDB = -23.01,
			aka = ['NBS-18', 'NBS_18', 'NBS 18'],
		),
	},
)

try:
	with open('carbonatix-defaults.toml') as fid:
		settings = tomlkit.load(fid).unwrap()
except FileNotFoundError:
	settings = default_settings

with tab_settings:

	col_left, col_right = st.columns(2)

	with col_left:
		if 'anchors' in settings:

			_anchors_list = []

			for a in settings['anchors']:
				_anchors_list.append(
					{'Sample': a} | {
						k: v for k,v in settings['anchors'][a].items()
						if k in ['d13C_VPDB', 'd18O_VPDB', 'd18O_VSMOW', 'aka']
					}
				)
		else:
			_anchors_list = [
				dict(Sample = 'IAEA603', d13C_VPDB =  '2.46', d18O_VPDB =  '-2.37', d18O_VSMOW = None, aka = None),
				dict(Sample = 'NBS18',    d13C_VPDB = '-5.01', d18O_VPDB = '-23.01', d18O_VSMOW = None, aka = None),
			]

		_anchors_df = pd.DataFrame(_anchors_list)

		st.markdown('### Reference Materials')

		anchors_df = st.data_editor(
			_anchors_df,
			width = 'content',
			# hide_index = True,
			num_rows = 'dynamic',
			column_config = {
				'd13C_VPDB': st.column_config.NumberColumn('δ13C_VPDB', format = '%+.3f'),
				'd18O_VPDB': st.column_config.NumberColumn('δ18O_VPDB', format = '%+.3f'),
				'd18O_VSMOW': st.column_config.NumberColumn('δ18O_VSMOW', format = '%+.3f'),
			},
		)

		anchors_df.d13C_VPDB = anchors_df.d13C_VPDB.astype(float)
		anchors_df.d18O_VPDB = anchors_df.d18O_VPDB.astype(float)

		settings['anchors'] = anchors_df.set_index('Sample').to_dict('index')

		with st.expander('Help on reference materials'):
			st.markdown("""
* Reference Materials (RMs) for isotopic standardization have nominal values in δ<sup>13</sup>C and/or δ<sup>18</sup>O.
* Carbonate RMs must have nominal values in δ<sup>13</sup>C<sub>VPDB</sub> and/or δ<sup>18</sup>O<sub>VPDB</sub>.
* CO<sub>2</sub> RMs must have nominal values in δ<sup>13</sup>C<sub>VPDB</sub> and/or δ<sup>18</sup>O<sub>VSMOW</sub>.
* No RM may have nominal values in both δ<sup>18</sup>O<sub>VPDB</sub> and/or δ<sup>18</sup>O<sub>VSMOW</sub>.
* The <u>Sample</u> entry (or one of the aliases defined in <u>aka</u>) must match exactly the <u>Sample</u> entry in your raw data (next tab)
""", unsafe_allow_html = True)

	with col_right:

		st.markdown('### Isotopic parameters')
		isoparams = {}
		for k, v in (
			('R13_VPDB',  0.01118),          # Chang & Li (1990)
			('R18_VSMOW', 0.0020052),        # Baertschi (1976)
			('R17_VSMOW', 0.00038475),       # Assonov & Brenninkmeijer (2003) rescaled to match R13_VPDB = 0.01118
			('LAMBDA_17', 0.528),            # Barkan & Luz (2005)
			('d18O_VSMOW_of_VPDB', 30.92),   # Coplen et al. (1983), corrected to match d18O_VPDB_of_NBS19 = -2.2 ‰
			('D17O_VSMOW_of_VPDB', 0.0),     # arbitrary convention
			('alpha18_acid', 1.008129),      # Kim et al. (2007) for calcite at 90 °C
		):
			if k in settings:
				isoparams[k] = settings[k]
			else:
				isoparams[k] = v

		updated_isoparams = st.data_editor(
			isoparams,
			width = 'content',
		)
		for k in isoparams:
			settings[k] = updated_isoparams[k]

		st.markdown('### Human operator')
		# with st.container(horizontal = True, horizontal_alignment = 'left'):

		settings['operator']['name'] = st.text_input(
			label = 'Name',
			value = settings['operator']['name']
				if 'operator' in settings and 'name' in settings['operator']
				else "",
			)
		settings['operator']['email'] = st.text_input(
			label = 'Email',
			value = settings['operator']['email']
				if 'operator' in settings and 'email' in settings['operator']
				else "",
			)

with tab_rawdata:
	st.markdown('### Input Data')
	uploaded_files = sorted(
		st.file_uploader(label = ' ', accept_multiple_files = True),
		key = lambda _: _.name,
	)

	error_placeholder = st.empty()

	if uploaded_files:
		rawdata_df = pd.DataFrame()

		for session_autonumber, uploaded_file in enumerate(uploaded_files):

			if uploaded_file.name.endswith('.xlsx'):
				_rawdata_df = pd.read_excel(
					uploaded_file,
					sheet_name = 'Batch Report',
					header = 3,
				).dropna(subset = ['Id'])
				if 'Session' not in _rawdata_df:
					_rawdata_df['Session'] = uploaded_file.name[:-5].split('_')[-1]
				_rawdata_df['UID'] = _rawdata_df['Id']
				_rawdata_df['Sample'] = _rawdata_df['Name'].str.split().str[0]
				_rawdata_df['d45'] = _rawdata_df['45/44 δ⁴⁵CO₂ (raw)']
				_rawdata_df['d46'] = _rawdata_df['46/44 δ⁴⁶CO₂ (raw)']
				_rawdata_df['Time'] = _rawdata_df.index + 0.
				_rawdata_df = _rawdata_df[['UID', 'Session', 'Time', 'Sample', 'd45', 'd46']].set_index('UID')
			elif uploaded_file.name.endswith('.xls'):
				_rawdata_df = pd.read_excel(
					uploaded_file,
					sheet_name = 'F1CO2-log',
					header = 1,
				)
				if 'Session' not in _rawdata_df:
					_rawdata_df['Session'] = 'nameless_session'
				_rawdata_df['UID'] = _rawdata_df['RunIndex']
				_rawdata_df['Sample'] = _rawdata_df['FileText']
				_rawdata_df['d45'] = _rawdata_df['RawDelta1']
				_rawdata_df['d46'] = _rawdata_df['RawDelta2']
				_rawdata_df['Time'] = _rawdata_df.Index + 0.
				_rawdata_df = _rawdata_df[['UID', 'Session', 'Time', 'Sample', 'd45', 'd46']].set_index('UID')

			elif uploaded_file.name.endswith('.csv'):
				_rawdata_df = pd.read_csv(
					uploaded_file,
				)
				if 'Session' not in _rawdata_df:
					_rawdata_df['Session'] = 'nameless_session'
				_rawdata_df = _rawdata_df[['UID', 'Session', 'Time', 'Sample', 'd45', 'd46']].set_index('UID')

			_rawdata_df = st.data_editor(
				_rawdata_df,
				height = 'content',
				# hide_index = True,
				num_rows = 'dynamic',
				column_config = {
					'Sample': st.column_config.TextColumn('Sample'),
					'd45': st.column_config.NumberColumn('d45', format = '%+.3f'),
					'd46': st.column_config.NumberColumn('d46', format = '%+.3f'),
				},
			)

			rawdata_df = pd.concat([rawdata_df, _rawdata_df], ignore_index = True)

with tab_constraints:

	if uploaded_files:

		st.markdown('### Additional Constraints')
		constraints = {}
		for session in sorted(rawdata_df.Session.unique()):
			constraints[f'd13C_VPDB_of_wg_{cx.co2irms.sanitize(session)}'] = None
			constraints[f'd18O_VSMOW_of_wg_{cx.co2irms.sanitize(session)}'] = None
			constraints[f'd45_scaling_{cx.co2irms.sanitize(session)}'] = None
			constraints[f'd46_scaling_{cx.co2irms.sanitize(session)}'] = None

		for sample in sorted(rawdata_df.Sample.unique()):
			constraints[f'd13C_VPDB_of_{cx.co2irms.sanitize(sample)}'] = None
			constraints[f'd18O_VSMOW_of_{cx.co2irms.sanitize(sample)}'] = None

		if 'constraints' in settings:
			for k in settings['constraints']:
				constraints[k] = settings['constraints'][k]

		updated_constraints = st.data_editor(
			constraints,
			height = 'content',
		)

		settings['constraints'] = {
			k: str(v)
			for k,v in updated_constraints.items()
			if v is not None
		}

if standardize_button:

	if uploaded_files:
		S = cx.process(
			rawdata_df = rawdata_df,
			settings = settings,
		)
		st.session_state['S'] = S
	else:
		error_placeholder.error('Missing raw data!')

if 'S' in st.session_state:
	S = st.session_state['S']

	with tab_results:
		st.markdown('### Results')


		st.markdown('#### Samples')
		st.data_editor(
			pd.read_csv(io.StringIO(S['csv_samples'])).set_index('Sample'),
			width = 'content',
			height = 'content',
		)

		st.markdown('#### Sessions')
		st.data_editor(
			pd.read_csv(io.StringIO(S['csv_sessions'])).set_index('Session').transpose(),
			width = 'content',
			height = 'content',
		)

		st.markdown('#### Analyses')
		st.data_editor(
			pd.read_csv(io.StringIO(S['csv_analyses'])).set_index('UID'),
			width = 'content',
			height = 'content',
		)

	with tab_plots:
		st.markdown('### Plots')
		plots = {}

		out = cx.co2irms.plot_residuals(S, 'd13C_VPDB', title = '$δ^{13}C_{VPDB}$ residuals (‰)')
		st.pyplot(out.fig, width = 600)
		plots['residuals_d13C'] = out

		out = cx.co2irms.plot_residuals(S, 'd18O_VPDB', title = '$δ^{18}O_{VPDB}$ residuals (‰)')
		st.pyplot(out.fig, width = 600)
		plots['residuals_d18O'] = out

		out = cx.co2irms.plot_correction(S, 'd45', 'd13C_VPDB')
		st.pyplot(out.fig, width = 400)
		plots['correction_d13C'] = out

		out = cx.co2irms.plot_correction(S, 'd46', 'd18O_VPDB')
		st.pyplot(out.fig, width = 400)
		plots['correction_d18O'] = out

		readme = f"""
This zipfile was generated by carbonatix on {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.

To reprocess the data, run the following command from the root of this directory:

uv run carbonatix -i 'input/??? -a 'input/anchors.csv'
"""

		here = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))
		buf = io.BytesIO()
		with zipfile.ZipFile(buf, 'x') as dl_zip:
			dl_zip.writestr('readme.txt', readme[1:])
			with open(here / 'assets/pyproject.toml') as fid:
				dl_zip.writestr('pyproject.toml', fid.read())
			with open(here / 'assets/uv.lock') as fid:
				dl_zip.writestr('uv.lock', fid.read())
			for uploaded_file in uploaded_files:
				uploaded_file.seek(0)
				dl_zip.writestr(f'input/{uploaded_file.name}', uploaded_file.read())
			with io.BytesIO() as fid:
				rawdata_df.to_csv(fid, index = False)
				fid.seek(0)
				dl_zip.writestr(f'input/rawdata.csv', fid.read())
			with io.StringIO() as fid:
				tomlkit.dump(settings, fid)
				fid.seek(0)
				dl_zip.writestr(f'input/carbonatix-defaults.toml', fid.read())
			for plotname in plots:
				pdf_buffer = io.BytesIO()
				plots[plotname].fig.savefig(pdf_buffer, format = "pdf")
				pdf_buffer.seek(0)
				dl_zip.writestr(f'results/plots/{plotname}.pdf', pdf_buffer.read())

			dl_zip.writestr(f'results/samples.csv', S['csv_samples'])
			dl_zip.writestr(f'results/sessions.csv', S['csv_sessions'])
			dl_zip.writestr(f'results/analyses.csv', S['csv_analyses'])

		downloadzip_button.download_button(
			label = 'Download zip',
			data = buf.getvalue(),
			file_name = 'carbonatix-results.zip',
			mime = 'application/zip',
			)
