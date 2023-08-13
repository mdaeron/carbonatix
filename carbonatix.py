import co2irms, zipfile
import streamlit as st
import numpy as np
import pandas as pd
from io import StringIO, BytesIO
from datetime import datetime as dt
from matplotlib.pyplot import style

__version__ = 1.0

st.set_page_config(
	page_title = 'Carbonatix',
	layout = 'wide',
	)
st.title('Carbonatix')

st.markdown('### Raw Data')
data = st.data_editor(
	pd.DataFrame([
		dict(Session = 'mysession', UID = uid, Sample = s, d45 = d45, d46 = d46)
		for uid,s,d45,d46 in [
			('30445', '603-1', 13.446903, 8.0692674),
			('30446', '603-1', 13.416918, 7.9028158),
			('30447', '603-1', 13.377629, 7.8887903),
			('30448', '603-1', 13.318248, 7.7940762),
			('30449', '603-1', 13.338062, 7.8398809),
			('30450', 'NBS-19', 12.874512, 8.0758038),
			('30451', 'NBS-19', 12.864055, 8.0118926),
			('30452', 'NBS-19', 12.721391, 7.8640607),
			('30453', 'NBS-19', 12.757975, 7.7891145),
			('30454', 'NBS-19', 12.88145, 8.0796388),
			('30455', 'WILEY', 10.493837, 3.0597417),
			('30456', 'WILEY', 10.404619, 2.9493359),
			('30457', 'WILEY', 10.476399, 3.0388496),
			('30458', 'WILEY', 10.382638, 2.9341232),
			('30459', 'WILEY', 10.443481, 2.8872835),
			('30460', 'NBS-18', 5.5641596, -12.6756),
			('30461', 'NBS-18', 5.6087321, -12.64079),
			('30462', 'NBS-18', 5.5814457, -12.78631),
			('30463', 'NBS-18', 5.5190723, -12.66599),
			('30464', 'NBS-18', 5.6119561, -12.62632),
			('30465', 'CAMORA', 13.468822, 7.6871945),
			('30466', 'CAMORA', 13.505947, 7.732545),
			('30467', 'CAMORA', 13.468423, 7.706157),
			('30468', 'CAMORA', 13.316551, 7.4070193),
			('30469', 'CAMORA', 13.341697, 7.4490296),
			('30470', '603-1', 13.098373, 7.5428612),
			('30471', '603-1', 13.245609, 7.6986331),
			('30472', '603-1', 12.894225, 7.222693),
			('30473', '603-1', 13.142402, 7.4075872),
			('30474', '603-1', 13.223274, 7.5159031),
			]
		]),
	num_rows = 'dynamic',
# 	use_container_width = True,
	)

st.markdown('### Anchors')
_anchors = st.data_editor(
	pd.DataFrame([
		dict(Sample = 'NBS-19', d13C_VPDB = 1.95, d18O_VPDB = -2.2),
		dict(Sample = 'NBS-18', d13C_VPDB = -5.014, d18O_VPDB = -23.01),
		]),
	num_rows = 'dynamic',
	)

if st.button('Standardize'):
	
	R45_VPDB, R46_VPDB = co2irms.deltas_to_ratios(0, 30.92)
	
	anchors = {a['Sample']: a for a in _anchors.to_dict('records')}
	
	data = data.to_dict('records')
	for r in data:
		if 'Session' not in r:
			r['Session'] = 'mysession'

	S = co2irms.standardize(data, anchors)

	st.markdown('### Results')

	st.markdown(f'''
- **Repeatability of δ13C:** {S["sigma_d13C"]:.3f} ‰ (1SD)
- **Repeatability of δ18O:** {S["sigma_d18O"]:.3f} ‰ (1SD)
'''[1:-1])
	
	st.markdown('#### Session parameters')
	st.data_editor(pd.read_csv(StringIO(S['csv_sessions'])), hide_index = True)
	
	st.markdown('#### Samples')
	st.data_editor(pd.read_csv(StringIO(S['csv_samples'])), hide_index = True)

	st.markdown('#### Analyses')
	df = pd.read_csv(StringIO(S['csv_analyses']))
	df['UID'] = df['UID'].apply(str)
	df['Sample'] = df['Sample'].apply(str)
	df['Session'] = df['Session'].apply(str)
	st.data_editor(df, hide_index = True)

	style.use('mydefault.mplstyle')
	fig = co2irms.plot_residuals(S)
	st.pyplot(fig, use_container_width = False, dpi = 120)

	buf = BytesIO()
	plotbuf = BytesIO()

	readme = f'''
Carbonatix version {__version__}, using co2irms version {co2irms.__version__}
Processed on {dt.now().strftime('%Y-%m-%d %H:%M:%S')}

Contents:

* analyses.csv  : table of analyses
* anchors.csv   : table of the anchors used to standardize δ13C, δ18O, Δ47, and Δ48 measurements
* residuals.pdf : plot of the δ13C and δ18O residuals for each analysis
* samples.csv   : table of the final values for each sample
* sessions.csv  : table of the standardization parameters for each session


Carbonatix [1] is an open-source web app based on the co2irms library [2],
using standardization and error propagation methods inspired by
Daëron (2021) [https://doi.org/10.1029/2020GC009592].

Please open an issue [3] for bug reports, questions and/or suggestions.

[1] https://github.com/mdaeron/carbonatix
[2] https://github.com/mdaeron/co2irms
[3] https://github.com/mdaeron/carbonatix/issues
'''[1:]

	with zipfile.ZipFile(buf, 'x') as dl_zip:
		dl_zip.writestr('analyses.csv', S['csv_analyses'])
		dl_zip.writestr('sessions.csv', S['csv_sessions'])
		dl_zip.writestr('samples.csv', S['csv_samples'])
		dl_zip.writestr('anchors.csv', _anchors.to_csv(index = False))
		dl_zip.writestr('readme.txt', readme)
		fig.savefig(plotbuf, format='pdf')
		plotbuf.seek(0)
		dl_zip.writestr('residuals.pdf', plotbuf.read())

	st.download_button(
		label = 'Download zip',
		data = buf.getvalue(),
		file_name = 'carbonatix-results.zip',
		mime = 'application/zip',
		)

st.write(f'''
<div style="color:#999999; margin-top:4ex;">
<a style="color:inherit" href="https://github.com/mdaeron/carbonatix">Carbonatix</a>
is an open-source web app based on the
<a style="color:inherit" href="https://github.com/mdaeron/co2irms">co2irms</a>
library, using standardization and error propagation methods inspired by
<a style="color:inherit" href="https://doi.org/10.1029/2020GC009592">Daëron (2021)</a>.

<p>
Please <a style="color:inherit" href="https://github.com/mdaeron/carbonatix/issues">open an issue</a>
for bug reports, questions and/or suggestions.
</div>
''', unsafe_allow_html = True)
