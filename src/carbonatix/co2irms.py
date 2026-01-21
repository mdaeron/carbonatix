"""
Process IRMS measurements of δ13C and δ18O in CO2 and carbonates
"""

__author__    = 'Mathieu Daëron'
__contact__   = 'daeron@lsce.ipsl.fr'
__copyright__ = 'Copyright (c) 2025 Mathieu Daëron'
__license__   = 'MIT License - https://opensource.org/licenses/MIT'
__date__      = '2025-12-21'
__version__   = 3.0

import lmfit
import warnings
import numpy as np
import pandas as pd
from types import SimpleNamespace
from scipy.optimize import fsolve
from scipy.stats import t as tstudent


class Config():
	def __init__(
		self,
		R13_VPDB           = 0.01118,    # Chang & Li (1990)
		R18_VSMOW          = 0.0020052,  # Baertschi (1976
		R17_VSMOW          = 0.00038475, # Assonov & Brenninkmeijer (2003) rescaled to match R13_VPDB = 0.01118
		LAMBDA_17          = 0.528,      # Barkan & Luz (2005)
		d18O_VSMOW_of_VPDB = 30.92,      # Coplen et al. (1983), corrected to match d18O_VPDB_of_NBS19 = -2.2 ‰
		D17O_VSMOW_of_VPDB = 0.,         # arbitrary convention
		alpha18_acid       = 1.008129,   #
	):
		self.R13_VPDB = R13_VPDB
		self.R18_VSMOW = R18_VSMOW
		self.R17_VSMOW = R17_VSMOW
		self.LAMBDA_17 = LAMBDA_17
		self.d18O_VSMOW_of_VPDB = d18O_VSMOW_of_VPDB
		self.D17O_VSMOW_of_VPDB = D17O_VSMOW_of_VPDB
		self.alpha18_acid = alpha18_acid

	@property
	def R18_VPDB(self):
		return self.R18_VSMOW * (1 + self.d18O_VSMOW_of_VPDB / 1000)

	@property
	def R17_VPDB(self):
		return (
			self.R18_VSMOW * np.exp(self.D17O_VSMOW_of_VPDB / 1000)
			* (1 + self.d18O_VSMOW_of_VPDB / 1000) ** self.LAMBDA_17
		)

config = Config()

def sanitize(x):
	return x.replace('-', '_').replace('.', '_')

def ratios_to_deltas(
	R45,
	R46,
	D17O = 0, # in permil
):

	R13_VPDB = config.R13_VPDB
	R18_VSMOW = config.R18_VSMOW
	R17_VSMOW = config.R17_VSMOW
	LAMBDA_17 = config.LAMBDA_17

	R45 = np.asarray(R45)
	R46 = np.asarray(R46)
	if R45.shape != R46.shape:
		raise ValueError('R45 and R46 must both be floats or both be arrays of the same shape.')

	def f(R18):
		K = np.exp(D17O/1e3) * R17_VSMOW / R18_VSMOW**LAMBDA_17
		return (-3 * K**2 * R18**(2*LAMBDA_17) + 2 * K * R45 * R18**LAMBDA_17 + 2 * R18 - R46)

	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		R18 = fsolve(f, R46/R18_VSMOW/2, xtol = 1e-16)
	R17 = np.exp(D17O/1e3) * R17_VSMOW * (R18 / R18_VSMOW) ** LAMBDA_17
	R13 = R45 - 2 * R17

	d13C_VPDB = (R13 / R13_VPDB - 1)*1e3
	d18O_VSMOW = (R18 / R18_VSMOW - 1)*1e3

	return (d13C_VPDB, d18O_VSMOW)

def deltas_to_ratios(
	d13C_VPDB,
	d18O_VSMOW,
	D17O = 0, # in permil
):

	R13_VPDB = config.R13_VPDB
	R18_VSMOW = config.R18_VSMOW
	R17_VSMOW = config.R17_VSMOW
	LAMBDA_17 = config.LAMBDA_17

	d13C_VPDB = np.asarray(d13C_VPDB)
	d18O_VSMOW = np.asarray(d18O_VSMOW)
	if d13C_VPDB.shape != d18O_VSMOW.shape:
		raise ValueError('d13C_VPDB and d18O_VSMOW must both be floats or both be arrays of the same shape.')

	R13 = R13_VPDB * (1 + d13C_VPDB/1e3)
	R18 = R18_VSMOW * (1 + d18O_VSMOW/1e3)
	R17 = np.exp(D17O/1e3) * R17_VSMOW * (1 + d18O_VSMOW/1e3)**LAMBDA_17

	R45 = 2 * R17 + R13
	R46 = 2 * R18 + 2 * R17 * R13 + R17**2

	return(R45, R46)

def standardize(
	data,
	anchors,
	constraints = {},
):
	print(f'constraints = {constraints}')
	alpha18_acid = config.alpha18_acid

	data = pd.DataFrame(data)
	if 'UID' not in data:
		data['UID'] = 'A' + (data.index + 1).astype(str)
	data = data.set_index('UID')

	if 'Session' not in data:
		data['Session'] = 'nameless_session'

	fitparams = lmfit.Parameters()

	for s in data.Session.unique():
		fitparams.add('d45_scaling_'+sanitize(s), value = 1.)
		fitparams.add('d46_scaling_'+sanitize(s), value = 1.)
		fitparams.add('d13C_VPDB_of_wg_'+sanitize(s), value = 0.)
		fitparams.add('d18O_VSMOW_of_wg_'+sanitize(s), value = 0.)

	for s in data.Sample.unique():
		fitparams.add('d13C_VPDB_of_'+sanitize(s), value = 0.)
		fitparams.add('d18O_VSMOW_of_'+sanitize(s), value = 0.)

	for a in anchors:
		if a in data.Sample.unique():
			if 'd13C_VPDB' in anchors[a]:
				fitparams[f'd13C_VPDB_of_'+sanitize(a)].expr = str(anchors[a]['d13C_VPDB'])
			if 'd18O_VPDB' in anchors[a]:
				fitparams[f'd18O_VSMOW_of_'+sanitize(a)].expr = str(
					(1000 + anchors[a]['d18O_VPDB']) * alpha18_acid * (1 + config.d18O_VSMOW_of_VPDB / 1000) - 1000
					)
			elif 'd18O_VSMOW' in anchors[a]:
				fitparams[f'd18O_VSMOW_of_'+sanitize(a)].expr = str(anchors[a]['d18O_VSMOW'])

	for p in fitparams:
		if p in constraints:
			fitparams[p].expr = constraints[p]

	def residuals(p, sigma45, sigma46):

		data['d13C_VPDB_of_wg'] = data.Session.map({
			s: p[f'd13C_VPDB_of_wg_{sanitize(s)}'].value
			for s in data.Session.unique()
		})

		data['d18O_VSMOW_of_wg'] = data.Session.map({
			s: p[f'd18O_VSMOW_of_wg_{sanitize(s)}'].value
			for s in data.Session.unique()
		})

		data['d45_scaling'] = data.Session.map({
			s: p[f'd45_scaling_{sanitize(s)}'].value
			for s in data.Session.unique()
		})

		data['d46_scaling'] = data.Session.map({
			s: p[f'd46_scaling_{sanitize(s)}'].value
			for s in data.Session.unique()
		})

		data['d13C_VPDB_true'] = data.Sample.map({
			s: p[f'd13C_VPDB_of_{sanitize(s)}'].value
			for s in data.Sample.unique()
		})

		data['d18O_VSMOW_true'] = data.Sample.map({
			s: p[f'd18O_VSMOW_of_{sanitize(s)}'].value
			for s in data.Sample.unique()
		})

		data['R45wg'], data['R46wg'] = deltas_to_ratios(
			data['d13C_VPDB_of_wg'],
			data['d18O_VSMOW_of_wg'],
		)

		data['R45_true'], data['R46_true'] = deltas_to_ratios(
			data['d13C_VPDB_true'],
			data['d18O_VSMOW_true'],
		)

		data['d45_model'] = data['d45_scaling'] * (data['R45_true'] / data['R45wg'] - 1) * 1000
		data['d46_model'] = data['d46_scaling'] * (data['R46_true'] / data['R46wg'] - 1) * 1000

		return np.hstack((
			(data['d45'] - data['d45_model']) / sigma45,
			(data['d46'] - data['d46_model']) / sigma46,
			))

	N = len(data)
	Nf13, Nf18 = N, N
	for p in fitparams:
		if fitparams[p].expr is not None:
			if p.startswith('d13') or p.startswith('d45'):
				Nf13 -= 1
			if p.startswith('d18') or p.startswith('d46'):
				Nf18 -= 1

	sigma45, sigma46 = 1., 1.

	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		for k in range(2):
			fitresult = lmfit.minimize(
				residuals,
				fitparams,
				method = 'least_squares',
				scale_covar = True,
				args = (sigma45, sigma46),
				x_scale = 'jac',
			)
			sigma45 *= ((fitresult.residual[:N]**2).sum() / Nf13)**.5
			sigma46 *= ((fitresult.residual[-N:]**2).sum() / Nf18)**.5

	p = {_: fitresult.params[_].value for _ in fitresult.params}

	data['d13C_VPDB_of_wg'] = data.Session.map({
		s: p[f'd13C_VPDB_of_wg_{sanitize(s)}']
		for s in data.Session.unique()
	})

	data['d18O_VSMOW_of_wg'] = data.Session.map({
		s: p[f'd18O_VSMOW_of_wg_{sanitize(s)}']
		for s in data.Session.unique()
	})

	data['d45_scaling'] = data.Session.map({
		s: p[f'd45_scaling_{sanitize(s)}']
		for s in data.Session.unique()
	})

	data['d46_scaling'] = data.Session.map({
		s: p[f'd46_scaling_{sanitize(s)}']
		for s in data.Session.unique()
	})

	data['d13C_VPDB_true'] = data.Sample.map({
		s: p[f'd13C_VPDB_of_{sanitize(s)}']
		for s in data.Sample.unique()
	})

	data['d18O_VSMOW_true'] = data.Sample.map({
		s: p[f'd18O_VSMOW_of_{sanitize(s)}']
		for s in data.Sample.unique()
	})
	data['d18O_VPDB_true'] = (1000 + data[f'd18O_VSMOW_true']) / (1 + config.d18O_VSMOW_of_VPDB / 1000) / alpha18_acid - 1000

	data['R45wg'], data['R46wg'] = deltas_to_ratios(
		data['d13C_VPDB_of_wg'],
		data['d18O_VSMOW_of_wg'],
	)

	data['R45_true'], data['R46_true'] = deltas_to_ratios(
		data['d13C_VPDB_true'],
		data['d18O_VSMOW_true'],
	)

	data['R45'] = (1 + data['d45'] / data['d45_scaling'] / 1000) * data['R45wg']
	data['R46'] = (1 + data['d46'] / data['d46_scaling'] / 1000) * data['R46wg']

	data['d13C_VPDB'], data['d18O_VSMOW'] = ratios_to_deltas(
		data['R45'],
		data['R46'],
	)
	data['d18O_VPDB'] = (1000 + data[f'd18O_VSMOW']) / (1 + config.d18O_VSMOW_of_VPDB / 1000) / alpha18_acid - 1000

	data['d13C_VPDB_residual'] = data['d13C_VPDB'] - data['d13C_VPDB_true']
	data['d18O_VSMOW_residual'] = data['d18O_VSMOW'] - data['d18O_VSMOW_true']
	data['d18O_VPDB_residual'] = data['d18O_VPDB'] - data['d18O_VPDB_true']

	data['d13C_anchor'] = data.Sample.map(lambda s: (s in anchors and 'd13C_VPDB' in anchors[s]))
	data['d18O_anchor'] = data.Sample.map(lambda s: (s in anchors and ('d18O_VPDB' in anchors[s] or 'd18O_VSMOW' in anchors[s])))

	out = {'data': data.copy()}

	out['bestfit'] = fitresult
	out['fitreport'] = lmfit.fit_report(fitresult)

	out['sigma_d45'] = sigma45
	out['sigma_d46'] = sigma46
	out['Nf_d45'] = Nf13
	out['Nf_d46'] = Nf18
	out['t95_d45'] = tstudent.ppf(1 - 0.05/2, Nf13)
	out['t95_d46'] = tstudent.ppf(1 - 0.05/2, Nf18)
	out['sigma_d13C_VPDB'] = ((data['d13C_VPDB_residual']**2).sum() / Nf13)**0.5
	out['sigma_d18O_VSMOW'] = ((data['d18O_VSMOW_residual']**2).sum() / Nf18)**0.5
	out['sigma_d18O_VPDB'] = ((data['d18O_VPDB_residual']**2).sum() / Nf18)**0.5
	out['Nf_d13C'] = Nf13
	out['Nf_d18O'] = Nf18
	out['t95_d13C'] = tstudent.ppf(1 - 0.05/2, Nf13)
	out['t95_d18O'] = tstudent.ppf(1 - 0.05/2, Nf18)

	out['sessions'] = {}
	for s in data.Session.unique():
		_s = sanitize(s)
		out['sessions'][s] = {}

		out['sessions'][s]['data'] = out['data'][out['data'].Session == s]
		out['sessions'][s]['N'] = len(out['sessions'][s]['data'])
		out['sessions'][s]['Na_d13C'] = int(data.apply(
			lambda r: r.Session == s and r.Sample in anchors and 'd13C_VPDB' in anchors[r.Sample],
			axis = 1,
		).sum())
		out['sessions'][s]['Na_d18O'] = int(data.apply(
			lambda r: r.Session == s and r.Sample in anchors and ('d18O_VPDB' in anchors[r.Sample] or 'd18O_VSMOW' in anchors[r.Sample]),
			axis = 1,
		).sum())
		out['sessions'][s]['Nu_d13C'] = out['sessions'][s]['N'] - out['sessions'][s]['Na_d13C']
		out['sessions'][s]['Nu_d18O'] = out['sessions'][s]['N'] - out['sessions'][s]['Na_d18O']

		out['sessions'][s]['d45_scaling'] = fitresult.params[f'd45_scaling_'+_s].value
		out['sessions'][s]['d46_scaling'] = fitresult.params[f'd46_scaling_'+_s].value
		out['sessions'][s]['d13C_VPDB_of_wg'] = fitresult.params['d13C_VPDB_of_wg_'+_s].value
		out['sessions'][s]['d18O_VSMOW_of_wg'] = fitresult.params['d18O_VSMOW_of_wg_'+_s].value

		out['sessions'][s]['SE_d45_scaling'] = fitresult.params[f'd45_scaling_'+_s].stderr
		out['sessions'][s]['SE_d46_scaling'] = fitresult.params[f'd46_scaling_'+_s].stderr
		out['sessions'][s]['SE_d13C_VPDB_of_wg'] = fitresult.params['d13C_VPDB_of_wg_'+_s].stderr
		out['sessions'][s]['SE_d18O_VSMOW_of_wg'] = fitresult.params['d18O_VSMOW_of_wg_'+_s].stderr

		out['sessions'][s][f'RMSE_d13C_VPDB'] = (out['sessions'][s]['data']['d13C_VPDB_residual']**2).mean()**0.5
		out['sessions'][s][f'RMSE_d18O_VSMOW'] = (out['sessions'][s]['data']['d18O_VSMOW_residual']**2).mean()**0.5
		out['sessions'][s][f'RMSE_d18O_VPDB'] = (out['sessions'][s]['data']['d18O_VPDB_residual']**2).mean()**0.5

	out['samples'] = {}
	for s in data.Sample.unique():
		_s = sanitize(s)

		out['samples'][s] = {}

		out['samples'][s]['data'] = data[data.Sample == s].copy()
		out['samples'][s]['N'] = len(out['samples'][s]['data'])


		out['samples'][s]['d13C_VPDB'] = fitresult.params['d13C_VPDB_of_'+_s].value
		out['samples'][s]['SD_d13C_VPDB'] = out['samples'][s]['data']['d13C_VPDB'].std()
		if fitresult.params['d13C_VPDB_of_'+_s].stderr:
			out['samples'][s]['is_d13C_anchor'] = False
			out['samples'][s]['SE_d13C_VPDB'] = fitresult.params['d13C_VPDB_of_'+_s].stderr
			out['samples'][s]['95CL_d13C_VPDB'] = out['samples'][s]['SE_d13C_VPDB'] * out['t95_d13C']
		else:
			out['samples'][s]['is_d13C_anchor'] = True

		out['samples'][s]['d18O_VSMOW'] = fitresult.params['d18O_VSMOW_of_'+_s].value
		out['samples'][s]['SD_d18O_VSMOW'] = out['samples'][s]['data']['d13C_VPDB'].std()
		if fitresult.params['d18O_VSMOW_of_'+_s].stderr:
			out['samples'][s]['is_d18O_anchor'] = False
			out['samples'][s]['SE_d18O_VSMOW'] = fitresult.params['d18O_VSMOW_of_'+_s].stderr
			out['samples'][s]['95CL_d18O_VSMOW'] = out['samples'][s]['SE_d18O_VSMOW'] * out['t95_d18O']
		else:
			out['samples'][s]['is_d18O_anchor'] = True

		out['samples'][s]['d18O_VPDB'] = (1000 + out['samples'][s]['d18O_VSMOW']) / (1 + config.d18O_VSMOW_of_VPDB / 1000) / alpha18_acid - 1000
		out['samples'][s]['SD_d18O_VPDB'] = out['samples'][s]['data']['d18O_VPDB'].std()
		if fitresult.params['d18O_VSMOW_of_'+_s].stderr:
			out['samples'][s]['SE_d18O_VPDB'] = out['samples'][s]['SE_d18O_VSMOW'] / (1 + config.d18O_VSMOW_of_VPDB / 1000) / alpha18_acid
			out['samples'][s]['95CL_d18O_VPDB'] = out['samples'][s]['SE_d18O_VPDB'] * out['t95_d18O']

	csv = f'Session,N,Na_d13C,Nu_d13C,Na_d18O,Nu_d18O,d45_scaling,SE_d45_scaling,d46_scaling,SE_d46_scaling,d13C_VPDB_of_wg,SE_d13C_VPDB_of_wg,d18O_VSMOW_of_wg,SE_d18O_VSMOW_of_wg,RMSE_d13C_VPDB,RMSE_d18O_VPDB'
	for s in out['sessions']:
		csv += f'\n{s}'
		csv += f',{out["sessions"][s]["N"]}'
		csv += f',{out["sessions"][s]["Na_d13C"]}'
		csv += f',{out["sessions"][s]["Nu_d13C"]}'
		csv += f',{out["sessions"][s]["Na_d18O"]}'
		csv += f',{out["sessions"][s]["Nu_d18O"]}'
		csv += f',{out["sessions"][s]["d45_scaling"]:.4f}'
		csv += f',{out["sessions"][s]["SE_d45_scaling"]:.4f}'
		csv += f',{out["sessions"][s]["d46_scaling"]:.4f}'
		csv += f',{out["sessions"][s]["SE_d46_scaling"]:.4f}'
		csv += f',{out["sessions"][s]["d13C_VPDB_of_wg"]:.3f}'
		csv += f',{out["sessions"][s]["SE_d13C_VPDB_of_wg"]:.3f}'
		csv += f',{out["sessions"][s]["d18O_VSMOW_of_wg"]:.3f}'
		csv += f',{out["sessions"][s]["SE_d18O_VSMOW_of_wg"]:.3f}'
		csv += f',{out["sessions"][s]["RMSE_d13C_VPDB"]:.3f}'
		csv += f',{out["sessions"][s]["RMSE_d18O_VPDB"]:.3f}'
	out['csv_sessions'] = csv

	csv = f'Sample,N,d13C_VPDB,SE_d13C_VPDB,95CL_d13C_VPDB,SD_d13C_VPDB,d18O_VPDB,SE_d18O_VPDB,95CL_d18O_VPDB,SD_d18O_VPDB'
	for s in out['samples']:
		csv += f'\n{s}'
		csv += f',{out["samples"][s]["N"]}'
		csv += f',{out["samples"][s]["d13C_VPDB"]:.3f}'
		try:
			csv += f',{out["samples"][s]["SE_d13C_VPDB"]:.3f}'
			csv += f',{out["samples"][s]["95CL_d13C_VPDB"]:.3f}'
		except KeyError:
			csv += ',,'
		csv += f',{out["samples"][s]["SD_d13C_VPDB"]:.3f}'
		csv += f',{out["samples"][s]["d18O_VPDB"]:.3f}'
		try:
			csv += f',{out["samples"][s]["SE_d18O_VPDB"]:.3f}'
			csv += f',{out["samples"][s]["95CL_d18O_VPDB"]:.3f}'
		except KeyError:
			csv += ',,'
		csv += f',{out["samples"][s]["SD_d18O_VPDB"]:.3f}'
	out['csv_samples'] = csv

	out['csv_analyses'] = data[[
		'Session',
		'Sample',
		'd45',
		'd46',
		'd13C_VPDB',
		'd18O_VSMOW',
		'd18O_VPDB',
		'd13C_VPDB_residual',
		'd18O_VSMOW_residual',
		'd18O_VPDB_residual',
	]].to_csv(float_format = '%.4f')

	return out

def plot_residuals(
	sdata,
	field,
	fig = None,
	figsize = (5,3),
	ax = None,
	s = 6,
	marker = 'x',
	alpha = 1,
	linewidths = 1,
	color_unknowns = 'k',
	color_anchors = 'r',
	show_sample_label = True,
	show_uid_label = True,
	sample_label_kw = dict(size = 6, alpha = 0.3),
	uid_label_kw = dict(size = 6, alpha = 0.2, color = 'k'),
):
	from matplotlib import pyplot as ppl
	out = SimpleNamespace()

	if fig is None and ax is None:
		out.fig = ppl.figure(figsize = figsize)
		out.ax = out.fig.add_subplot(111)
	elif ax is None:
		out.fig = fig
		out.ax = out.fig.add_subplot(111)
	elif fig is None:
		out.fig = ax.get_figure()
		out.ax = ax
	else:
		if fig == ax.get_figure():
			out.fig = fig
			out.ax = ax
		else:
			raise ValueError("Both `fig` and `ax` were specified, but `fig` is not the parent figure of `ax`.")

	N = len(sdata['data'])
	Ns = len(sdata['sessions'])
	session_offsets = {s: k for k,s in enumerate(sdata['sessions'])}
	X = sdata['data']['Session'].map(session_offsets) + range(N) + 1
	Y = sdata['data'][f'{field}_residual']
	anchorfield = f'{field.split('_')[0]}_anchor'
	colors = sdata['data'][anchorfield].map(
		{
			True: color_anchors,
			False: color_unknowns,
		}
	)
	out.ax.scatter(X, Y, s = s, c = colors, marker = marker, alpha = alpha, linewidths = linewidths)
	if show_sample_label:
		for x,y,s,c in zip(X, Y, sdata['data']['Sample'], colors):
			out.ax.text(
				x, y, s + '\n',
				**(dict(color = c, va = 'center', ha = 'center') | sample_label_kw),
			)
	if show_uid_label:
		for x,y,u,c in zip(X, Y, sdata['data'].index, colors):
			out.ax.text(
				x, y, '\n' + u,
				**(dict(color = c, va = 'center', ha = 'center') | uid_label_kw),
			)


	return out
