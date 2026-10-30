import pandas as pd
import numpy as np 
import re

axial_tm = pd.read_excel('descriptor.xlsx', sheet_name='axial', index_col=0)
lateral_tm = pd.read_excel('descriptor.xlsx', sheet_name='lateral', index_col=0)
backbone = pd.read_excel('descriptor.xlsx', sheet_name='propensity', index_col=0)

class PredictXY:
	def __init__(self, seq, m):
		if isinstance(seq, str):
			self.seq = {'a':seq, 'b':seq, 'c':seq}
		elif isinstance(seq, list):
			self.seq = {'a':seq[0], 'b':seq[1], 'c':seq[2]}
		else:
			print("please input the right seq formation")
		self.m = str(m)

	def pro(self):
		m = self.m

		config_map = {
			"1": ('GPY', 'GPY', 'GPY', 'GPY', 'GXO', 'GXO', 'GXO', 'GXO', 'GXO', 'GXO'),
			"2": ('GPY-KD', 'GPY', 'GPY-KD', 'GPY', 'GXO-HB', 'GXO', 'GXO-HB', 'GXO', 'GXO', 'GXO-HB'),
			"3": ('GPY', 'GPY', 'GY', 'GY', 'GXO', 'GXO', 'GXP', 'GXP', 'GX', 'GX'),
			"4": ('GPY-KD', 'GPY', 'GY-KD', 'GY', 'GXO-HB', 'GXO', 'GXP', 'GXP', 'GX', 'GX')
		}
		else_keys = ('GPY-KD', 'GPY', 'GY-KD', 'GY', 'GXO-HB', 'GXO', 'GXP', 'GXP', 'GX', 'GX')

		keys = config_map.get(m, else_keys)

		return tuple(backbone[k] for k in keys)
		#GPY_KD, GPY_no, GY_KD, GY_no, GXO_KD, GXO_no, GXP_no, GX_no, GX_KD


	def GXY(self):
		m = self.m
		seq = self.seq
		xy = dict()
		for chain in ['a', 'b', 'c']:
			xyg = list(seq[chain])
			if seq['a'][0] == 'G':
				del (xyg[2::3])
			else:
				del (xyg[2::3])
			xy[chain] = xyg

		alltris = {}
		matriXY = np.transpose(np.array([['Z','Z'] + xy['c'][:-2],
										xy['a'], xy['b'], xy['c'], 
										xy['a'][2:]+['Z','Z']]))

		GPY_KD, GPY_no, GY_KD, GY_no, GXO_KD, GXO_no, GXP_KD, GXP_no, GX_no, GX_KD = self.pro()
		proXY = {'GXO_KD':GXO_KD, 'GXO_no':GXO_no,
		 	     'GXP_KD':GXP_KD, 'GXP_no':GXP_no}

		for row in range(0, len(xy['a'])//2):
			for column in range(1, 4):
				triplets = matriXY[row*2:row*2+2, column]
				if triplets[0] == 'P': # X=P
					if triplets[1] in ['O', 'P']: # GPO
						name = triplets[0]+triplets[1]
						alltris.setdefault('G'+name, []).append('%s'%name)
					else: # GPY
						try:
							if matriXY[row*2+1, column-1]+matriXY[row*2+2, column] in ['KD', 'KE']:  
								alltris.setdefault('GPY_KD', []).append(GPY_KD[triplets[1]])
							else: 
								alltris.setdefault('GPY_no', []).append(GPY_no[triplets[1]])
						except: alltris.setdefault('GPY_no', []).append(GPY_no[triplets[1]])
				else:  # X != P
					if triplets[1] == 'O':  #GXO
						key_str = 'GX%s_KD'%triplets[1] if matriXY[row*2, column+1] != 'P' else 'GX%s_KD'%triplets[1]
						alltris.setdefault(key_str, []).append(proXY[key_str][triplets[0]])
					elif triplets[1] == 'P':
						alltris.setdefault('GX%s_no'%triplets[1], []).append(proXY['GX%s_no'%triplets[1]][triplets[0]])

					else:#GXY_GX
						try:
							key_str = 'GX_KD' if matriXY[row*2, column+1] != 'P' else 'GX_no'
							alltris.setdefault(key_str, []).append(GX_KD[triplets[0]])
						except: alltris.setdefault('GX_no', []).append(GX_no[triplets[0]])
						# GXY_GY
						try:
							gy_key = ('GY_KD' if matriXY[row*2+1, column-1]+matriXY[row*2+2, column] in ['KD', 'KE'] else 'GY_no')
							alltris.setdefault(gy_key, []).append(GY_KD[triplets[1]] if gy_key=='GY_KD' else GY_no[triplets[1]])
						except: alltris.setdefault('GY_no', []).append(GY_no[triplets[1]])

		newgxy = dict.fromkeys(['GPO', 'GPP', 
								'GPY_no', 'GPY_KD', 
								'GXO_no', 'GXO_KD', 
								'GXP_no', 'GXP_KD',
								'GX_no', 'GY_no', 
								'GX_KD', 'GY_KD'], 0)
		trangle_num = dict.fromkeys(['GPO', 'GPP', 
								'GPY_no', 'GPY_KD', 
								'GXO_no', 'GXO_KD', 
								'GXP_no', 'GXP_KD',
								'GX_no', 'GY_no', 
								'GX_KD', 'GY_KD'], 0)
		for k, v in alltris.items():
			if isinstance(v[0], str):
				if k == 'GPP':
					newgxy[k] = v.count(k[1:])*-0.6
			else:
				newgxy[k] = sum(v)

		return newgxy, matriXY


	def Side(self):
		seq = self.seq
		yaa, xaa = axial_tm.columns.values.tolist(), axial_tm.index.tolist()
		
		pairwise = {'axial': {'01':{}, '12': {}, '23': {}},
					'lateral': {'01':{}, '12':{}, '23':{}}}
		matriXY = self.GXY()[1]

		for col in range(3):
			for row in range(len(matriXY)//2):
				aa_y = matriXY[row*2+1, col]
				aa_x2 = matriXY[row*2, col+1]
				if (aa_x2 in xaa) and (aa_y in yaa):
					pairwise['lateral'][str(col)+str(col+1)][str(row+1)] = lateral_tm[aa_y][aa_x2]
				else: pairwise['lateral'][str(col)+str(col+1)][str(row+1)] = 0
				try:
					aa_x = matriXY[(row+1)*2, col+1]
					if (aa_x in xaa) and (aa_y in yaa):
						pairwise['axial'][str(col)+str(col+1)][str(row+1)] = axial_tm[aa_y][aa_x]
					else: pairwise['axial'][str(col)+str(col+1)][str(row+1)] = 0
				except IndexError: continue
		return pairwise

	def wSide(self, *args):
		pairwise = self.Side()

		def calc_sum(data_dict):
			d_list = [list(data_dict[k].values()) for k in ['01', '12', '23']]	
			
			return sum(x * w for x, w in zip(d_list[0], args))  + sum(sum(other) for other in d_list[1:-1]) + sum(x * args[-1] for x in d_list[-1])	
		
		return calc_sum(pairwise['axial']), calc_sum(pairwise['lateral'])
		
	def weights(self):
		weight_map = {
			'1':[0.26, 0.26, 0.26, 0.64, 0.64, 0.64, 0.64, 0.64, 0.26, 0.64, 0.26, 0.26, 0.01, 0.82],
			'2':[0.17, 0.17, 0.33, 0.79, 0.86, 0.79, 0.86, 0.79, 0.17, 0.86, 0.33, 0.26, 0.01, 0.82],
			'3':[0.39, 0.39, 0.39, 0.38, 0.38, 0.60, 0.60, 0.50, 0.01, 0.50, 0.01, 0.26, 0.01, 0.82],
			'4':[0.35, 0.35, 0.36, 0.38, 0.86, 0.71, 0.71, 0.39, 0.16, 0.39, 0.03, 0.26, 0.01, 0.82],}
		return weight_map.get(self.m, [0.82, 0.82, 0.61, 0.13, 0.90, 0.75, 0.75, 0.26, 0.30, 0.26, 0.19, 0.34, 0.43, 0.89],) 




