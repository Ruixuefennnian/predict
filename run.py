import pandas as pd
import numpy as np 
import re
import predict
'''
axial_tm = pd.read_excel('descriptor.xlsx', sheet_name='axial', index_col=0)
lateral_tm = pd.read_excel('descriptor.xlsx', sheet_name='lateral', index_col=0)
backbone = pd.read_excel('descriptor.xlsx', sheet_name='propensity', index_col=0)

if model == "1":
	GPY_KD = backbone['GPY']  
	GPY_no = backbone['GPY']
	GY_KD = backbone['GPY']
	GY_no = backbone['GPY']
	GXO_KD = backbone['GXO']
	GXO_no = backbone['GXO']
	GXP_KD = backbone['GXO']
	GXP_no = backbone['GXO']
	GX_no = backbone['GXO']   	
	GX_KD = backbone['GXO']
elif model == "2":
	GPY_KD = backbone['GPY-KD']  
	GPY_no = backbone['GPY']
	GY_KD = backbone['GPY-KD']
	GY_no = backbone['GPY']
	GXO_KD = backbone['GXO-HB']
	GXO_no = backbone['GXO']
	GXP_KD = backbone['GXO-HB']
	GXP_no = backbone['GXO']
	GX_no = backbone['GXO']   	
	GX_KD = backbone['GXO-HB'] 
elif model == "3":
	GPY_KD = backbone['GPY']  
	GPY_no = backbone['GPY']
	GY_KD = backbone['GY']
	GY_no = backbone['GY']
	GXO_KD = backbone['GXO']
	GXO_no = backbone['GXO']
	GXP_KD = backbone['GXP']
	GXP_no = backbone['GXP']
	GX_no = backbone['GX']   	
	GX_KD = backbone['GX']
else:
	GPY_KD = backbone['GPY-KD']  
	GPY_no = backbone['GPY']
	GY_KD = backbone['GY-KD']
	GY_no = backbone['GY']
	GXO_KD = backbone['GXO-HB']
	GXO_no = backbone['GXO']
	GXP_KD = backbone['GXP']
	GXP_no = backbone['GXP']
	GX_no = backbone['GX']   	
	GX_KD = backbone['GX']
'''
if __name__ == "__main__":
	# Input paramers 4
	"""
	example:
	seq = ['PPGPKGDPGPKGDPGPPGARGQAGVMGFOGPPGPPG',
			'PKGPPGDKGPPGDPGPPGARGEOGNIGFOGPPGPKG',
			'PPGPKGPKGDPGDPGPPGARGQAGVMGFOGPPGPKG']
	"""
	chain_a = input("please input the chain A sequence: ").strip()
	chain_b = input("please input the chain B sequence: ").strip()
	chain_c = input("please input the chain C sequence: ").strip()

	paramers1 = [chain_a, chain_b, chain_c]

	paramers2 = input("which model select (e.g 1, 2, 3, 4):").strip()
	paramers3 = int(input("Include N terminal (yes: 1 or No: 0):").strip())
	paramers4 = int(input("Include C terminal (yes: 1 or No: 0):").strip())

	obj = predict.PredictXY(paramers1, paramers2)
	axi, lat = obj.wSide(0.3,0.7,0.9)
	Hbond = 66*np.log(len(paramers1[0]))-156.33
	weights = obj.weights() 
	trangle = [v for v in obj.GXY()[0].values()][1:]+[axi, lat, Hbond]
	
	Nterminal = 3 if paramers3 == 1 else 0
	Cterminal = 8 if paramers4 == 1 else 0

	winTm = np.dot(np.array(trangle), np.array(weights).T)+Nterminal+Cterminal
	
	print('Tm: ', '%0.1f'%winTm)

