# -*- coding: utf-8 -*-
"""
Created on Sun Jan  3 15:45:34 2021

@author: batuhan
"""

import pandas as pd
import numpy as np
def read_excel():
    excel_data = pd.read_excel('beam_moment.xlsx') 
    
    
    n_ex=len(excel_data)
    list_beam = np.zeros((n_ex,4),dtype=float)
    
    beam_list = excel_data['Frame'].tolist()
    mi_list = excel_data['Mi(kNm)'].tolist()
    mspan_list = excel_data['Mspan(kNm)'].tolist()
    mj_list = excel_data['Mj(kNm)'].tolist()
    
    for i in range(0,n_ex,1):
        list_beam[i][0]=beam_list[i] # Label of beams
        list_beam[i][1]=mi_list[i] # List of left support moment 
        list_beam[i][2]=mspan_list[i] # Span moment list
        list_beam[i][3]=mj_list[i]   # List of right support moment 
        
    return n_ex,list_beam
    
    
    
    
