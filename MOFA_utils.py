import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
from sklearn.preprocessing import StandardScaler

Time10= ['10']*21 + ['24']*21
Time6= ['6']*21 + ['24']*21
avg_col_name10 = ['feature','view','T10_T10C0','T10_T10C001','T10_T10C01','T10_T10C1','T10_T10C10','T10_T10C100','T10_T10C1000',\
                    'T24_T24C0','T24_T24C001','T24_T24C01','T24_T24C1','T24_T24C10','T24_T24C100','T24_T24C1000']
avg_col_name6 = ['feature','view','T6_T6C0','T6_T6C001','T6_T6C01','T6_T6C1',\
                'T6_T6C10','T6_T6C100','T6_T6C1000','T24_T24C0','T24_T24C001',\
                'T24_T24C01','T24_T24C1','T24_T24C10','T24_T24C100','T24_T24C1000']

# Preprocess the omics data for Multi-Omics Factor Analysis (MOFA)
def pre_process4MOFA(omics, id_colname, omics_name):

    # Drop rows with fewer than 15 non-null values
    omics.dropna(thresh=15, inplace=True)

    # Transpose the data, set column names, and convert to float data type
    omicsT = omics.T
    omicsT.columns = omicsT.iloc[0]
    omicsT.drop([id_colname], inplace= True)
    omicsT.reset_index(drop = True, inplace = True)
    omicsT = omicsT.rename_axis(None, axis=1)
    omicsT = omicsT.astype('float')

    # Assign time points and average column names based on the omics name
    if omics_name == 'Transcriptome':
        Time = Time6
        avg_col_name = avg_col_name6
    elif omics_name == 'Proteome' or omics_name == 'Phosphoproteome':
        Time = Time10
        avg_col_name = avg_col_name10
    else:
        print('Please enter valid omics layer')

    # Assign concentrations, time points, and extract feature list
    Con = (['0']*3+['0.01']*3+['0.1']*3+['1']*3+['10']*3+['100']*3+['1000']*3)*2
    omicsT['Time'] = Time
    omicsT['Concen'] = Con
    feature_list = list(omicsT.columns)[0:-2]

    # Perform ANOVA for each feature and identify features with p-values < 0.05
    changed_feature = []
    for feature in feature_list:
        model = ols('{} ~ C(Time) + C(Concen) + C(Time):C(Concen)'.format(feature), data=omicsT).fit()
        df_results = sm.stats.anova_lm(model, typ=2)
        if (np.array(df_results.iloc[0:3, 3])<0.05).any():
            changed_feature.append(feature)
    
    # Subset the omics data for the identified features
    de_features = omics[omics[id_colname].isin(changed_feature)]
    de_features.reset_index(drop= True, inplace = True)
    
    # Compute the average for three replicates
    omics_avg = pd.DataFrame(columns = avg_col_name)
    omics_avg['feature'] = de_features[id_colname]
    omics_avg['view'] = omics_name
    for i in range(14):
        omics_avg.iloc[:,i+2] = np.mean(de_features.iloc[:,i*3+1:i*3+4],axis = 1)

    # # Scale the averaged values
    scaler = StandardScaler()
    scaler.fit(omics_avg.iloc[:,2:].T)
    omics_avg_scaled = scaler.transform(omics_avg.iloc[:,2:].T) 
    omics_avg_scaled = pd.concat([omics_avg.iloc[:,:2],pd.DataFrame(omics_avg_scaled).T],axis=1) 
    omics_avg_scaled.columns = avg_col_name
    
    return omics_avg_scaled