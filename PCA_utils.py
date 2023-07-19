import numpy as np
import matplotlib.lines as mlines
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def pcaplt(data,omics, early_time= "10 min",con_loc = "upper right", t_loc = "lower right"):
    # Define a function called pcaplt to perform PCA analysis and plot the results
    
    # Drop rows with missing values from the data
    data = data.dropna(how = 'any')
    # Standardize the data 
    data_scaled = StandardScaler().fit(data.iloc[:,1:].T).transform(data.iloc[:,1:].T)
    # Perform PCA on the scaled data and retain the top 2 principal components
    pca_data = PCA(n_components=2).fit(data_scaled)
    data_pca = pca_data.transform(data_scaled)
    print('Original shape:{}'.format(str(data_scaled.shape)))
    print('Reduced shape:{}'.format(str(data_pca.shape)))
    
    # Get the explained variance ratios of the principal components
    pca_ratio= pca_data.explained_variance_ratio_
    
    # Create a new figure for the plot
    title_font = {'family': 'sans-serif','fontname': 'Helvetica','size': 16}
    axtitle_font = {'family': 'sans-serif','fontname': 'Helvetica','size': 18}
    
    plt.figure(figsize=(9.5,7.4))
    
    colors = cm.Blues(np.linspace(0.1, 1, 7))
    concen = ['0 μM','0.01 μM','0.1 μM','1 μM','10 μM','100 μM','1000 μM']
    
    for i in range(7):
        plt.scatter(data_pca[i*3:i*3+3,0],data_pca[i*3:i*3+3,1],
                    color = colors[i], s= 90)
    # Create a legend for the caffeine concentration
    first_legend = plt.legend(concen, title = 'Caffeine concentration',\
                                title_fontsize=14, loc=con_loc, fontsize="14")


    for i in range(7):
        plt.scatter(data_pca[i*3+21:i*3+24,0],data_pca[i*3+21:i*3+24,1],
                    color = colors[i], marker = "^", s = 90)

    # Create the legend to indicate marker shape
    circle = mlines.Line2D([], [], color='black', marker='o', linestyle='None',\
                                markersize=10, label=early_time)
    triangle = mlines.Line2D([], [], color='black', marker='^', linestyle='None',\
                                markersize=10, label='24 h')
    second_legend = plt.legend(handles=[circle, triangle], loc=t_loc, title="Time",\
                                title_fontsize=14,fontsize="14")

    # Add the first legend manually to the current Axes.
    plt.gca().add_artist(first_legend)
    
    # Set the axis labels and plot title
    plt.xlabel('PCA Component 1 (Variance ratio={:.2f})'.format(pca_ratio[0]),fontdict=axtitle_font)
    plt.ylabel('PCA Component 2 (Variance ratio={:.2f})'.format(pca_ratio[1]),fontdict=axtitle_font)
    plt.title(omics,fontdict=title_font)

    plt.show()