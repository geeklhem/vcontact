
def plot_niches(subset,metadata,title="Niches"):
    for i,(q1,q3) in enumerate(zip(metadata.loc[:,cols].quantile(0.25),metadata.loc[:,cols].quantile(0.75))):
        plt.axhspan(xmin=q1,xmax=q3,ymin=i,ymax=i+1, hatch="/", alpha=0.1)  
    for i,(q1,q3) in enumerate(zip(subset.ix[:,cols].quantile(0.25),subset.ix[:,cols].quantile(0.75))):
        plt.axhspan(xmin=q1,xmax=q3,ymin=i,ymax=i+1, facecolor='g', alpha=0.3)
    for i,(mx,mi) in enumerate(zip(subset.ix[:,cols].max(),subset.ix[:,cols].min())):
        plt.axhspan(xmin=mi,xmax=mi+0.005,ymin=i,ymax=i+1, alpha=0.5,facecolor="y")
        plt.axhspan(xmin=mx,xmax=mx+0.005,ymin=i,ymax=i+1, alpha=0.5,facecolor="y")
    for i,mean in enumerate(subset.ix[:,cols].mean()):
        plt.axhspan(xmin=mean,xmax=mean+0.005,ymin=i,ymax=i+1, facecolor='g', alpha=1)
    for i,mean in enumerate(metadata.loc[:,cols].mean()):
        plt.axhspan(xmin=mean,xmax=mean+0.005,ymin=i,ymax=i+1, facecolor='k', alpha=1)
        
    s_range = plt.Rectangle((0, 0), 1, 1, fc="g", alpha=0.3,)
    ref_q50 = plt.Rectangle((0, 0), 1, 1, hatch="/", alpha=0.1)
    ref_m = plt.Line2D((0,), (1,), 2,color='k')
    s_m = plt.Line2D((0,), (1,), 2,color='g')
    s_mm = plt.Line2D((0,), (1,), 2,color='y')
    
    
    leg = plt.legend([s_range,ref_q50,ref_m,s_m,s_mm],["q50 of subset",
                                       "q50 of reference",
                                       "Mean of reference",
                                   "Mean of subset","Subset min/max"],loc='best',fancybox=True)
    leg.get_frame().set_alpha(0.5)
        
        
    plt.ylim(0,i+1)
    plt.xlim(0,1)
    plt.xticks([0, 1],
           ["min","max"])
    plt.yticks(np.arange(len(cols))+0.5,
           cols)
    plt.title(title)
