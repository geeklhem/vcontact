""" Feature selection functions """
import numpy as np 
def gram_smith(matrix):
    """ Do the gram smith orthogonalisation of a matrix """
    output = [matrix[:,0]]
    for i in range(1,matrix.shape[1]):
        substract = [(np.dot(matrix[:,i],u) / float(np.dot(u,u))) * u for u in output]
        print substract
        output.append(matrix[:,i] - sum(substract))
    return  np.transpose(np.array(output))


def gs_ranking(features,target,threshold):
    """Use the gram-smith algorithm to perform a feature ranking
    INPUT : 
    - features (array) data (one observation by line) 
    - target (array) predicted variable (one observation by line) """

    ranking = np.zeros(features.shape[1])
    rscore = np.zeros(features.shape[1])
    i = 0
    while risk < threshold:
        i += 1 
        # Compute the score for each feature
        scores = pow(np.dot(target,features)/((norm(target)*norm(features,axis=0))) ,2)

        # Get the best ranking one 
        ranking[i] = np.argmax(scores)
        rscore[i] = scores[ranking[i]]

        # Project the features on the null-space of the best one
        #        features = 


        # Compute the risk 
        risk = risk + (cumulative(scores[ranking[i]]) * (1-risk))

def cumulative(x,v):
    if v%2: #if v is odd
        if v == 3:
            gamma = 1
        else:
            gamma = 0  
        return np.sqrt(x) * gamma
    else : # if v is even
        if v== 2:
            phi = 0
        elif v==4:
            phi = 1
        else:
            phi = 1 + sum([(2**k *  np.math.fact(k))/np.fact(np.math.factorial(np.math.factorial((2*k+1)))) * (1-x)**k])
        return 2/np.pi * (np.arcsin(np.sqrt(x))+np.sqrt(x*(1-x))*phi)














