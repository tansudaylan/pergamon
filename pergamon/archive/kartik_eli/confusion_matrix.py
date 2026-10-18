from tdpy.verbosity import print
import numpy as np

def confusion_matrix(actual, predicted, positive=1, negative=0, prt=False):
    '''
    calculates a confusion matrix based on actual labels and predicted labels
    for binary classification problem

    actual, predicted are 1D numpy arrays or lists of the same length

    returns: 2x2 numpy array representing the confusion matrix for the data
    '''

    result = np.zeros((2,2))
    num_vals = len(actual)

    for i in range(num_vals):
        a = actual[i]
        p = predicted[i]
        if a == negative:           # Actual Negative
            if p == negative:       # True Negative
                result[1][1] += 1
            elif p == positive:     # False Positive
                result[0][1] += 1        
        else:                       # Actual Positive
            if p == positive:       # True Positive
                result[0][0] += 1           
            elif p == negative:     # False Negative
                result[1][0] += 1   

    accuracy = (result[0][0] + result[1][1]) / num_vals
    recall = result[0][0] / (result[0][0] + result[1][0])
    precision = result[0][0] / (result[0][0] + result[0][1])
    F1 = 2 * precision * recall / (precision + recall)
    if prt:
        print("Confusion Matrix: \n", result)
        print("Accuracy: ", round(accuracy, 2))
        print("Recall: ", round(recall, 2))
        print("Precision: ", round(precision, 2))
        print("F1: ", round(F1, 2))

    return result, accuracy, precision, recall, F1

