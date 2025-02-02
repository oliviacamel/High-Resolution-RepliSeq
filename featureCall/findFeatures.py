import numpy as np
import argparse
import sys
import os


def cluster_by_birch(array: np.ndarray, n_clusters: int = None, threshold: float = 0.5) -> np.ndarray:
    """
    Birch cluster a repli-seq array of shape (fractions, genomic bins).

    Parameters:
    - array (np.ndarray): Input 2D array of shape (rows, columns)
    - n_clusters (int, optional): Number of clusters. If None, Birch will determine it.
    - threshold (float, optional): The radius threshold for clustering.

    Returns:
    - np.ndarray: Birch model
    """

    from sklearn.cluster import Birch
    # Ensure input is valid
    if not isinstance(array, np.ndarray) or array.ndim != 2:
        raise ValueError("Input must be a 2D NumPy array.")

    # Check for NaN or infinite values
    if not np.isfinite(array).all():
        raise ValueError("Input array contains NaN or infinite values.")

    try:
        birch = Birch(n_clusters=n_clusters, threshold=threshold)
        arrayT = array.T
        model = birch.fit(arrayT)
        return model
    except Exception as e:
        print(f"Error during Birch clustering: {e}")
        return None

'''
def _emergingTime(arr, threshold):
    """
    Find the first S phase fraction where the value in the arr exceeds an arbitrary threshold.

    Parameters:
    - arr (array-like): Input 1D array representing a genomic bin's signal across S phase.
    - threshold (float): An arbitrary threshold for signal.

    Returns:
        int: Index of the first element exceeding the threshold.
    """
    arr = np.array(arr)
    exceed_indices = np.where(arr > threshold)[0]

    if exceed_indices.size > 0:
        return exceed_indices[0]
    else:
        return len(arr) - 1  # If no value exceeds, return last fraction

'''

def find_peaks(ArgMaxArr):
    """
    Identifies plateau peaks in an index array of argmax of birch model subcluster centers
    
    Parameters:
    - ArgMaxArr (array-like): Input 1D array representing the S phase fraction of max signal for a given genomic bin.
    
    Returns:
    - list: list of indices for IZs
    """
    ArgMaxArr = np.array(ArgMaxArr)
    plateaus = []
    start = None

    for i in range(1, len(ArgMaxArr)):
        if ArgMaxArr[i] == ArgMaxArr[i - 1]:
            if start is None:
                start = i - 1
        else:
            if start is not None:
                end = i - 1
                left_lower = start == 0 or ArgMaxArr[start - 1] < ArgMaxArr[start]
                right_lower = i == len(ArgMaxArr) or ArgMaxArr[i] < ArgMaxArr[start]

                if left_lower and right_lower:
                    plateaus.append([start,end])

                start = None
        #allow single-bin peaks
        if (i == 0 or ArgMaxArr[i - 1] < ArgMaxArr[i]) and (i == len(ArgMaxArr) - 1 or ArgMaxArr[i + 1] < ArgMaxArr[i]):
            plateaus.append([i, i])
    if start is not None:
        left_lower = start == 0 or ArgMaxArr[start - 1] < ArgMaxArr[start]
        right_lower = ArgMaxArr[-1] < ArgMaxArr[start]
        if left_lower and right_lower:
            plateaus.append([start,len(ArgMaxArr)-1])

    return plateaus

def find_rightward_slopes(ArgMaxArr):
    """
    Identifies rightward slopes in a 1D NumPy array.
    Parameters:
    - ArgMaxArr (array-like): Input 1D array representing the S phase fraction of max signal for a given genomic bin.

    Returns:
    - list: list of indices for rightward TTRs
    """
    ArgMaxArr = np.array(ArgMaxArr)  
    slopes = []
    n = len(ArgMaxArr)
    
    i = 0
    while i < n - 1:
        while i < n - 1 and ArgMaxArr[i] == ArgMaxArr[i + 1]:  
            i += 1

        left_slope_start = i  

        slope_start = i + 1
        while slope_start < n - 1:
            if ArgMaxArr[slope_start] < ArgMaxArr[slope_start - 1]:  
                slope_start += 1
            #allowing breakages
            elif ArgMaxArr[slope_start] == ArgMaxArr[slope_start - 1]:  
                temp = slope_start
                while temp < n - 1 and ArgMaxArr[temp] == ArgMaxArr[temp + 1]:
                    temp += 1
                if temp < n - 1 and ArgMaxArr[temp + 1] < ArgMaxArr[temp]:  
                    slope_end = slope_start - 1
                    if slope_end > left_slope_start:
                        slopes.append([left_slope_start, slope_end])
                    left_slope_start = temp 
                    slope_start = temp + 1
                else:
                    break  
            else:
                break  

        slope_end = slope_start  -1
        if slope_end > left_slope_start:
            slopes.append([left_slope_start, slope_end])

        i = slope_end + 1  

    return slopes


def main(filename, n_clusters=None, threshold=0.5,feature_type = 'IZ'):
    """
    load data, perform Birch clustering, and find feature indices

    Parameters:
    - filename (str): Path to the .npy file containing the input array, assuming the array is a scaled repli-seq array for a single chromosome with shape (n Sphase fractions, genomic bins)
    - n_clusters (int, optional): Number of clusters for Birch clustering.
    - threshold (float, optional): Threshold for Birch clustering.
    - feature_type (str, optional): Type of feature to call. Options are 'IZ' or 'rightwardTTR'.
    """
    try:
        data = np.load(filename)
        data = np.nan_to_num(data).copy()

        # Perform Birch clustering
        model = cluster_by_birch(data, n_clusters=n_clusters, threshold=threshold)
        if model is None:
            raise RuntimeError("Birch clustering failed.")

        # Extract argmax indices from subcluster centers
        ArgMaxArr = np.array([np.argmax(model.subcluster_centers_, axis=1)[i] * -1 for i in model.labels_[:]])

        # Find plateau peaks
        if feature_type == 'IZ':
            features = find_peaks(ArgMaxArr)        
        elif feature_type == 'rightwardTTR':
            features = find_rightward_slopes(ArgMaxArr)
        print(f"{feature_type} indices:", features)
        output_filename = os.path.splitext(filename)[0] + f"_{feature_type}Indices.npy"
        np.save(output_filename,np.array(features))
        print(f"Plateau indices written to {output_filename}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform Birch clustering and identify replication features in a scaled and normed Repli-seq array.")
    parser.add_argument("filename", type=str, help="Path to the .npy file containing the input array.")
    parser.add_argument("--n_clusters", type=int, default=None, help="Number of clusters for Birch clustering.")
    parser.add_argument("--threshold", type=float, default=0.5, help="Threshold for Birch clustering.")
    parser.add_argument("--feature_type", type=str, default='IZ', help="Features type to call, options are IZ, rightwardTTR.")

    args = parser.parse_args()
    if args.feature_type not in ['IZ', 'rightwardTTR']:
        print(f"Error: Invalid feature_type '{args.feature_type}'. Options are 'IZ' or 'rightwardTTR'.", file=sys.stderr)
        sys.exit(1)
main(args.filename, n_clusters=args.n_clusters, threshold=args.threshold, feature_type=args.feature_type)