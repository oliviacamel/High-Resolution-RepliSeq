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

def find_plateau_peaks(ArgMaxArr):
    """
    Identifies plateau peaks in an index array of argmax of birch model subcluster centers

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

    if start is not None:
        left_lower = start == 0 or ArgMaxArr[start - 1] < ArgMaxArr[start]
        right_lower = ArgMaxArr[-1] < ArgMaxArr[start]
        if left_lower and right_lower:
            plateaus.append([start,len(ArgMaxArr)-1])

    return plateaus

def main(filename, n_clusters=None, threshold=0.5):
    """
    load data, perform Birch clustering, and find IZ peak indices

    Parameters:
    - filename (str): Path to the .npy file containing the input array, assuming the array is a scaled repli-seq array for a single chromosome with shape (n Sphase fractions, genomic bins)
    - n_clusters (int, optional): Number of clusters for Birch clustering.
    - threshold (float, optional): Threshold for Birch clustering.
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
        IZs = find_plateau_peaks(ArgMaxArr)
        print("IZ indices:", IZs)
        output_filename = os.path.splitext(filename)[0] + "_IZIndices.npy"
        np.save(output_filename,np.array(IZs))
        print(f"Plateau indices written to {output_filename}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform Birch clustering and identify plateau peaks.")
    parser.add_argument("filename", type=str, help="Path to the .npy file containing the input array.")
    parser.add_argument("--n_clusters", type=int, default=None, help="Number of clusters for Birch clustering.")
    parser.add_argument("--threshold", type=float, default=0.5, help="Threshold for Birch clustering.")

    args = parser.parse_args()

    main(args.filename, n_clusters=args.n_clusters, threshold=args.threshold)
