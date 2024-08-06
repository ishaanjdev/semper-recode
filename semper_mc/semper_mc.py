import pandas as pd
import numpy as np

import seaborn as sns
from sklearn.datasets import load_iris
import matplotlib.pyplot as plt


def calc_mRNA_sim(number_cells, load_cell_dist, params, **kwargs):
    """takes in the number of cells (without replicates) that you want to load mRNA with as well as with the loading
    distribution. Suggests a number of mRNAs to simulate
    
    Parameters
    -----------
    number_cells: int
        Number of cells (wihtout replicates) that you want to load with mRNA

    load_cell_dist: (callable)
        A np.random distribution to sample from. Size in the load distribution is the number of 
        times you sample from distribtion. This will be entered in **kwargs.

    Returns
    -------
    int
        Number of mRNAs to simulate in your class instance.
    
    """
    
    samples=load_cell_dist(*params, **kwargs)

    number_mRNA_sim = int(np.round(number_cells * (np.mean(samples) + np.std(samples))))
    return number_mRNA_sim

class SemperPredict:
    def __init__(self, pTIS, number_mRNA):
        """
        Initializes the MCmRNA object.

        Parameters
        ----------
        pTIS : list
            List of floats containing the probabilities of scanning ribosome 
            initiating translation at TIS 1, TIS 2, TIS 3, and so on.

        number_mRNA
            number of mRNA you wish to simulate for your mRNA Array.

        """
        self.pTIS = pTIS
        self.number_mRNA = number_mRNA
        self.number_ORFs = len(pTIS)
        
        
    def sample_ribosomes(self, rib_dist_function, params, **kwargs):
        """Sample the number of ribosomes from a distribution that will traverse the mRNA
        
        Parameters
        ----------
        distribution_function: np.random distribution function
            np.random distribution function (i.e. np.random.normal())

        **kwargs: 
            Keyword arguments required for to sample from the particular np.random distribution

        Returns
        -------
        np.array
            1xnumber_mRNA array of ribosome numbers that will traverse each mRNA in 1:number_mRNA
        
        """
        # sample the number of ribosomes that will cross each mRNA that you will simulate.
        random_ribosomes = rib_dist_function(*params, **kwargs, size = self.number_mRNA)
    
        # Round the samples to the nearest integer to discretize and ensure they are not less than 0.
        # If they are 0, set them to 0.
        # This can be improved. Resample negative values until positive for instance.
        random_ribosomes = np.maximum(np.round(random_ribosomes), 0).astype(int)
        
        self.random_ribosomes = random_ribosomes

        return 
    
    def make_mRNA_array(self):
        """Create an array of the results of each MCmRNA simulation

        Parameters:
        mc_sim: callable
            The monte carlo simulation method you want to use to simulate each mRNA.
        """
     
        mRNA_array = np.zeros((self.number_mRNA, self.number_ORFs))
    
        for i, r in enumerate(self.random_ribosomes):
            mRNA_array[i, :] = self.simulate_MCmRNA(r)
        
        self.mRNA_array = mRNA_array
        return 
    
    
    def simulate_MCmRNA(self, trav_rib):
        """ribosomes_i will traverse an mRNA and produce either ORF 1 based on a bernoulli success probability
        
        Parameters
        ----------
        trav_rib: int
            Number of ribosomes to simulate that will traverse the mRNA
        
        Returns
        -------

        np.array
            2x1 array of with amount of protein 1 and protein 2 made.
        """
        
        proteins = np.zeros([1, self.number_ORFs])
        
        # for every r ribosome, simulate if the first, second, third, none, etc. ORF
        # is successfully translated.
        for r in range(0, trav_rib):
            
            # for each value in pTIS (indicative of TIS in front of ORF from 5'-->3')
            # sample from a bernoulli distribution using the success probability pTIS_prob
            for i, pTIS_prob in enumerate(self.pTIS):
                
                if np.random.uniform(low=0, high=1) <= pTIS_prob:
                    proteins[0,i] += 1
                    break # if criteria met, simulate next ribosome

        return proteins
    
    
    def load_cells(self, number_cells, number_replicates, describe, load_cell_dist, params, **kwargs):
        """Given the number of cells needed to be predicted, the function will 
        take in an mRNA_array and fill the cells with j mRNA where j is a number sampled from a normal distribution
        
        Parameters
        ----------
        number_cells: int
            The number of cells you want to load with mRNA. You can increase the number of replicates
            to resample the simulated mRNAs and increase this number without having to simulate new mRNA

        number_replicates: int
            Number of replicates to produce by resampling the mRNA array

        describe: str
            String describing the type of gene circuit being tested.

        load_cell_dist: (callable)
            Np.random distribution function for sampling number of mRNA results to go into each cell

        Returns
        -------
        pd.DataFrame
            Realizations of predicted protein quantities for single cells.
        """

        cell_array = []

        for r in range(0, number_replicates):
            mRNA_array_rep = np.copy(self.mRNA_array)
            np.random.shuffle(mRNA_array_rep)

            # sample the number of mRNA results to load into each cell
            random_mRNA = load_cell_dist(*params, **kwargs, size = number_cells)
            
            # round the number of mRNA to the nearest integer and make sure the value is positive, 
            # if negative set the value to 0
            random_mRNA = np.maximum(np.round(random_mRNA), 0).astype(int)
            
            # load m mRNAs into a cell, and then delete them from the shuffled array of mRNA simulation results.
            for i, m in enumerate(random_mRNA):
                cell_array.append(mRNA_array_rep[0:m,:].sum(axis=0))
                mRNA_array_rep = np.delete(mRNA_array_rep, slice(0,m), axis=0)
                
                
        output_df = pd.DataFrame(cell_array)
        output_df["description"] = describe
        
        self.simulated_cells = output_df
        return 
        
                