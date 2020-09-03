from __future__ import division, print_function


import os
import numpy as np
import pandas as pd
import datetime
import itertools
from collections import defaultdict


from src import calibration


class Hymod:
    def __init__(self):
        return

    # Get random parameter set
    @staticmethod
    def get_random_params():

        param_bounds = defaultdict(
            cmax = defaultdict(lower = 1.0, upper = 100),
            bexp = defaultdict(lower = 0.0, upper = 2.0),
            alpha = defaultdict(lower = 0.2, upper = 0.99),
            ks = defaultdict(lower = 0.01, upper = 0.5),
            kq = defaultdict(lower = 0.5, upper = 1.2)
        )

        out = defaultdict()
        minV,maxV,p = None,None,None

        for k in param_bounds.keys():
            minV = param_bounds[k]["lower"]
            maxV = param_bounds[k]["upper"]
            p = np.random.uniform(minV,maxV)
            out[k] = p

        return out

    @staticmethod
    def _power(X,Y):
        X=abs(X) # Needed to capture invalid overflow with netgative values
        return X**Y


    @staticmethod
    def _excess(x_loss,cmax,bexp,Pval,PETval):
        # this function calculates excess precipitation and evaporation
        xn_prev = x_loss
        ct_prev = cmax * (1 - Hymod._power((1 - ((bexp + 1) * (xn_prev) / cmax)), (1 / (bexp + 1))))
        # Calculate Effective rainfall 1
        ER1 = max((Pval - cmax + ct_prev), 0.0)
        Pval = Pval - ER1
        dummy = min(((ct_prev + Pval) / cmax), 1)
        xn = (cmax / (bexp + 1)) * (1 - Hymod._power((1 - dummy), (bexp + 1)))

        # Calculate Effective rainfall 2
        ER2 = max(Pval - (xn - xn_prev), 0)

        # Alternative approach
        evap = (1 - (((cmax / (bexp + 1)) - xn) / (cmax / (bexp + 1)))) * PETval  # actual ET is linearly related to the soil moisture state
        xn = max(xn - evap, 0)  # update state

        return ER1,ER2,xn


    @staticmethod
    def linearReservoir(x_slow,inflow,Rs):
        # Linear reservoir
        x_slow = (1 - Rs) * x_slow + (1 - Rs) * inflow
        outflow = (Rs / (1 - Rs)) * x_slow
        return x_slow,outflow

    @staticmethod
    def simulate(forcings, precipCol="precip", petCol="pet", cmax=None,bexp=None,alpha=None,ks=None,kq=None,initFlow=True):
        """
        Implementation of the Hymod lumped hydrologic model
        See https://www.proc-iahs.net/368/180/2015/piahs-368-180-2015.pdf for a scientific paper.
        Args: precip (pandas.DataFrame): 1-column dataframe with time series of daily precipitation values
              pet (pandas.DataFrame):  1-column dataframe with time series of daily potential evapotranspiration values
        Kwargs: cmax (float): cmax parameter
                bexp (float): bexp parameter
                alpha (float): alpha parameter
                Ks (float): Ks parameter
                Kq (float): Kq parameter
        Returns: outDf (pandas.DataFrame): resulting discharge from the model
        """

        #TODO: add warnings about n columns
        # if len(pcolumns) > 1:

        p = forcings[precipCol].values
        e = forcings[petCol].values

        lt_to_m = 0.001

        # HYMOD PROGRAM IS SIMPLE RAINFALL RUNOFF MODEL
        x_loss = 0.0
        # Initialize slow tank state
        # value of 0 init flow works ok if calibration data starts with low discharge
        x_slow = 2.3503 / (ks * 22.5) if initFlow else 0
        # Initialize state(s) of quick tank(s)
        x_quick = np.zeros(3)
        t = 0
        outflow = np.zeros_like(p)
        output = np.zeros_like(p)
        # START PROGRAMMING LOOP WITH DETERMINING RAINFALL - RUNOFF AMOUNTS

        while t <= len(p)-1:
            Pval = p[t]
            PETval = e[t]
            # Compute excess precipitation and evaporation
            ER1, ER2, x_loss = Hymod._excess(x_loss, cmax, bexp, Pval, PETval)
            # Calculate total effective rainfall
            ET = ER1 + ER2
            #  Now partition ER between quick and slow flow reservoirs
            UQ = alpha * ET
            US = (1 - alpha) * ET
            # Route slow flow component with single linear reservoir
            x_slow, QS = Hymod.linearReservoir(x_slow, US, ks)
            # Route quick flow component with linear reservoirs
            inflow = UQ

            for i in range(3):
                # Linear reservoir
                x_quick[i], outflow = Hymod.linearReservoir(x_quick[i], inflow, kq)
                inflow = outflow

            # Compute total flow for timestep
            output[t] = ((QS + outflow)/lt_to_m)
            t = t+1

        return output


    @staticmethod
    def hargreaves(forcings, tminCol="tmin",tmaxCol="tmax",dtCol="datetime"):
        """
        accepts panda Series
        returns panda Series
        """
        Gsc = 367
        lhov = 2.257

        dts = forcings[dtCol]
        tmin = forcings[tminCol]
        tmax = forcings[tmaxCol]
        n = len(tmax)
        doy = [x.timetuple().tm_yday for x in dts]

        tavg = pd.concat([tmin,tmax],axis=1).mean(axis=1).rename('tavg')

        eto = np.zeros(n)

        for i,t in enumerate(doy):
            b = 2 * np.pi * (t/365)
            Rav = 1.00011 + 0.034221*np.cos(b) + 0.00128*np.sin(b) + 0.000719*np.cos(2*b) + 0.000077*np.sin(2*b)
            Ho = ((Gsc * Rav) * 86400)/1e6

            eto[i] = (0.0023 * Ho * (tmax[i]-tmin[i])**0.5 * (tavg[i]+17.8))

        return eto

    @staticmethod
    def calibrate(forcings,obs,paramSpace, nsamples, precipCol="precip", petCol="pet", saveResults=False):
        samples = dict()

        keyList = list(paramSpace.keys())

        for k in keyList:
            p = paramSpace[k]
            minV = p["lower"]
            maxV = p["upper"]
            samples[k] = np.random.uniform(minV,maxV,size=nsamples)

        losses = np.zeros(nsamples)

        print(f"Running {nsamples} iterations...")
        for i in range(nsamples):
            pars = {k:samples[k][i] for i,k in enumerate(keyList)}
            q = Hymod.simulate(forcings, precipCol=precipCol, petCol=petCol, **pars)
            losses[i] = nse(q,obs)

        loss,idx = losses.max(),losses.argmax()
        finalPars = {k:samples[k][idx] for k in keyList}
        finalQ = Hymod.simulate(forcings, precipCol=precipCol, petCol=petCol, **finalPars)
        finalLoss = nse(finalQ,obs)

        return finalQ, finalPars, loss


def nse(sim,obs):
    numerator = np.nansum((obs-sim)**2)
    denominator = np.nansum((obs-np.nanmean(obs))**2)
    return 1 - (numerator/denominator)
