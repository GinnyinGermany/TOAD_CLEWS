import numpy as np
from netCDF4 import Dataset
import warnings
import os

os.environ['HDF5_USE_FILE_LOCKING']='FALSE'


class global_functions():
    def Rain(datasets, rain_fact):
        """
        Computation of rainfall values; all datasets are necessary and a potential rain factor
        """
        rain = []
        for data in datasets:
            net_data = Dataset(data)
            #get mean value for rainfall
            rain_dataset = np.multiply(rain_fact, net_data.variables["prec"][:])
            if len(rain) == 0:
                rain = rain_dataset
            else:
                rain = np.add(rain, rain_dataset)

        rain = np.array(rain)
        return rain


    def Rain_crit(data_crit, data_crit_std, adapt_fact):
        """
        Computation of critical rainfall values; all datasets are necessary and a potential rain factor
        """
        rain_crit = []
        for i in range(0, len(data_crit)):
            net_data = Dataset(data_crit[i])
            net_data_std = Dataset(data_crit_std[i])

            #get mean value for rainfall
            if type(adapt_fact) is np.ndarray:
                rain_dataset = np.subtract(net_data.variables["prec"][:], np.multiply(net_data_std.variables["prec"][:], adapt_fact[:]))
            else:
                rain_dataset = np.subtract(net_data.variables["prec"][:], np.multiply(net_data_std.variables["prec"][:], adapt_fact))

            rain_dataset[rain_dataset < 0.0] = 0.0
            if len(rain_crit) == 0:
                rain_crit = rain_dataset
            else:
                rain_crit = np.add(rain_crit, rain_dataset)

        rain_crit = np.array(rain_crit)

        global_limit_rain = 0.0
        rain_crit[rain_crit < global_limit_rain] = global_limit_rain
        return rain_crit


    def Mcwd(datasets, rain_fact):
        """
        Computation of MCWD; all datasets are necessary and a potential rain factor
        """
        mcwd_array = []
        for data in datasets:
            net_data = Dataset(data)
            #get mean value for rainfall
            rain = np.multiply(rain_fact, net_data.variables["prec"][:])
            evap = net_data.variables["evap"][:]
            #computation of MCWD
            diff = np.subtract(evap, rain)
            mcwd_array.append(diff)
        mcwd_array = np.array(mcwd_array)

        mcwd = []
        for i in range(0, len(mcwd_array[0])): #necessary for assessing an mcwd value for each cell
            mcwd_probe = []
            for j in mcwd_array:
                mcwd_probe.append(j[i])
            mcwd_probe = np.array(mcwd_probe)
            mcwd_sign = np.sign(mcwd_probe)

            mcwd_real = []
            diff = 0.
            for j in range(0, len(mcwd_probe)):
                if mcwd_sign[j] < 0.:
                    diff = 0.
                else:
                    diff += mcwd_probe[j]
                mcwd_real.append(diff)
            mcwd_real = np.array(mcwd_real)

            mcwd.append(np.amax(mcwd_real))
        mcwd = np.array(mcwd)
        return mcwd


    def Mcwd_crit(data_crit, data_crit_std, adapt_fact):
        """
        Computation of critical MCWD values; all datasets are necessary and a potential rain factor
        """
        mcwd_array = []
        for i in range(0, len(data_crit)):
            net_data = Dataset(data_crit[i])
            net_data_std = Dataset(data_crit_std[i])
            #get mean value for rainfall
            if type(adapt_fact) is np.ndarray:
                rain = np.subtract(net_data.variables["prec"][:], np.multiply(net_data_std.variables["prec"][:], adapt_fact[:]))
            else:
                rain = np.subtract(net_data.variables["prec"][:], np.multiply(net_data_std.variables["prec"][:], adapt_fact))
            #values smaller zero do not make sense
            rain[rain<0.0] = 0.0

            evap = net_data.variables["evap"][:]

            #computation of MCWD
            diff = np.subtract(evap, rain)
            mcwd_array.append(diff)
        mcwd_array = np.array(mcwd_array)


        mcwd_crit = []
        for i in range(0, len(mcwd_array[0])): #necessary for assessing an mcwd_crit value for each cell
            mcwd_probe = []
            for j in mcwd_array:
                mcwd_probe.append(j[i])
            mcwd_probe = np.array(mcwd_probe)
            mcwd_sign = np.sign(mcwd_probe)

            mcwd_real = []
            diff = 0.
            for j in range(0, len(mcwd_probe)):
                if mcwd_sign[j] < 0.:
                    diff = 0.
                else:
                    diff += mcwd_probe[j]
                mcwd_real.append(diff)
            mcwd_real = np.array(mcwd_real)

            mcwd_crit.append(np.amax(mcwd_real))
        mcwd_crit = np.array(mcwd_crit)

        #here also absolute limits can be set
        global_limit_mcwd = 0.0
        mcwd_crit[mcwd_crit < global_limit_mcwd] = global_limit_mcwd
        return mcwd_crit        
      
    def Amazon_CUSPc(a, b, rain_mean, rain_critical, rain_current, mcwd_mean, mcwd_critical, mcwd_current):
        if mcwd_critical < mcwd_mean:
            print("""Error: The mean seasonality value is below the average.
                This does not make sense since on average the Amazon rainforest would be tipped - this is not what is observed...""")
            die

        c_rain = np.sqrt((4*np.abs(b)**3) / (27*np.abs(a)))/(rain_critical - rain_mean)*(rain_current - rain_mean)

        if mcwd_critical == 0.0 and mcwd_mean == 0.0:
            c_mcwd = 0.0
        else:
            c_mcwd = np.sqrt((4*np.abs(b)**3) / (27*np.abs(a)))/(mcwd_critical - mcwd_mean)*(mcwd_current - mcwd_mean)


        #case distinction to find "correct" criitcal value c
        if 0 < c_rain < np.sqrt(4/27) and 0 < c_mcwd < np.sqrt(4/27):
            c = np.amax([c_rain, c_mcwd]) + (np.sqrt(4/27) - np.amax([c_rain, c_mcwd]))/(np.sqrt(4/27))*np.amin([c_rain, c_mcwd])
        else:
            c = np.amax([c_rain, c_mcwd])

        return c


    def Amazon_cpl(a, b, rain_mean, rain_critical, rain_current, delta_rain, mcwd_mean, mcwd_critical, mcwd_current, delta_mcwd):
        c_rain = np.sqrt((4*np.abs(b)**3) / (27*np.abs(a)))/(rain_critical - rain_mean)*(rain_current - rain_mean)
        cpl_rain = np.sqrt((4*np.abs(b)**3) / (27*np.abs(a)))/(rain_mean - rain_critical)*(1/2)*(-1)*delta_rain

        if mcwd_critical == 0.0 and mcwd_mean == 0.0:
            c_mcwd = 0.0
            cpl_mcwd = 0.0 
        else:
            c_mcwd = np.sqrt((4*np.abs(b)**3) / (27*np.abs(a)))/(mcwd_critical - mcwd_mean)*(mcwd_current - mcwd_mean)
            cpl_mcwd = np.sqrt((4*np.abs(b)**3) / (27*np.abs(a)))/(mcwd_mean - mcwd_critical)*(1/2)*(-1)*delta_mcwd

        if 0 < c_rain + cpl_rain < np.sqrt(4/27) and 0 < c_mcwd + cpl_mcwd < np.sqrt(4/27):
            index = np.argmax([c_rain, c_mcwd])
            if index  == 0:
                cpl = cpl_rain + (np.sqrt(4/27) - cpl_rain)/(np.sqrt(4/27))*cpl_mcwd
            elif index == 1:
                cpl = cpl_mcwd + (np.sqrt(4/27) - cpl_mcwd)/(np.sqrt(4/27))*cpl_rain
            else:
                print("Wrong index!!")
                die
        else:
            index = np.argmax([c_rain, c_mcwd])
            if index == 0:
                cpl = cpl_rain
            elif index == 1:
                cpl = cpl_mcwd
            else:
                print("Wrong index!!")
                die

        if cpl < 0.0:
            print("Coupling strengths below 0.0 are not allowed")
            die
    
        return cpl



    def Rain_moisture_delta_only(moist_rec_val):
        rain_moist = - np.sum(moist_rec_val)
        return rain_moist


    def Mcwd_moisture(datasets, rain_fact, id_receiving_cell, moist_rec_val):
        mcwd_array = []
        for i in range(0, len(datasets)):
            net_data = Dataset(datasets[i])
            #get mean value for rainfall
            rain = np.multiply(rain_fact, net_data.variables["prec"][id_receiving_cell])
            evap = net_data.variables["evap"][id_receiving_cell]

            #computation of MCWD
            diff = np.subtract(evap, np.subtract(rain, moist_rec_val[i]))
            mcwd_array.append(diff)
        mcwd_array = np.array(mcwd_array)
        mcwd_sign = np.sign(mcwd_array)


        mcwd_real = []
        diff = 0.
        for j in range(0, len(mcwd_array)):
            if mcwd_sign[j] < 0.:
                diff = 0.
            else:
                diff += mcwd_array[j]
            mcwd_real.append(diff)
        mcwd_real = np.array(mcwd_real)
        

        mcwd = np.amax(mcwd_real)
        return mcwd