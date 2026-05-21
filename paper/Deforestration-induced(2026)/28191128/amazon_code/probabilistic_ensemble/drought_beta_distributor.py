import numpy as np
import matplotlib.pyplot as plt

#1x1°lon-lat grid
size = 416

#upper and lower limit of standard deviations that should be computed
upper_limit = 1.25
lower_limit = 0.75

beta_randomizer = 2.53

for i in range(0, 100):
    print(i)
    beta_sample = (upper_limit - lower_limit)*np.random.beta(beta_randomizer, beta_randomizer, size) + lower_limit
    np.savetxt("drought_start_sample/{}.txt".format(str(i).zfill(3)), np.array(beta_sample))
    
    #get image
    plt.hist(np.array(beta_sample), 50)
    plt.savefig("drought_start_sample/{}.png".format(str(i).zfill(3)))
    plt.close()
    plt.clf()

print("Finish")