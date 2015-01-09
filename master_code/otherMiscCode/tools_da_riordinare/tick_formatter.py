import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker 

class MyFormatter(ticker.LogFormatter): 
    def __call__(self, val, pos=None): 
        fx = int(np.floor(np.log(abs(val))/np.log(self._base) +0.5))  
        isDecade = self.is_decade(fx) 
        if not isDecade and self.labelOnlyBase: 
            return '' 
        if (fx%2)==1:
            return '' 
        return '%d'%fx 
