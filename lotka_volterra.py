import numpy as np
import matplotlib.pyplot as plt

class LV:

    def __init__(self,init_state=None,payoff=None,n_timestep=1000,dt = 0.1):

        self.n_timestep = n_timestep

        if init_state is None:
            init_state = {'p':0.5,'r':0.5}
        
        self.state = init_state
        self.init_state = init_state

        if payoff is None:
            payoff = {'a':0.5,
                      'b':0,
                      'c':0,
                      'd':0.5}

        self.payoff = payoff

        self.timecourse = {'p':np.zeros(n_timestep+1),
                           'r':np.zeros(n_timestep+1)}

        self.timecourse['p'][0] = self.init_state['p']
        self.timecourse['r'][0] = self.init_state['r']

        self.dt = dt

    
    def iterate(self):

        p = self.state['p']
        r = self.state['r']

        dp,dr = self.growth_rate(p,r,self.payoff,which='both')

        p_t = p + dp*self.dt
        r_t = r + dr*self.dt

        if r_t < 0:
            r_t = 0
        
        if p_t < 0:
            p_t = 0

        tot = p_t + r_t

        p_t = p_t/tot
        r_t = r_t/tot

        self.state['p'] = p_t
        self.state['r'] = r_t

    def simulate(self):

        for t in range(self.n_timestep):
            self.iterate()
            self.timecourse['p'][t+1] = self.state['p']
            self.timecourse['r'][t+1] = self.state['r']
    

    def growth_rate(self,p,r,payoff,which='both'):
        """growth rate equations

        Args:
            p (float): proportion parental
            r (float): proportion resistant
            payoff (dict): payoff matrix
            which (str): determines output type. 'resistant' returns resistant growth rate. 
                'parental' returns parental growth rate. 'both' returns both in the form of (parental,resistant) tuple

        Raises:
            ValueError: _description_

        Returns:
            _type_: _description_
        """
        a = payoff['a']
        b = payoff['b']
        c = payoff['c']
        d = payoff['d']

        if which == 'parental':
            return a*p + b*r
        elif which == 'resistant':
            return c*p + d*r
        elif which == 'both':
            return a*p + b*r, c*p + d*r
        else:
            raise ValueError('Unknown argument given for which.')
    
    def plot_timecourse(self,title=None):

        fig,ax = plt.subplots()
        
        for key in self.timecourse.keys():
            ax.plot(self.timecourse[key],label=key)
        
        ax.set_title(title)

        ax.set_ylim(-0.1,1.1)
        
        ax.legend(frameon=False)
