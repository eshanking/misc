from lotka_volterra import LV
#%% Kaznatcheev et al

i_s = {'p':0.8,
       'r':0.2}

# DMSO + CAF
payoff = {'a':2.6,
            'b':3.6,
            'c':3.1,
            'd':3.0}

lv = LV(payoff=payoff,init_state=i_s,dt = 0.001)
lv.simulate()
lv.plot_timecourse(title = 'DMSO + CAF')

i_s = {'p':0.8,
       'r':0.2}

# DMSO
payoff = {'a':2.5,
            'b':2.4,
            'c':4.0,
            'd':2.7}
lv = LV(payoff=payoff,init_state=i_s,dt = 0.001)
lv.simulate()
lv.plot_timecourse(title = 'DMSO')

i_s = {'p':0.8,
       'r':0.2}

# Alectinib + CAF
payoff = {'a':0.5,
            'b':-0.4,
            'c':3.8,
            'd':2.4}
lv = LV(payoff=payoff,init_state=i_s,dt = 0.001)
lv.simulate()
lv.plot_timecourse(title = 'Alectinib + CAF')

i_s = {'p':0.8,
       'r':0.2}

# Alectinib

payoff = {'a':-1.0,
            'b':-1.3,
            'c':4.3,
            'd':2.3}
lv = LV(payoff=payoff,init_state=i_s,dt = 0.001)
lv.simulate()
lv.plot_timecourse(title = 'Alectinib')